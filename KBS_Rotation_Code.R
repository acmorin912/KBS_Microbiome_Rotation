# =============================================================================
#  KBS LTER MULTI-KINGDOM MICROBIOME, 2018-2020
#  Combined analysis code for Morin et al.
#
#  The project's fifteen working scripts placed end to end in one file, in the
#  order they need to run. The analysis is unchanged. Two kinds of edit were
#  made, and nothing else:
#
#    1  the repeated library() calls and colour definitions were deleted from
#       each script, because both are now handled once at the top of this file
#    2  one bug fix in section 5, documented in full where it appears
#
#  HOW TO RUN
#      # restart R first, so nothing is left over from a previous session
#      source("KBS_combined_analysis.R", echo = TRUE)
#  Use source() rather than pasting blocks by hand. That is what guarantees the
#  order below is respected, and most of the confusing errors in this kind of
#  script come from running things out of order in a session that already has
#  half-built objects in it.
#
#  ORDER
#    - section  1 builds every phyloseq object the rest of the file uses
#    - section  5 builds seq_bac and seq_fungi, which section 6 needs, which is
#                 why Figure 1 comes before the PERMANOVA
#    - sections 5, 12, 13 and 14 write the CSVs that section 15 reads
#
#  WORKING DIRECTORY
#  Section 1 calls setwd() three times: bacterial inputs, fungal inputs, then
#  the output folder. Everything written afterwards goes to the output folder,
#  with one exception flagged in place at section 13.
#
#  CONTENTS
#     1  Pipelines, subsetting and verification          (CombinedRunning.R)
#     2  Rarefaction curves                              (rarecurves.R)
#     3  Community composition tables                    (communitycompositiontables.R)
#     4  Alpha diversity                                 (alpha.R)
#     5  Figure 1 and indicator species analysis         (Fig1.R)
#     6  PERMANOVA                                       (PERMAFINAL.R)
#     7  Pairwise PERMANOVA heatmap                      (r_2.R)
#     8  Distance-based redundancy analysis              (db-rda.R)
#     9  Compartment overlap and cultivar confounding    (Cultivar.R)
#    10  Niche-level clustering                          (NicheLevel.R)
#    11  Management and growth stage PCoA and UpSet      (supplementalpcoa.R)
#    12  Cross-kingdom co-occurrence networks            (HUB.R)
#    13  Random forest                                   (RF_Man.R)
#    14  Differential abundance tables and volcano plots (deseqvolcano.R)
#    15  Four-way analysis overlap                       (ASV_overlap.R)
# =============================================================================


# =============================================================================
#  PACKAGES
#  Every package any section below needs, loaded once.
# =============================================================================
suppressPackageStartupMessages({
  # data handling
  library(dplyr)
  library(tidyverse)
  library(tidyr)
  library(tibble)
  library(data.table)
  library(openxlsx)
  library(yaml)
  
  # microbiome
  library(phyloseq)
  library(Biostrings)
  library(decontam)
  library(vegan)
  library(hillR)
  library(indicspecies)
  library(DESeq2)
  library(S4Vectors)
  
  # modelling
  library(ranger)
  library(randomForest)
  library(caret)
  library(compositions)
  library(igraph)
  
  # plotting
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(ggExtra)
  library(ggplotify)
  library(ggh4x)
  library(patchwork)
  library(cowplot)
  library(scales)
  library(grid)
  library(gridExtra)
  library(eulerr)
  library(ComplexUpset)
  library(ComplexHeatmap)
  library(circlize)
  library(colorspace)
  
  library(devtools)
})

# cowplot masks ggplot2::margin, and a package loaded after it masks the mask,
# so calls to margin() in section 7 are written as ggplot2::margin().


# =============================================================================
#  COLOURS
#  Every palette used below. Same hex values as the original scripts, so no
#  figure changes.
#
#  Two colour vectors are deliberately NOT here, and stay where they were:
#    - mgmt_cols in section 9, because section 5 uses that same name for the
#      multipatt output columns, which are not colours at all
#    - the two vectors in section 14, which carry an extra "NS" entry that the
#      volcano panels need and nothing else should inherit
# =============================================================================
mgmt_order          <- c("Conventional", "No-Till", "Organic")

crop_labels <- c("2018" = "2018\nSoybean",
                 "2019" = "2019\nWheat",
                 "2020" = "2020\nMaize")

ORIGIN_COL    <- "Origin"
origin_levels <- c("Soil","Stem","Leaf","Stem+Leaf","Root")

mgmt_levels  <- make.names(c("Conventional", "No-Till", "Organic"))
mgmt_display <- c("Conventional", "No-Till", "Organic")
name_map     <- setNames(mgmt_display, mgmt_levels)

year_labels <- c(
  "2018" = "2018 — Soybean",
  "2019" = "2019 — Wheat",
  "2020" = "2020 — Maize")

kingdom_labels <- c(
  "bac"      = "Bacterial Microbiome",
  "fun"      = "Fungal Microbiome",
  "combined" = "Combined Microbiome")

compartment_labels <- c(
  "pooled"      = "All Compartments",
  "aboveground" = "Aboveground",
  "belowground" = "Belowground")

years        <- c("2018", "2019", "2020")
kingdoms     <- c("bac", "fun", "combined")
compartments <- c("pooled", "aboveground", "belowground")
comp_suffix  <- c("all", "a", "b")

# Management, growth stage and year. Used by most sections.
management_colors   <- c("Conventional" = "#bd461d",
                         "No-Till"      = "#0d3660",
                         "Organic"      = "#005b40")
growth_stage_colors <- c("Vegetative"    = "#7c2529",
                         "Inflorescence" = "#5d325d",
                         "Reproductive"  = "#ee8868")
year_colors         <- c("2018" = "#6f8740",
                         "2019" = "#bc8400",
                         "2020" = "#ffba2c")

# Section 11 uses the growth stage palette under a second name.
stage_colors <- growth_stage_colors

# Section 2, rarefaction curves, coloured by compartment within each kingdom.
fungi_colors    <- c("aboveground" = "#84BA5B", "belowground" = "#8B4513")
bacteria_colors <- c("aboveground" = "#6A9BD1", "belowground" = "#4A2C7A")

# Section 7, pairwise PERMANOVA heatmap. The same four colours, keyed by
# kingdom and compartment together, so the two figures read as one scheme.
kc_colors <- c("Fungi aboveground"    = "#84BA5B",
               "Fungi belowground"    = "#8B4513",
               "Bacteria aboveground" = "#6A9BD1",
               "Bacteria belowground" = "#4A2C7A")

# Section 9, the aboveground against belowground Euler panels.
col_ab <- "#5b4b8a"
col_bg <- "#b3622d"

# Section 10, sampling niche.
origin_colors <- c("Soil"      = "#6b4226", "Stem" = "#2a9d8f", "Leaf" = "#3f8f3f",
                   "Stem+Leaf" = "#8a6d3b", "Root" = "#c08a3e")

# Section 13, random forest summary figures.
kingdom_colors     <- c("Bacteria" = "#83331D", "Fungi" = "#22569F",
                        "Combined" = "#7A6E3C")
compartment_colors <- c("All Compartments" = "#4A4A4A",
                        "Aboveground"      = "#7B5EA7",
                        "Belowground"      = "#2E7D7B")

# Section 15, the four-way method overlap.
set_cols <- c(ISA = "#6f8740", DA = "#bc8400", RF = "#bd461d", Hub = "#4A2C7A")


kc_order <- names(kc_colors)
# =============================================================================
# =============================================================================
#  SECTION 1
#  Pipelines, subsetting and verification
#
#  Reads the raw bacterial and fungal data, processes each year in its own
#  object, and builds every psr_* and ps_unrare_* object used below.
#  Sets the working directory three times: bacterial inputs, fungal inputs,
#  then the output folder.

# STAGE 1 — BACTERIAL PIPELINE (DADA2 taxonomy)
# PER-YEAR ARCHITECTURE — every step spelled out, once per year.
# ============================================================


# WORKING DIRECTORY

setwd("~/Desktop/Paper/NEW_R_Code/combined_working")


# DADA2 TAXONOMY LOADER  (unchanged)

ranks_7 <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
load_taxonomy <- function(path, ranks = ranks_7) {
  tx <- read.delim(path, header = TRUE, stringsAsFactors = FALSE,
                   na.strings = c("NA", ""))
  stopifnot("ASV_ID" %in% names(tx),
            "Sequence" %in% names(tx),
            all(ranks %in% names(tx)))
  tax <- tx[, ranks, drop = FALSE]
  rownames(tax) <- tx$ASV_ID
  seq <- Biostrings::DNAStringSet(tx$Sequence)
  names(seq) <- tx$ASV_ID
  list(tax = tax, seq = seq)
}

# within-YEAR collision guard (used ONLY for 2019 ag + roots merge)
prefix_asv <- function(ps, tag) {
  taxa_names(ps) <- paste0(taxa_names(ps), "_", tag)   # renames otu + tax + refseq together
  ps
}


# PER-YEAR RAREFACTION DEPTH  (edit any year independently)

rare_depth_bac <- c("2018" = 5000, "2019" = 5000, "2020" = 5000)

# RAREFACTION HELPER  (unchanged math; depth passed per call)

manual_rarefy_bac <- function(physeq, depth) {
  set.seed(123)
  otu_bac           <- as(otu_table(physeq), "matrix")
  taxa_are_rows_val <- taxa_are_rows(otu_table(physeq))
  if (taxa_are_rows_val) otu_bac <- t(otu_bac)
  keep_samples_bac <- rowSums(otu_bac) >= depth
  otu_subset_bac   <- otu_bac[keep_samples_bac, , drop = FALSE]
  rarefied_bac     <- t(vegan::rrarefy(otu_subset_bac, depth))
  otu_rarefied <- if (taxa_are_rows_val) {
    otu_table(rarefied_bac, taxa_are_rows = TRUE)
  } else {
    otu_table(t(rarefied_bac), taxa_are_rows = FALSE)
  }
  ps_out <- phyloseq(
    otu_rarefied,
    tax_table(physeq),
    sample_data(physeq)[rownames(otu_subset_bac), ],
    refseq(physeq)
  )
  prune_taxa(taxa_sums(ps_out) > 0, ps_out)
}


# RAREFACTION DIAGNOSTIC PLOT HELPERS (unchanged)

RareStats_bac <- function(ps) {
  asv_bac  <- as.data.frame(ps@otu_table)
  meta_bac <- as.data.frame(as.matrix(ps@sam_data))
  findoutlier_bac <- function(x) {
    x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x)
  }
  asv_bac %>%
    rownames_to_column("ASV_ID") %>%
    pivot_longer(-ASV_ID, names_to = "Sample_ID", values_to = "Seq_No") %>%
    group_by(Sample_ID) %>%
    summarize(
      Read_No      = sum(Seq_No),
      Singleton_No = sum(Seq_No == 1),
      Goods        = 100 * (1 - Singleton_No / Read_No)
    ) %>%
    mutate(outlier = ifelse(findoutlier_bac(log10(Read_No)), Read_No, NA)) %>%
    ungroup() %>%
    left_join(meta_bac %>% rownames_to_column("Sample_ID"), by = "Sample_ID")
}

PloRareStats_bac <- function(df, year_label = "") {
  require(ggrepel)
  p <- ggarrange(
    df %>% ggplot(aes(x = Read_No)) +
      geom_histogram(binwidth = 5000, fill = "steelblue", color = "steelblue") +
      labs(title = paste(year_label, "Histogram")),
    df %>% ggplot(aes(x = Read_No)) +
      geom_histogram(binwidth = 1000, fill = "steelblue", color = "steelblue") +
      coord_cartesian(xlim = c(0, 10000)) +
      labs(title = "Histogram Zoom"),
    df %>% ggplot(aes(x = Read_No, y = Goods)) +
      geom_point(shape = 1, color = "steelblue") +
      labs(title = "Good's Coverage"),
    df %>% ggplot(aes(x = 1, y = Read_No)) +
      geom_jitter(shape = 1, color = "steelblue") +
      scale_y_log10() + labs(title = "Log10 Jitter"),
    df %>% ggplot(aes(x = 1, y = Read_No)) +
      geom_boxplot(color = "steelblue") + scale_y_log10() +
      geom_text_repel(
        data    = dplyr::filter(df, !is.na(outlier)),
        mapping = aes(x = 1, y = Read_No, label = outlier),
        max.overlaps = 100, size = 3
      ) + labs(title = "Log10 Boxplot"),
    df %>% arrange(Read_No) %>%
      ggplot(aes(x = 1:nrow(.), y = Read_No)) +
      geom_bar(stat = "identity", color = "steelblue") +
      labs(title = "Ranked"),
    ncol = 3, nrow = 2, align = "hv", labels = c("A","B","C","D","E","F")
  )
  print(p)
  p
}

# TAXONOMY LABEL FUNCTION (for downstream plotting; unchanged)

make_tax_labels <- function(ps, prefix) {
  if (is.null(ps)) return(NULL)
  tax     <- as.data.frame(tax_table(ps))
  tax$ASV <- paste0(prefix, "_", rownames(tax))
  kingdom_tag <- ifelse(prefix == "bac", "[B] ", ifelse(prefix == "fun", "[F] ", ""))
  clean <- function(x) {
    x <- sub("_[0-9.]+$", "", x); x <- gsub("_", " ", x)
    x[is.na(x) | x == ""] <- NA; x
  }
  genus <- clean(tax$Genus); fam <- clean(tax$Family); ord <- clean(tax$Order)
  best_name <- ifelse(!is.na(genus), genus, ifelse(!is.na(fam), fam, ifelse(!is.na(ord), ord, NA)))
  tax$label <- ifelse(!is.na(best_name),
                      paste0(kingdom_tag, best_name, " (", rownames(tax), ")"),
                      paste0(kingdom_tag, rownames(tax)))
  tax[, c("ASV", "label")]
}


####  YEAR 2018 — SOYBEAN — full pipeline in isolation  ###

# ----- STEP 1 (2018): LOAD RAW DATA -----
cat("\n=== [2018] Loading Bacterial Data ===\n")
ASV_bac_18  <- read.delim("otutable_18_UNOISE_180bp.txt", row.names = 1)
Meta_bac_18 <- read.delim("Bac_18_Mapping.txt")
tax18       <- load_taxonomy("asv18_180bp_taxonomy.tsv")
Tax_bac_18  <- tax18$tax
fasta_bac_18 <- tax18$seq

common_samples_bac_18 <- intersect(colnames(ASV_bac_18), Meta_bac_18$SampleID)
Meta_bac_filtered_18  <- Meta_bac_18[Meta_bac_18$SampleID %in% common_samples_bac_18, ]
rownames(Meta_bac_filtered_18) <- Meta_bac_filtered_18$SampleID
ASV_bac_filtered_18   <- ASV_bac_18[, common_samples_bac_18]
cat("[2018] tax coverage:",
    length(intersect(rownames(ASV_bac_filtered_18), rownames(Tax_bac_18))),
    "/", nrow(ASV_bac_filtered_18), "ASVs\n")

ps_18 <- phyloseq(
  otu_table(ASV_bac_filtered_18, taxa_are_rows = TRUE),
  sample_data(Meta_bac_filtered_18),
  tax_table(as.matrix(Tax_bac_18)),
  fasta_bac_18
)
ps_18 <- subset_samples(ps_18, Experiment %in% c("obj_1"))
cat("[2018] loaded:", nsamples(ps_18), "samples,", ntaxa(ps_18), "ASVs\n")

# ----- STEP 4 (2018): FIX METADATA -----
cat("\n=== [2018] Fixing Metadata ===\n")
sd <- as.data.frame(sample_data(ps_18))
# remove control growth stages (none expected in 2018, but mirror the rule)
remove_ctrl <- sd$Growth_stage %in% c("C1", "CO", "NE", "PO")
cat("[2018] removing control-stage samples:", sum(remove_ctrl, na.rm = TRUE), "\n")
ps_18 <- prune_samples(!remove_ctrl, ps_18)
ps_18 <- prune_taxa(taxa_sums(ps_18) > 0, ps_18)

sd <- as.data.frame(sample_data(ps_18))
# spelling
sd$Growth_Stage_Description <- gsub("Inflorescense", "Inflorescence",
                                    sd$Growth_Stage_Description, ignore.case = TRUE)
# 2018 crop-code -> common stage
is_2018 <- !is.na(sd$Year) & sd$Year == "2018"
if (sum(is_2018) > 0) {
  sd$Growth_Stage_Description[is_2018] <- dplyr::recode(
    sd$Growth_stage[is_2018],
    "V2" = "Vegetative", "R2" = "Inflorescence", "R6" = "Reproductive",
    .default = sd$Growth_Stage_Description[is_2018]
  )
  cat("[2018] growth stages rebuilt (V2/R2/R6)\n")
}
sample_data(ps_18) <- sd

# drop off-design management + control stages
sd <- as.data.frame(sample_data(ps_18))
keep_samples <- !(sd$Management %in% c("Control", "Low-Input", "Other")) &
  !(sd$Growth_Stage_Description %in% c("CONTROL", "Control_or_Check", "N/A"))
ps_18 <- prune_samples(keep_samples, ps_18)
ps_18 <- prune_taxa(taxa_sums(ps_18) > 0, ps_18)

# standardize compartment + factor levels
sd <- as.data.frame(sample_data(ps_18))
sd$Compartment <- gsub("-", "", tolower(as.character(sd$Compartment)))
sd$Management  <- factor(sd$Management, levels = mgmt_order)
sd$Year        <- factor(sd$Year, levels = c("2018", "2019", "2020"))
sample_data(ps_18) <- sd
cat("[2018] after metadata cleanup:", nsamples(ps_18), "samples\n")
print(table(Compartment = sd$Compartment, Stage = sd$Growth_Stage_Description))

# ----- STEP 5 (2018): REMOVE NON-BACTERIAL TAXA -----
cat("\n=== [2018] Removing Non-Bacterial Taxa ===\n")
cat("[2018] ASVs before:", ntaxa(ps_18), "\n")
ps_18 <- subset_taxa(
  ps_18,
  !(Order  %in% c("Chloroplast", "Mitochondria")) &
    !(Family %in% c("Chloroplast", "Mitochondria")) &
    !(Phylum %in% c("Chloroplast", "Mitochondria")) &
    !(Kingdom %in% c("Anthophyta", "Alveolata", "Ichthyosporia",
                     "Protista", "Metazoa", "Rhizaria", "Viridiplantae")) &
    Kingdom %in% c("Bacteria", "Archaea")
)
cat("[2018] ASVs after:", ntaxa(ps_18), "\n")

# ----- STEP 6 (2018): CONTAMINANT REMOVAL & ABUNDANCE FILTER -----
cat("\n=== [2018] Contaminant Removal ===\n")
sample_data(ps_18)$is.neg <- sample_data(ps_18)$Sample_or_Control == "Control_Sample"
if (any(sample_data(ps_18)$is.neg, na.rm = TRUE)) {
  contamdf_18 <- isContaminant(ps_18, method = "prevalence", neg = "is.neg")
  cat("[2018] contaminants identified:", sum(contamdf_18$contaminant), "\n")
  ps_18 <- prune_taxa(!contamdf_18$contaminant, ps_18)
} else {
  cat("[2018] no negative controls present — skipping decontam\n")
}
ps_18 <- subset_samples(ps_18, Sample_or_Control == "True_Sample")
cat("[2018] samples after removing controls:", nsamples(ps_18), "\n")

# node removal
sd_check <- as.data.frame(sample_data(ps_18))
if ("node" %in% sd_check$Origin) {
  cat("[2018] removing", sum(sd_check$Origin == "node", na.rm = TRUE), "node samples\n")
  ps_18 <- subset_samples(ps_18, Origin != "node")
  ps_18 <- prune_taxa(taxa_sums(ps_18) > 0, ps_18)
} else cat("[2018] no node samples\n")

# abundance filtering (WITHIN 2018)
cat("\n=== [2018] Abundance Filtering (within-year) ===\n")
ps_18 <- prune_samples(sample_sums(ps_18) >= 1000, ps_18)
cat("[2018] samples after <1000-read drop:", nsamples(ps_18), "\n")
m18 <- as(otu_table(ps_18), "matrix"); if (!taxa_are_rows(ps_18)) m18 <- t(m18)
keep_asv_18 <- (rowSums(m18 >= 0) >= 3) & (rowSums(m18) >= 10)
cat("[2018] rule prev>=3 & reads>=10 | kept:", sum(keep_asv_18),
    "removed:", sum(!keep_asv_18),
    sprintf("(%.1f%% removed)\n", 100 * mean(!keep_asv_18)))
ps_18 <- prune_taxa(keep_asv_18, ps_18)
ps_18 <- subset_samples(ps_18,
                        !Management %in% c("Control", "Low-Input") &
                          !Growth_Stage_Description %in% c("CONTROL"))
cat("[2018] after final cleanup:", nsamples(ps_18), "samples,", ntaxa(ps_18), "ASVs\n")

# ----- STEP 7 (2018): KNOWN OUTLIER SAMPLES -----
cat("\n=== [2018] Removing Known Outlier Samples ===\n")
outlier_samples <- c("Sample1435", "Sample1534")
pres_18 <- outlier_samples[outlier_samples %in% sample_names(ps_18)]
if (length(pres_18) > 0) {
  cat("[2018] removing:", paste(pres_18, collapse = ", "), "\n")
  ps_18 <- prune_samples(!sample_names(ps_18) %in% pres_18, ps_18)
  ps_18 <- prune_taxa(taxa_sums(ps_18) > 0, ps_18)
} else cat("[2018] no known outliers present\n")

# ----- STEP 8 (2018): COMPOSITIONAL OUTLIERS (SF > 10) -----
cat("\n=== [2018] Removing Compositional Outliers ===\n")
sample_data(ps_18)$Management <- factor(sample_data(ps_18)$Management, levels = mgmt_order)
dds_18 <- phyloseq_to_deseq2(ps_18, ~ Management)
dds_18 <- estimateSizeFactors(dds_18, type = "poscounts")
sf_18  <- sizeFactors(dds_18)
out_18 <- names(sf_18)[sf_18 > 10]
cat("[2018] total samples:", nsamples(ps_18), "| SF>10:", length(out_18), "\n")
if (length(out_18) > 0) {
  md <- data.frame(sample_data(ps_18))
  info <- md[out_18, c("Year","Management","Compartment","Growth_Stage_Description")]
  info$SizeFactor <- round(sf_18[out_18], 2); print(info)
  ps_18 <- prune_samples(!sample_names(ps_18) %in% out_18, ps_18)
  ps_18 <- prune_taxa(taxa_sums(ps_18) > 0, ps_18)
  cat("[2018] after removal:", nsamples(ps_18), "samples\n")
} else cat("[2018] no compositional outliers\n")

# ----- STEP 9 (2018): RAREFACTION -----
cat("\n=== [2018] Rarefaction Statistics ===\n")
rare_stats_18 <- RareStats_bac(ps_18)
PloRareStats_bac(rare_stats_18, "2018")
d18 <- rare_depth_bac[["2018"]]
cat("\n=== [2018] Rarefying to", d18, "reads ===\n")
ps_18 <- prune_samples(sample_sums(ps_18) >= d18, ps_18)
# NOTE: ps_18 is now the cleaned UNRAREFIED, depth-matched object. It
# persists in the global env (manual_rarefy_bac returns a NEW object and
# never overwrites ps_18), so it is the bacterial DESeq2 DA input.
physeq_rare_bac_2018 <- manual_rarefy_bac(ps_18, depth = d18)
sample_data(physeq_rare_bac_2018)$Management <- factor(sample_data(physeq_rare_bac_2018)$Management, levels = mgmt_order)
sample_data(physeq_rare_bac_2018)$Year       <- factor(sample_data(physeq_rare_bac_2018)$Year, levels = c("2018","2019","2020"))
dep18 <- sample_sums(physeq_rare_bac_2018)
cat("[2018] Min:", min(dep18), "Max:", max(dep18), "all==", d18, ":", all(dep18 == d18),
    "|", nsamples(physeq_rare_bac_2018), "samples,", ntaxa(physeq_rare_bac_2018), "ASVs\n")


# ###  YEAR 2019 — WHEAT — full pipeline in isolation    ###
# ###  (aboveground + roots = two runs of THIS year)     ###

# ----- STEP 1 (2019): LOAD ABOVEGROUND RUN -----
cat("\n=== [2019] Loading Aboveground Data ===\n")
ASV_bac_19  <- read.delim("otutable_19_UNOISE_180bp.txt", row.names = 1)
Meta_bac_19 <- read.delim("bac_19_mapping_2.txt")
tax19       <- load_taxonomy("asv19_180bp_taxonomy.tsv")
Tax_bac_19  <- tax19$tax
fasta_bac_19 <- tax19$seq

common_samples_bac_19 <- intersect(colnames(ASV_bac_19), Meta_bac_19$SampleID)
Meta_bac_filtered_19  <- Meta_bac_19[Meta_bac_19$SampleID %in% common_samples_bac_19, ]
rownames(Meta_bac_filtered_19) <- Meta_bac_filtered_19$SampleID
ASV_bac_filtered_19   <- ASV_bac_19[, common_samples_bac_19]
cat("[2019 ag] tax coverage:",
    length(intersect(rownames(ASV_bac_filtered_19), rownames(Tax_bac_19))),
    "/", nrow(ASV_bac_filtered_19), "ASVs\n")

ps_19 <- phyloseq(
  otu_table(ASV_bac_filtered_19, taxa_are_rows = TRUE),
  sample_data(Meta_bac_filtered_19),
  tax_table(as.matrix(Tax_bac_19)),
  fasta_bac_19
)
ps_19 <- subset_samples(ps_19, Experiment %in% c("obj_1"))
cat("[2019 ag] loaded:", nsamples(ps_19), "samples,", ntaxa(ps_19), "ASVs\n")

# ----- STEP 1 (2019): LOAD ROOTS RUN -----
cat("\n=== [2019] Loading Root Data ===\n")
ASV_19_roots_bac  <- read.delim("root_19_otutable_UNOISE_180bp.txt", row.names = 1)
Meta_19_roots_bac <- read.delim("root_bac_19_mapping.txt")
tax19r            <- load_taxonomy("asv_180bp_19root_taxonomy.tsv")
Tax_19_roots_bac  <- tax19r$tax
fasta_19_roots_bac <- tax19r$seq

common_samples_roots   <- intersect(colnames(ASV_19_roots_bac), Meta_19_roots_bac$SampleID)
Meta_19_roots_filtered <- Meta_19_roots_bac[Meta_19_roots_bac$SampleID %in% common_samples_roots, ]
rownames(Meta_19_roots_filtered) <- Meta_19_roots_filtered$SampleID
ASV_19_roots_filtered  <- ASV_19_roots_bac[, common_samples_roots]
cat("[2019 root] tax coverage:",
    length(intersect(rownames(ASV_19_roots_filtered), rownames(Tax_19_roots_bac))),
    "/", nrow(ASV_19_roots_filtered), "ASVs\n")

ps_19_roots <- phyloseq(
  otu_table(ASV_19_roots_filtered, taxa_are_rows = TRUE),
  sample_data(Meta_19_roots_filtered),
  tax_table(as.matrix(Tax_19_roots_bac)),
  fasta_19_roots_bac
)
cat("[2019 root] loaded:", nsamples(ps_19_roots), "samples,", ntaxa(ps_19_roots), "ASVs\n")

# ----- STEP 2 (2019 ONLY): MERGE THE TWO 2019 RUNS, COLLISION-GUARDED -----
# These are two UNOISE runs of the SAME year -> ids both restart at Zotu1,
# so tag _ag / _root before merging so the two runs cannot fuse.
cat("\n=== [2019] Merging aboveground + roots (within-year, prefixed) ===\n")
ps_19       <- prefix_asv(ps_19,       "ag")
ps_19_roots <- prefix_asv(ps_19_roots, "root")

main_cols <- colnames(sample_data(ps_19))
root_cols <- colnames(sample_data(ps_19_roots))
miss_root <- setdiff(main_cols, root_cols)
if (length(miss_root) > 0) {
  x <- as.data.frame(sample_data(ps_19_roots)); for (c in miss_root) x[[c]] <- NA
  sample_data(ps_19_roots) <- x
}
miss_main <- setdiff(root_cols, main_cols)
if (length(miss_main) > 0) {
  x <- as.data.frame(sample_data(ps_19)); for (c in miss_main) x[[c]] <- NA
  sample_data(ps_19) <- x
}
ps_19 <- merge_phyloseq(ps_19, ps_19_roots)
cat("[2019] merged (ag + roots):", nsamples(ps_19), "samples,", ntaxa(ps_19), "ASVs\n")

## within-year collision verification: ids must be ag- or root-unique,
## never appearing in BOTH ag and root sample sets at once.
.m <- as(otu_table(ps_19), "matrix"); if (!taxa_are_rows(ps_19)) .m <- t(.m)
.src <- ifelse(grepl("_ag$", rownames(.m)), "ag",
               ifelse(grepl("_root$", rownames(.m)), "root", "??"))
cat("[2019] id source tags (want only ag/root, no ??):\n"); print(table(.src)); rm(.m, .src)

# ----- STEP 4 (2019): FIX METADATA -----
cat("\n=== [2019] Fixing Metadata ===\n")
sd <- as.data.frame(sample_data(ps_19))
remove_ctrl <- sd$Growth_stage %in% c("C1", "CO", "NE", "PO")
cat("[2019] removing control-stage samples (C1/CO/NE/PO):", sum(remove_ctrl, na.rm = TRUE), "\n")
ps_19 <- prune_samples(!remove_ctrl, ps_19)
ps_19 <- prune_taxa(taxa_sums(ps_19) > 0, ps_19)

sd <- as.data.frame(sample_data(ps_19))
sd$Growth_Stage_Description <- gsub("Inflorescense", "Inflorescence",
                                    sd$Growth_Stage_Description, ignore.case = TRUE)
is_2019 <- !is.na(sd$Year) & sd$Year == "2019"
if (sum(is_2019) > 0) {
  sd$Growth_Stage_Description[is_2019] <- dplyr::recode(
    sd$Growth_stage[is_2019],
    "C2" = "Vegetative", "C3" = "Inflorescence", "C4" = "Reproductive",
    .default = sd$Growth_Stage_Description[is_2019]
  )
  cat("[2019] growth stages rebuilt (C2/C3/C4)\n")
}
sample_data(ps_19) <- sd

sd <- as.data.frame(sample_data(ps_19))
keep_samples <- !(sd$Management %in% c("Control", "Low-Input", "Other")) &
  !(sd$Growth_Stage_Description %in% c("CONTROL", "Control_or_Check", "N/A"))
ps_19 <- prune_samples(keep_samples, ps_19)
ps_19 <- prune_taxa(taxa_sums(ps_19) > 0, ps_19)

sd <- as.data.frame(sample_data(ps_19))
sd$Compartment <- gsub("-", "", tolower(as.character(sd$Compartment)))
sd$Management  <- factor(sd$Management, levels = mgmt_order)
sd$Year        <- factor(sd$Year, levels = c("2018", "2019", "2020"))
sample_data(ps_19) <- sd
cat("[2019] after metadata cleanup:", nsamples(ps_19), "samples\n")
print(table(Compartment = sd$Compartment, Stage = sd$Growth_Stage_Description))

# ----- STEP 5 (2019): REMOVE NON-BACTERIAL TAXA -----
cat("\n=== [2019] Removing Non-Bacterial Taxa ===\n")
cat("[2019] ASVs before:", ntaxa(ps_19), "\n")
ps_19 <- subset_taxa(
  ps_19,
  !(Order  %in% c("Chloroplast", "Mitochondria")) &
    !(Family %in% c("Chloroplast", "Mitochondria")) &
    !(Phylum %in% c("Chloroplast", "Mitochondria")) &
    !(Kingdom %in% c("Anthophyta", "Alveolata", "Ichthyosporia",
                     "Protista", "Metazoa", "Rhizaria", "Viridiplantae")) &
    Kingdom %in% c("Bacteria", "Archaea")
)
cat("[2019] ASVs after:", ntaxa(ps_19), "\n")

# ----- STEP 6 (2019): CONTAMINANT REMOVAL & ABUNDANCE FILTER -----
cat("\n=== [2019] Contaminant Removal ===\n")
sample_data(ps_19)$is.neg <- sample_data(ps_19)$Sample_or_Control == "Control_Sample"
if (any(sample_data(ps_19)$is.neg, na.rm = TRUE)) {
  contamdf_19 <- isContaminant(ps_19, method = "prevalence", neg = "is.neg")
  cat("[2019] contaminants identified:", sum(contamdf_19$contaminant), "\n")
  ps_19 <- prune_taxa(!contamdf_19$contaminant, ps_19)
} else {
  cat("[2019] no negative controls present — skipping decontam\n")
}
ps_19 <- subset_samples(ps_19, Sample_or_Control == "True_Sample")
cat("[2019] samples after removing controls:", nsamples(ps_19), "\n")

sd_check <- as.data.frame(sample_data(ps_19))
if ("node" %in% sd_check$Origin) {
  cat("[2019] removing", sum(sd_check$Origin == "node", na.rm = TRUE), "node samples\n")
  ps_19 <- subset_samples(ps_19, Origin != "node")
  ps_19 <- prune_taxa(taxa_sums(ps_19) > 0, ps_19)
} else cat("[2019] no node samples\n")

cat("\n=== [2019] Abundance Filtering (within-year) ===\n")
ps_19 <- prune_samples(sample_sums(ps_19) >= 1000, ps_19)
cat("[2019] samples after <1000-read drop:", nsamples(ps_19), "\n")
m19 <- as(otu_table(ps_19), "matrix"); if (!taxa_are_rows(ps_19)) m19 <- t(m19)
keep_asv_19 <- (rowSums(m19 >= 0) >= 3) & (rowSums(m19) >= 10)
cat("[2019] rule prev>=3 & reads>=10 | kept:", sum(keep_asv_19),
    "removed:", sum(!keep_asv_19),
    sprintf("(%.1f%% removed)\n", 100 * mean(!keep_asv_19)))
ps_19 <- prune_taxa(keep_asv_19, ps_19)
ps_19 <- subset_samples(ps_19,
                        !Management %in% c("Control", "Low-Input") &
                          !Growth_Stage_Description %in% c("CONTROL"))
cat("[2019] after final cleanup:", nsamples(ps_19), "samples,", ntaxa(ps_19), "ASVs\n")

# ----- STEP 7 (2019): KNOWN OUTLIER SAMPLES -----
cat("\n=== [2019] Removing Known Outlier Samples ===\n")
pres_19 <- outlier_samples[outlier_samples %in% sample_names(ps_19)]
if (length(pres_19) > 0) {
  cat("[2019] removing:", paste(pres_19, collapse = ", "), "\n")
  ps_19 <- prune_samples(!sample_names(ps_19) %in% pres_19, ps_19)
  ps_19 <- prune_taxa(taxa_sums(ps_19) > 0, ps_19)
} else cat("[2019] no known outliers present\n")

# ----- STEP 8 (2019): COMPOSITIONAL OUTLIERS (SF > 10) -----
cat("\n=== [2019] Removing Compositional Outliers ===\n")
sample_data(ps_19)$Management <- factor(sample_data(ps_19)$Management, levels = mgmt_order)
dds_19 <- phyloseq_to_deseq2(ps_19, ~ Management)
dds_19 <- estimateSizeFactors(dds_19, type = "poscounts")
sf_19  <- sizeFactors(dds_19)
out_19 <- names(sf_19)[sf_19 > 10]
cat("[2019] total samples:", nsamples(ps_19), "| SF>10:", length(out_19), "\n")
if (length(out_19) > 0) {
  md <- data.frame(sample_data(ps_19))
  info <- md[out_19, c("Year","Management","Compartment","Growth_Stage_Description")]
  info$SizeFactor <- round(sf_19[out_19], 2); print(info)
  ps_19 <- prune_samples(!sample_names(ps_19) %in% out_19, ps_19)
  ps_19 <- prune_taxa(taxa_sums(ps_19) > 0, ps_19)
  cat("[2019] after removal:", nsamples(ps_19), "samples\n")
} else cat("[2019] no compositional outliers\n")

# ----- STEP 9 (2019): RAREFACTION -----
cat("\n=== [2019] Rarefaction Statistics ===\n")
rare_stats_19 <- RareStats_bac(ps_19)
PloRareStats_bac(rare_stats_19, "2019")
d19 <- rare_depth_bac[["2019"]]
cat("\n=== [2019] Rarefying to", d19, "reads ===\n")
ps_19 <- prune_samples(sample_sums(ps_19) >= d19, ps_19)
# NOTE: ps_19 is now the cleaned UNRAREFIED, depth-matched object (DA input).
physeq_rare_bac_2019 <- manual_rarefy_bac(ps_19, depth = d19)
sample_data(physeq_rare_bac_2019)$Management <- factor(sample_data(physeq_rare_bac_2019)$Management, levels = mgmt_order)
sample_data(physeq_rare_bac_2019)$Year       <- factor(sample_data(physeq_rare_bac_2019)$Year, levels = c("2018","2019","2020"))
dep19 <- sample_sums(physeq_rare_bac_2019)
cat("[2019] Min:", min(dep19), "Max:", max(dep19), "all==", d19, ":", all(dep19 == d19),
    "|", nsamples(physeq_rare_bac_2019), "samples,", ntaxa(physeq_rare_bac_2019), "ASVs\n")


# ###  YEAR 2020 — MAIZE — full pipeline in isolation    ###

# ----- STEP 1 (2020): LOAD RAW DATA -----
cat("\n=== [2020] Loading Bacterial Data ===\n")
ASV_bac_20  <- read.delim("otutable_20_UNOISE_180bp.txt", row.names = 1)
Meta_bac_20 <- read.delim("bac_20_mapping_2.txt")
tax20       <- load_taxonomy("asv20_180bp_taxonomy.tsv")
Tax_bac_20  <- tax20$tax
fasta_bac_20 <- tax20$seq

common_samples_bac_20 <- intersect(colnames(ASV_bac_20), Meta_bac_20$SampleID)
Meta_bac_filtered_20  <- Meta_bac_20[Meta_bac_20$SampleID %in% common_samples_bac_20, ]
rownames(Meta_bac_filtered_20) <- Meta_bac_filtered_20$SampleID
ASV_bac_filtered_20   <- ASV_bac_20[, common_samples_bac_20]
cat("[2020] tax coverage:",
    length(intersect(rownames(ASV_bac_filtered_20), rownames(Tax_bac_20))),
    "/", nrow(ASV_bac_filtered_20), "ASVs\n")

ps_20 <- phyloseq(
  otu_table(ASV_bac_filtered_20, taxa_are_rows = TRUE),
  sample_data(Meta_bac_filtered_20),
  tax_table(as.matrix(Tax_bac_20)),
  fasta_bac_20
)
ps_20 <- subset_samples(ps_20, Experiment %in% c("obj_1"))
cat("[2020] loaded:", nsamples(ps_20), "samples,", ntaxa(ps_20), "ASVs\n")

# ----- STEP 4 (2020): FIX METADATA -----
cat("\n=== [2020] Fixing Metadata ===\n")
sd <- as.data.frame(sample_data(ps_20))
remove_ctrl <- sd$Growth_stage %in% c("C1", "CO", "NE", "PO")
cat("[2020] removing control-stage samples:", sum(remove_ctrl, na.rm = TRUE), "\n")
ps_20 <- prune_samples(!remove_ctrl, ps_20)
ps_20 <- prune_taxa(taxa_sums(ps_20) > 0, ps_20)

sd <- as.data.frame(sample_data(ps_20))
sd$Growth_Stage_Description <- gsub("Inflorescense", "Inflorescence",
                                    sd$Growth_Stage_Description, ignore.case = TRUE)
is_2020 <- !is.na(sd$Year) & sd$Year == "2020"
if (sum(is_2020) > 0) {
  sd$Growth_Stage_Description[is_2020] <- dplyr::recode(
    sd$Growth_stage[is_2020],
    "V4" = "Vegetative", "V5" = "Vegetative",
    "VT" = "Inflorescence", "R4" = "Reproductive",
    .default = sd$Growth_Stage_Description[is_2020]
  )
  cat("[2020] growth stages rebuilt (V4+V5=Vegetative, VT, R4)\n")
}
sample_data(ps_20) <- sd

sd <- as.data.frame(sample_data(ps_20))
keep_samples <- !(sd$Management %in% c("Control", "Low-Input", "Other")) &
  !(sd$Growth_Stage_Description %in% c("CONTROL", "Control_or_Check", "N/A"))
ps_20 <- prune_samples(keep_samples, ps_20)
ps_20 <- prune_taxa(taxa_sums(ps_20) > 0, ps_20)

sd <- as.data.frame(sample_data(ps_20))
sd$Compartment <- gsub("-", "", tolower(as.character(sd$Compartment)))
sd$Management  <- factor(sd$Management, levels = mgmt_order)
sd$Year        <- factor(sd$Year, levels = c("2018", "2019", "2020"))
sample_data(ps_20) <- sd
cat("[2020] after metadata cleanup:", nsamples(ps_20), "samples\n")
print(table(Compartment = sd$Compartment, Stage = sd$Growth_Stage_Description))

# ----- STEP 5 (2020): REMOVE NON-BACTERIAL TAXA -----
cat("\n=== [2020] Removing Non-Bacterial Taxa ===\n")
cat("[2020] ASVs before:", ntaxa(ps_20), "\n")
ps_20 <- subset_taxa(
  ps_20,
  !(Order  %in% c("Chloroplast", "Mitochondria")) &
    !(Family %in% c("Chloroplast", "Mitochondria")) &
    !(Phylum %in% c("Chloroplast", "Mitochondria")) &
    !(Kingdom %in% c("Anthophyta", "Alveolata", "Ichthyosporia",
                     "Protista", "Metazoa", "Rhizaria", "Viridiplantae")) &
    Kingdom %in% c("Bacteria", "Archaea")
)
cat("[2020] ASVs after:", ntaxa(ps_20), "\n")

# ----- STEP 6 (2020): CONTAMINANT REMOVAL & ABUNDANCE FILTER -----
cat("\n=== [2020] Contaminant Removal ===\n")
sample_data(ps_20)$is.neg <- sample_data(ps_20)$Sample_or_Control == "Control_Sample"
if (any(sample_data(ps_20)$is.neg, na.rm = TRUE)) {
  contamdf_20 <- isContaminant(ps_20, method = "prevalence", neg = "is.neg")
  cat("[2020] contaminants identified:", sum(contamdf_20$contaminant), "\n")
  ps_20 <- prune_taxa(!contamdf_20$contaminant, ps_20)
} else {
  cat("[2020] no negative controls present — skipping decontam\n")
}
ps_20 <- subset_samples(ps_20, Sample_or_Control == "True_Sample")
cat("[2020] samples after removing controls:", nsamples(ps_20), "\n")

sd_check <- as.data.frame(sample_data(ps_20))
if ("node" %in% sd_check$Origin) {
  cat("[2020] removing", sum(sd_check$Origin == "node", na.rm = TRUE), "node samples\n")
  ps_20 <- subset_samples(ps_20, Origin != "node")
  ps_20 <- prune_taxa(taxa_sums(ps_20) > 0, ps_20)
} else cat("[2020] no node samples\n")

cat("\n=== [2020] Abundance Filtering (within-year) ===\n")
ps_20 <- prune_samples(sample_sums(ps_20) >= 1000, ps_20)
cat("[2020] samples after <1000-read drop:", nsamples(ps_20), "\n")
m20 <- as(otu_table(ps_20), "matrix"); if (!taxa_are_rows(ps_20)) m20 <- t(m20)
keep_asv_20 <- (rowSums(m20 >= 0) >= 3) & (rowSums(m20) >= 10)
cat("[2020] rule prev>=3 & reads>=10 | kept:", sum(keep_asv_20),
    "removed:", sum(!keep_asv_20),
    sprintf("(%.1f%% removed)\n", 100 * mean(!keep_asv_20)))
ps_20 <- prune_taxa(keep_asv_20, ps_20)
ps_20 <- subset_samples(ps_20,
                        !Management %in% c("Control", "Low-Input") &
                          !Growth_Stage_Description %in% c("CONTROL"))
cat("[2020] after final cleanup:", nsamples(ps_20), "samples,", ntaxa(ps_20), "ASVs\n")

# ----- STEP 7 (2020): KNOWN OUTLIER SAMPLES -----
cat("\n=== [2020] Removing Known Outlier Samples ===\n")
pres_20 <- outlier_samples[outlier_samples %in% sample_names(ps_20)]
if (length(pres_20) > 0) {
  cat("[2020] removing:", paste(pres_20, collapse = ", "), "\n")
  ps_20 <- prune_samples(!sample_names(ps_20) %in% pres_20, ps_20)
  ps_20 <- prune_taxa(taxa_sums(ps_20) > 0, ps_20)
} else cat("[2020] no known outliers present\n")

# ----- STEP 8 (2020): COMPOSITIONAL OUTLIERS (SF > 10) -----
cat("\n=== [2020] Removing Compositional Outliers ===\n")
sample_data(ps_20)$Management <- factor(sample_data(ps_20)$Management, levels = mgmt_order)
dds_20 <- phyloseq_to_deseq2(ps_20, ~ Management)
dds_20 <- estimateSizeFactors(dds_20, type = "poscounts")
sf_20  <- sizeFactors(dds_20)
out_20 <- names(sf_20)[sf_20 > 10]
cat("[2020] total samples:", nsamples(ps_20), "| SF>10:", length(out_20), "\n")
if (length(out_20) > 0) {
  md <- data.frame(sample_data(ps_20))
  info <- md[out_20, c("Year","Management","Compartment","Growth_Stage_Description")]
  info$SizeFactor <- round(sf_20[out_20], 2); print(info)
  ps_20 <- prune_samples(!sample_names(ps_20) %in% out_20, ps_20)
  ps_20 <- prune_taxa(taxa_sums(ps_20) > 0, ps_20)
  cat("[2020] after removal:", nsamples(ps_20), "samples\n")
} else cat("[2020] no compositional outliers\n")

# ----- STEP 9 (2020): RAREFACTION -----
cat("\n=== [2020] Rarefaction Statistics ===\n")
rare_stats_20 <- RareStats_bac(ps_20)
PloRareStats_bac(rare_stats_20, "2020")
d20 <- rare_depth_bac[["2020"]]
cat("\n=== [2020] Rarefying to", d20, "reads ===\n")
ps_20 <- prune_samples(sample_sums(ps_20) >= d20, ps_20)
# NOTE: ps_20 is now the cleaned UNRAREFIED, depth-matched object (DA input).
physeq_rare_bac_2020 <- manual_rarefy_bac(ps_20, depth = d20)
sample_data(physeq_rare_bac_2020)$Management <- factor(sample_data(physeq_rare_bac_2020)$Management, levels = mgmt_order)
sample_data(physeq_rare_bac_2020)$Year       <- factor(sample_data(physeq_rare_bac_2020)$Year, levels = c("2018","2019","2020"))
dep20 <- sample_sums(physeq_rare_bac_2020)
cat("[2020] Min:", min(dep20), "Max:", max(dep20), "all==", d20, ":", all(dep20 == d20),
    "|", nsamples(physeq_rare_bac_2020), "samples,", ntaxa(physeq_rare_bac_2020), "ASVs\n")


# ###  STEP 10 — SUBSETTING ###

cat("\n=== Creating Subsets (per year) ===\n")
year_objs <- list("2018" = physeq_rare_bac_2018,
                  "2019" = physeq_rare_bac_2019,
                  "2020" = physeq_rare_bac_2020)

compartments_bac <- c("aboveground", "belowground"); comp_suffix_bac <- c("a", "b")
managements_bac  <- c("Conventional", "No-Till", "Organic"); m_suffix_bac <- c("C", "N", "O")

for (yr in names(year_objs)) {
  psr_year_bac <- year_objs[[yr]]
  assign(paste0("psr_", yr, "_bac"), psr_year_bac)
  
  for (i in seq_along(compartments_bac)) {
    psr_comp_bac <- subset_samples(psr_year_bac, Compartment == compartments_bac[i])
    psr_comp_bac <- prune_taxa(taxa_sums(psr_comp_bac) > 0, psr_comp_bac)
    assign(paste0("psr_", yr, "_", comp_suffix_bac[i], "_bac"), psr_comp_bac)
    
    for (j in seq_along(managements_bac)) {
      psr_mgmt_bac <- subset_samples(psr_comp_bac, Management == managements_bac[j])
      psr_mgmt_bac <- prune_taxa(taxa_sums(psr_mgmt_bac) > 0, psr_mgmt_bac)
      assign(paste0("psr_", yr, "_", m_suffix_bac[j], "_", comp_suffix_bac[i], "_bac"),
             psr_mgmt_bac)
    }
  }
  gc()
}

cat("\n===========================================================\n")
cat("  BACTERIAL PIPELINE COMPLETE (PER-YEAR)\n")
cat("===========================================================\n")
for (yr in names(year_objs)) {
  ps <- year_objs[[yr]]
  cat(sprintf("\n[%s] %d samples, %d ASVs\n", yr, nsamples(ps), ntaxa(ps)))
  print(table(Year = sample_data(ps)$Year, Compartment = sample_data(ps)$Compartment))
  cat("  Growth stage:\n")
  print(table(sample_data(ps)$Growth_Stage_Description))
}

## PROOF 1 — no cross-year ASV id collisions are even POSSIBLE.
cat("\n--- PROOF: cross-year id overlap is irrelevant (years never merged) ---\n")
cat("2018 vs 2019 shared ids:", length(intersect(taxa_names(physeq_rare_bac_2018), taxa_names(physeq_rare_bac_2019))), "\n")
cat("2018 vs 2020 shared ids:", length(intersect(taxa_names(physeq_rare_bac_2018), taxa_names(physeq_rare_bac_2020))), "\n")
cat("2019 vs 2020 shared ids:", length(intersect(taxa_names(physeq_rare_bac_2019), taxa_names(physeq_rare_bac_2020))), "\n")
cat("  (these can be >0 and it does NOT matter — no object contains two years)\n")

## PROOF 2 — refseq survived into the slices (sequence-based Euler needs it).
cat("\n--- PROOF: refseq present in slices ---\n")
for (nm in c("psr_2018_a_bac","psr_2019_a_bac","psr_2020_a_bac",
             "psr_2018_b_bac","psr_2019_b_bac","psr_2020_b_bac")) {
  ps <- get(nm)
  rs <- tryCatch(length(refseq(ps)), error = function(e) NA)
  cat(sprintf("  %-18s seqs: %s\n", nm, rs))
}

## PROOF 3 — every slice is at its year's rarefaction depth (flat).
cat("\n--- PROOF: per-slice depth is flat at the year's target ---\n")
for (yr in names(year_objs)) {
  d <- rare_depth_bac[[yr]]
  for (suf in c("a","b")) {
    ps <- get(paste0("psr_", yr, "_", suf, "_bac"))
    if (nsamples(ps) > 0) {
      dep <- sample_sums(ps)
      cat(sprintf("  psr_%s_%s_bac: n=%d depth all==%d: %s\n",
                  yr, suf, nsamples(ps), d, all(dep == d)))
    } else {
      cat(sprintf("  psr_%s_%s_bac: EMPTY\n", yr, suf))
    }
  }
}

cat("\nObjects available:\n")
cat("  physeq_rare_bac_2018 / _2019 / _2020   (per-year rarefied finals)\n")
cat("  psr_<yr>_bac, psr_<yr>_<a|b>_bac, psr_<yr>_<C|N|O>_<a|b>_bac\n")
cat("  ps_18 / ps_19 / ps_20  (cleaned UNRAREFIED per-year, DESeq2 DA input)\n")
cat("  NO combined physeq_rare_bac — years never merged.\n")

# ============================================================
# STAGE 2 — FUNGAL PIPELINE (ITS, CONSTAX + LCA taxonomy)
# PER-YEAR ARCHITECTURE — every step spelled out, once per year.

# WORKING DIRECTORY
# =========================================================
setwd("~/Desktop/Paper/NEW_R_Code/fungi_working")

# PER-YEAR RAREFACTION DEPTH

rare_depth_fun <- c("2018" = 7500, "2019" = 7500, "2020" = 7500)

# within-source collision guard (used ONLY for 2019 multi-source merge)
prefix_asv <- function(ps, tag) {
  taxa_names(ps) <- paste0(taxa_names(ps), "_", tag)   # renames otu + tax + refseq together
  ps
}

# RAREFACTION HELPER

manual_rarefy_fun <- function(physeq, depth) {
  set.seed(123)
  otu <- as(otu_table(physeq), "matrix")
  taxa_are_rows_val <- taxa_are_rows(otu_table(physeq))
  if (taxa_are_rows_val) otu <- t(otu)
  keep_samples <- rowSums(otu) >= depth
  otu_subset   <- otu[keep_samples, , drop = FALSE]
  rarefied     <- t(vegan::rrarefy(otu_subset, depth))
  otab <- if (taxa_are_rows_val) otu_table(rarefied, taxa_are_rows = TRUE) else otu_table(t(rarefied), taxa_are_rows = FALSE)
  ps_out <- phyloseq(otab, tax_table(physeq),
                     sample_data(physeq)[rownames(otu_subset), ], refseq(physeq))
  prune_taxa(taxa_sums(ps_out) > 0, ps_out)
}

# RAREFACTION DIAGNOSTIC PLOT HELPERS

RareStats_fun <- function(ps) {
  asv  <- as.data.frame(ps@otu_table)
  meta <- as.data.frame(as.matrix(ps@sam_data))
  findoutlier <- function(x) x < quantile(x,0.25)-1.5*IQR(x) | x > quantile(x,0.75)+1.5*IQR(x)
  asv %>%
    rownames_to_column("ASV_ID") %>%
    pivot_longer(-ASV_ID, names_to = "Sample_ID", values_to = "Seq_No") %>%
    group_by(Sample_ID) %>%
    summarize(Read_No = sum(Seq_No), Singleton_No = sum(Seq_No == 1),
              Goods = 100 * (1 - Singleton_No / Read_No)) %>%
    mutate(outlier = ifelse(findoutlier(log10(Read_No)), Read_No, NA)) %>%
    ungroup() %>%
    left_join(meta %>% rownames_to_column("Sample_ID"), by = "Sample_ID")
}

PloRareStats_fun <- function(df, year_label = "") {
  require(ggrepel)
  p <- ggarrange(
    df %>% ggplot(aes(x = Read_No)) + geom_histogram(binwidth = 5000, fill = "firebrick", color = "firebrick") + labs(title = paste(year_label, "Histogram")),
    df %>% ggplot(aes(x = Read_No)) + geom_histogram(binwidth = 1000, fill = "firebrick", color = "firebrick") + coord_cartesian(xlim = c(0,10000)) + labs(title = "Histogram Zoom"),
    df %>% ggplot(aes(x = Read_No, y = Goods)) + geom_point(shape = 1, color = "firebrick") + labs(title = "Good's Coverage"),
    df %>% ggplot(aes(x = 1, y = Read_No)) + geom_jitter(shape = 1, color = "firebrick") + scale_y_log10() + labs(title = "Log10 jitter"),
    df %>% ggplot(aes(x = 1, y = Read_No)) + geom_boxplot(color = "firebrick") + scale_y_log10() +
      geom_text_repel(data = dplyr::filter(df, !is.na(outlier)), mapping = aes(x = 1, y = Read_No, label = outlier), max.overlaps = 100, size = 3) + labs(title = "Log10 boxplot"),
    df %>% arrange(Read_No) %>% ggplot(aes(x = 1:nrow(.), y = Read_No)) + geom_bar(stat = "identity", color = "firebrick") + labs(title = "Ranked"),
    ncol = 3, nrow = 2, align = "hv", labels = c("A","B","C","D","E","F"))
  print(p); p
}


# PHYLUM CONSOLIDATION HELPER

consolidate_phyla <- function(ps) {
  phyla_mapping <- list(
    Mucoromycota    = c("Glomeromycota","Mortierellomycota","Diversisporales","Entrophosporales","Calcarisporiellomycota"),
    Basidiomycota   = c("Tremellomycetes","Agaricomycetes","Thelephorales"),
    Ascomycota      = c("Dothideomycetes","Pleosporales","Archaeorhizomycetes"),
    Zoopagomycota   = c("Kickxellomycota","Zoopagales","Entomophthoromycota","Zoopagomycetes","Basidiobolomycota"),
    Chytridiomycota = c("Monoblepharomycota","Olpidiomycota"),
    Cryptomycota    = c("Rozellomycota"))
  tax_mat <- tax_table(ps)
  for (target in names(phyla_mapping))
    for (source in phyla_mapping[[target]])
      tax_mat[tax_mat[,"Phylum"] == source, "Phylum"] <- target
  low_count <- c("Aphelidiomycota","Sanchytriomycota","Neocallimastigomycota","Endogonomycetes","Entorrhizomycota")
  tax_mat[tax_mat[,"Phylum"] %in% low_count, "Phylum"] <- ""
  tax_table(ps) <- tax_mat
  ps
}


# CONSTAX taxonomy loader for the 200bp/180bp SUPPLEMENTS

load_constax_consensus <- function(tax_file) {
  tx <- read.delim(tax_file, stringsAsFactors = FALSE, row.names = 1)
  data.frame(
    Kingdom = tx$Kingdom_Consensus, Phylum = tx$Phylum_Consensus,
    Class   = tx$Class_Consensus,   Order  = tx$Order_Consensus,
    Family  = tx$Family_Consensus,  Genus  = tx$Genus_Consensus,
    Species = tx$Species_Consensus, row.names = rownames(tx),
    stringsAsFactors = FALSE)
}


# ###  STEP 1 — LOAD MAIN RUN (all years) + LCA FILL     ###

cat("\n=== Loading Main Fungal Data (180bp, all years) ===\n")
ASV_fungi   <- read.delim("Fungi_ASVtable_180bp.txt", row.names = 1)
Meta_fungi  <- read.delim("fungi_metadata.txt", row.names = 1)
Tax_fungi   <- read.delim("fungi_constax_taxonomy_180.txt", header = TRUE, row.names = 1)
fasta_fungi <- readDNAStringSet("asv_180bp.fasta", format = "fasta",
                                seek.first.rec = TRUE, use.names = TRUE)

ps_main <- phyloseq(otu_table(ASV_fungi, taxa_are_rows = TRUE),
                    sample_data(Meta_fungi),
                    tax_table(as.matrix(Tax_fungi)),
                    fasta_fungi)
ps_main <- subset_samples(ps_main, is.na(Origin) | Origin != "Node")
ps_main <- subset_samples(ps_main, Experiment %in% c("obj_1"))
ps_main <- subset_samples(ps_main, Year != "2021")
ps_main <- prune_taxa(taxa_sums(ps_main) > 0, ps_main)
cat("Main loaded:", nsamples(ps_main), "samples,", ntaxa(ps_main), "ASVs\n")

# drop main-set copies of the 2019 C2 roots (preferred from the 200bp run)
cat("\n=== Removing main-set 2019 C2 root copies (200bp preferred) ===\n")
meta_200bp_check <- read.delim("2019_root_meta.txt", stringsAsFactors = FALSE)
samples_in_200bp <- meta_200bp_check$SampleID[
  meta_200bp_check$Growth_stage == "C2" &
    meta_200bp_check$Sample_or_Control == "True_Sample" &
    meta_200bp_check$Management %in% c("Conventional", "No-Till", "Organic") &
    meta_200bp_check$Origin == "root"]
cat("Removing", length(samples_in_200bp), "main-set samples coming from 200bp\n")
ps_main <- prune_samples(!sample_names(ps_main) %in% samples_in_200bp, ps_main)
ps_main <- prune_taxa(taxa_sums(ps_main) > 0, ps_main)
cat("Main after removal:", nsamples(ps_main), "samples\n")

# ----- LCA GAP-FILL (main only, ORIGINAL Zotu names, pre-prefix) -----
cat("\n=== LCA Gap-Fill (main, pre-prefix) ===\n")
tax_initial   <- as.data.frame(tax_table(ps_main))
uncult_before <- grepl("uncultured|unidentified|metagenome", tax_initial$Genus, ignore.case = TRUE)
empty_before  <- is.na(tax_initial$Genus) | tax_initial$Genus == ""
cat("  Total ASVs:", nrow(tax_initial),
    "| Well-ID:", sum(!uncult_before & !empty_before, na.rm = TRUE),
    "| Uncult:", sum(uncult_before, na.rm = TRUE),
    "| Empty:", sum(empty_before, na.rm = TRUE), "\n")

lca_csv <- "updating_tax/lca_majority70_toGenus.csv"
if (file.exists(lca_csv)) {
  lca_df  <- read.csv(lca_csv, stringsAsFactors = FALSE, check.names = FALSE)
  ranks   <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  missing <- setdiff(c("ASV", ranks), names(lca_df))
  if (length(missing) == 0) {
    tax <- as.data.frame(tax_table(ps_main), stringsAsFactors = FALSE)
    tax[] <- lapply(tax, function(x){ x <- as.character(x); x[is.na(x)] <- ""; trimws(x) })
    filled <- 0L
    for (i in seq_len(nrow(lca_df))) {
      asv <- lca_df$ASV[i]
      if (!asv %in% rownames(tax)) next
      changed <- FALSE
      for (rk in ranks) {
        if (rk == "Species") next
        val <- lca_df[[rk]][i]
        if (!is.null(val) && nzchar(val) && tax[asv, rk] == "") { tax[asv, rk] <- val; changed <- TRUE }
      }
      if (changed) filled <- filled + 1L
    }
    cat("LCA filled", filled, "additional ASVs\n")
    tax_table(ps_main) <- as.matrix(tax)
    uncult_final <- grepl("uncultured|unidentified|metagenome", tax$Genus, ignore.case = TRUE)
    empty_final  <- is.na(tax$Genus) | tax$Genus == ""
    cat("After LCA — Well-ID:", sum(!uncult_final & !empty_final),
        "| Uncult:", sum(uncult_final), "| Empty:", sum(empty_final), "\n")
  } else cat("LCA file missing columns:", paste(missing, collapse=", "), "\n")
} else cat("No LCA file found — skipping fill\n")

# split main by year (one run = consistent id space, collision-free split)
cat("\n=== Splitting main run by year ===\n")
sm <- as.data.frame(sample_data(ps_main))
main_2018 <- prune_samples(rownames(sm)[sm$Year == "2018"], ps_main)
main_2019 <- prune_samples(rownames(sm)[sm$Year == "2019"], ps_main)
main_2020 <- prune_samples(rownames(sm)[sm$Year == "2020"], ps_main)
main_2018 <- prune_taxa(taxa_sums(main_2018) > 0, main_2018)
main_2019 <- prune_taxa(taxa_sums(main_2019) > 0, main_2019)
main_2020 <- prune_taxa(taxa_sums(main_2020) > 0, main_2020)
cat("  main 2018:", nsamples(main_2018), "| 2019:", nsamples(main_2019),
    "| 2020:", nsamples(main_2020), "samples\n")

# ###  PER-YEAR PROCESSOR ###

process_year_fun <- function(ps, yr) {
  cat("\n##### PROCESSING", yr, "(fungi) #####\n")
  
  # --- STEP 4: metadata ---
  sd <- as.data.frame(sample_data(ps))
  
  remove_ctrl <- sd$Year == yr & sd$Growth_stage %in% c("C1","CO","NE","PO")
  if (sum(remove_ctrl, na.rm = TRUE) > 0) {
    cat("  removing control-stage samples:", sum(remove_ctrl, na.rm = TRUE), "\n")
    ps <- prune_samples(!remove_ctrl, ps); ps <- prune_taxa(taxa_sums(ps) > 0, ps)
    sd <- as.data.frame(sample_data(ps))
  }
  
  sd$Growth_Stage_Description <- gsub("Infloresence", "Inflorescence",
                                      sd$Growth_Stage_Description, ignore.case = TRUE)
  
  # 2019 rebuilds from C-codes; 2018/2020 already carry correct descriptions
  if (yr == "2019") {
    is_y <- !is.na(sd$Year) & sd$Year == "2019"
    sd$Growth_Stage_Description[is_y] <- dplyr::recode(
      sd$Growth_stage[is_y],
      "C2" = "Vegetative", "C3" = "Inflorescence", "C4" = "Reproductive",
      .default = sd$Growth_Stage_Description[is_y])
    cat("  2019 growth stages rebuilt (C2/C3/C4)\n")
  }
  sample_data(ps) <- sd
  
  sd <- as.data.frame(sample_data(ps))
  keep <- !(sd$Management %in% c("Control","Low-Input","Other")) &
    !(sd$Growth_Stage_Description %in% c("CONTROL","Control_or_Check","N/A"))
  ps <- prune_samples(keep, ps); ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  
  sd <- as.data.frame(sample_data(ps))
  sd$Compartment <- gsub("-", "", tolower(as.character(sd$Compartment)))
  sd$Management  <- factor(sd$Management, levels = mgmt_order)
  sd$Year        <- factor(sd$Year, levels = c("2018","2019","2020"))
  sample_data(ps) <- sd
  cat("  after metadata cleanup:", nsamples(ps), "samples\n")
  print(table(Compartment = sd$Compartment, Stage = sd$Growth_Stage_Description))
  
  # --- STEP 5: remove non-fungal taxa + consolidate phyla (per year) ---
  cat("  non-fungal removal + phylum consolidation\n")
  cat("    ASVs before:", ntaxa(ps), "\n")
  ps <- subset_taxa(ps,
                    !Phylum %in% c("Chloroplast","Mitochondria") &
                      !Kingdom %in% c("Anthophyta","Alveolata","Ichthyosporia",
                                      "Protista","Metazoa","Rhizaria","Viridiplantae"))
  ps <- consolidate_phyla(ps)
  cat("    ASVs after:", ntaxa(ps), "\n")
  
  # --- STEP 6: decontam + abundance filter (WITHIN year) ---
  sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control_Sample"
  if (any(sample_data(ps)$is.neg, na.rm = TRUE)) {
    cont <- isContaminant(ps, method = "prevalence", neg = "is.neg")
    cat("  contaminants:", sum(cont$contaminant), "\n")
    ps <- prune_taxa(!cont$contaminant, ps)
  } else cat("  no negative controls — skipping decontam\n")
  ps <- subset_samples(ps, Sample_or_Control == "True_Sample")
  cat("  after removing controls:", nsamples(ps), "samples\n")
  
  ps <- prune_samples(sample_sums(ps) >= 1000, ps)
  cat("  after <1000-read drop:", nsamples(ps), "samples\n")
  m <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) m <- t(m)
  keep_asv <- (rowSums(m >= 0) >= 3) & (rowSums(m) >= 10)
  cat("  abundance filter prev>=3 & reads>=10 | kept:", sum(keep_asv),
      "removed:", sum(!keep_asv), sprintf("(%.1f%%)\n", 100*mean(!keep_asv)))
  ps <- prune_taxa(keep_asv, ps)
  ps <- subset_samples(ps,
                       !Management %in% c("Control","Low-Input") & !Growth_Stage_Description %in% c("CONTROL"))
  
  # --- STEP 7: known outliers ---
  outlier_samples <- c()
  pres <- outlier_samples[outlier_samples %in% sample_names(ps)]
  if (length(pres) > 0) {
    ps <- prune_samples(!sample_names(ps) %in% pres, ps)
    cat("  removed known outliers:", paste(pres, collapse=", "), "\n")
  } else cat("  no known outliers\n")
  
  # --- STEP 8: compositional outliers (SF > 10), within year ---
  sample_data(ps)$Management <- factor(sample_data(ps)$Management, levels = mgmt_order)
  sf <- tryCatch({
    dd <- phyloseq_to_deseq2(ps, ~ Management)
    sizeFactors(estimateSizeFactors(dd, type = "poscounts"))
  }, error = function(e) { cat("  SF calc skipped:", conditionMessage(e), "\n"); NULL })
  if (!is.null(sf)) {
    bad <- names(sf)[sf > 10]
    cat("  compositional outliers (SF>10):", length(bad), "\n")
    if (length(bad) > 0) {
      md <- data.frame(sample_data(ps))
      info <- md[bad, c("Year","Management","Compartment","Growth_Stage_Description")]
      info$SizeFactor <- round(sf[bad], 2); print(info)
      ps <- prune_samples(!sample_names(ps) %in% bad, ps)
      ps <- prune_taxa(taxa_sums(ps) > 0, ps)
    }
  }
  
  # --- STEP 9: rarefy at this year's depth ---
  d <- rare_depth_fun[[yr]]
  cat("  rarefaction diagnostics:\n")
  PloRareStats_fun(RareStats_fun(ps), yr)
  ps <- prune_samples(sample_sums(ps) >= d, ps)
  
  assign(paste0("ps_unrare_", yr), ps, envir = .GlobalEnv)
  cat("  saved ps_unrare_", yr, " (", nsamples(ps), " samples, ", ntaxa(ps), " ASVs)\n", sep = "")
  
  ps_rare <- manual_rarefy_fun(ps, depth = d)
  sample_data(ps_rare)$Management <- factor(sample_data(ps_rare)$Management, levels = mgmt_order)
  sample_data(ps_rare)$Year       <- factor(sample_data(ps_rare)$Year, levels = c("2018","2019","2020"))
  dep <- sample_sums(ps_rare)
  cat("  rarefied to", d, "-> all==", d, ":", all(dep == d),
      "|", nsamples(ps_rare), "samples,", ntaxa(ps_rare), "ASVs\n")
  ps_rare
}


# ###  YEAR 2018 — SOYBEAN — main run only, no supplement ###

physeq_rare_2018 <- process_year_fun(main_2018, "2018")

# ###  YEAR 2020 — MAIZE — main run only, no supplement   ###

physeq_rare_2020 <- process_year_fun(main_2020, "2020")

# ###  YEAR 2019 — WHEAT — main + roots + stemleaf        ###

# ----- load 2019 C2 root supplement (200bp) -----
cat("\n=== [2019] Loading C2 Root Supplement (200bp) ===\n")
root_otu_file <- "otutable_UNOISE_200bp (1).txt"
root_tax_file <- "asv_200bp_taxonomy.txt"
root_fasta    <- "asv_200bp.fasta"
ps_19_roots <- NULL
if (all(file.exists(c("2019_root_meta.txt", root_otu_file, root_tax_file, root_fasta)))) {
  meta_2019_root <- read.delim("2019_root_meta.txt", stringsAsFactors = FALSE,
                               na.strings = c("NA","#N/A",""))
  meta_C2 <- meta_2019_root %>%
    filter(!is.na(SampleID) & SampleID != "#N/A", Growth_stage == "C2",
           Sample_or_Control == "True_Sample",
           Management %in% c("Conventional","No-Till","Organic"), Origin == "root") %>%
    mutate(Compartment = "belowground", Year = as.character(Year),
           Experiment = "obj_1", is.neg = FALSE)
  rownames(meta_C2) <- meta_C2$SampleID
  
  otu_2019_root <- read.delim(root_otu_file, row.names = 1, check.names = FALSE)
  sk  <- intersect(meta_C2$SampleID, colnames(otu_2019_root))
  otu_C2 <- otu_2019_root[, sk, drop = FALSE]; otu_C2 <- otu_C2[rowSums(otu_C2) > 0, ]
  meta_C2 <- meta_C2[sk, ]
  tax_C2 <- load_constax_consensus(root_tax_file)
  ak <- intersect(rownames(otu_C2), rownames(tax_C2)); tax_C2 <- tax_C2[ak, ]; otu_C2 <- otu_C2[ak, ]
  fr <- readDNAStringSet(root_fasta, format = "fasta"); fr <- fr[names(fr) %in% rownames(otu_C2)]; fr <- fr[rownames(otu_C2)]
  ps_19_roots <- phyloseq(otu_table(as.matrix(otu_C2), taxa_are_rows = TRUE),
                          sample_data(meta_C2), tax_table(as.matrix(tax_C2)), fr)
  ps_19_roots <- subset_taxa(ps_19_roots, Kingdom == "Fungi")
  ps_19_roots <- prune_taxa(taxa_sums(ps_19_roots) > 0, ps_19_roots)
  cat("  2019 roots:", nsamples(ps_19_roots), "samples,", ntaxa(ps_19_roots), "ASVs\n")
} else cat("  no 2019 root files found\n")

# ----- load 2019 stem/leaf supplement (180bp) -----
cat("\n=== [2019] Loading Stem/Leaf Supplement (180bp) ===\n")
sl_otu_file <- "2019_Stemleafadd_otutable_UNOISE_180bp.txt"
sl_tax_file <- "2019_Stemleafadd_taxonomy_180bp_new.txt"
sl_fasta    <- "2019_Stemleafadd_asv_180bp.fasta"
if (!all(file.exists(c("leafstem_fun_map_wheat_2019_finished.txt", sl_otu_file, sl_tax_file, sl_fasta))))
  stop("Missing 2019 stem/leaf supplemental files. Cannot proceed.")

meta_2019_sl <- read.delim("leafstem_fun_map_wheat_2019_finished.txt",
                           stringsAsFactors = FALSE, na.strings = c("NA","#N/A",""))
meta_sl <- meta_2019_sl %>%
  filter(!is.na(SampleID) & SampleID != "#N/A", Sample_or_Control == "True_Sample",
         Management %in% c("Conventional","No-Till","Organic")) %>%
  mutate(Compartment = "aboveground", Year = as.character(Year),
         Experiment = "obj_1", is.neg = FALSE)
rownames(meta_sl) <- meta_sl$SampleID

otu_2019_sl <- read.delim(sl_otu_file, row.names = 1, check.names = FALSE)
sk2 <- intersect(meta_sl$SampleID, colnames(otu_2019_sl))
otu_sl <- otu_2019_sl[, sk2, drop = FALSE]; otu_sl <- otu_sl[rowSums(otu_sl) > 0, ]
meta_sl <- meta_sl[sk2, ]
tax_sl <- load_constax_consensus(sl_tax_file)
ak2 <- intersect(rownames(otu_sl), rownames(tax_sl)); tax_sl <- tax_sl[ak2, ]; otu_sl <- otu_sl[ak2, ]
fs <- readDNAStringSet(sl_fasta, format = "fasta"); fs <- fs[names(fs) %in% rownames(otu_sl)]; fs <- fs[rownames(otu_sl)]
ps_19_sl <- phyloseq(otu_table(as.matrix(otu_sl), taxa_are_rows = TRUE),
                     sample_data(meta_sl), tax_table(as.matrix(tax_sl)), fs)
ps_19_sl <- subset_taxa(ps_19_sl, Kingdom == "Fungi")
ps_19_sl <- prune_taxa(taxa_sums(ps_19_sl) > 0, ps_19_sl)
cat("  2019 stem/leaf:", nsamples(ps_19_sl), "samples,", ntaxa(ps_19_sl), "ASVs\n")

# ----- PREFIX-GUARDED MERGE of the three 2019 sources -----

cat("\n=== [2019] Prefixing + merging three sources ===\n")
main_2019 <- prefix_asv(main_2019,  "main")

align_cols <- function(a, b) {
  ca <- colnames(sample_data(a)); cb <- colnames(sample_data(b))
  if (length(setdiff(ca, cb))) { x <- as.data.frame(sample_data(b)); for (c in setdiff(ca, cb)) x[[c]] <- NA; sample_data(b) <- x }
  if (length(setdiff(cb, ca))) { x <- as.data.frame(sample_data(a)); for (c in setdiff(cb, ca)) x[[c]] <- NA; sample_data(a) <- x }
  list(a, b)
}

ps_2019 <- main_2019
if (!is.null(ps_19_roots)) {
  ps_19_roots <- prefix_asv(ps_19_roots, "2019root200")
  pr <- align_cols(ps_2019, ps_19_roots); ps_2019 <- merge_phyloseq(pr[[1]], pr[[2]])
}
ps_19_sl <- prefix_asv(ps_19_sl, "2019sl180")
# drop any duplicate sample ids before folding stemleaf in
ov <- intersect(sample_names(ps_2019), sample_names(ps_19_sl))
if (length(ov) > 0) { ps_2019 <- prune_samples(!sample_names(ps_2019) %in% ov, ps_2019); ps_2019 <- prune_taxa(taxa_sums(ps_2019) > 0, ps_2019) }
ps2 <- align_cols(ps_2019, ps_19_sl); ps_2019 <- merge_phyloseq(ps2[[1]], ps2[[2]])
cat("  2019 merged (main+roots+stemleaf):", nsamples(ps_2019), "samples,", ntaxa(ps_2019), "ASVs\n")

# within-2019 source-tag verification
.m <- as(otu_table(ps_2019), "matrix"); if (!taxa_are_rows(ps_2019)) .m <- t(.m)
.src <- ifelse(grepl("_main$", rownames(.m)), "main",
               ifelse(grepl("_2019root200$", rownames(.m)), "root",
                      ifelse(grepl("_2019sl180$", rownames(.m)), "stemleaf", "??")))
cat("  [2019] id source tags (want only main/root/stemleaf, no ??):\n"); print(table(.src)); rm(.m, .src)

# validate stem/leaf actually backfilled aboveground
sdp <- as.data.frame(sample_data(ps_2019))
n_ag <- sum(sdp$Compartment == "aboveground", na.rm = TRUE)
cat("  [2019] aboveground samples after merge:", n_ag, "\n")
if (nsamples(ps_19_sl) == 0) stop("Stem/leaf added 0 samples. Cannot proceed.")

# now process 2019 through the same per-year pipeline
physeq_rare_2019 <- process_year_fun(ps_2019, "2019")

# ###  STEP 10 — SUBSETTING (per year, no cross-year obj) ###

cat("\n=== Creating Subsets (per year) ===\n")
year_objs <- list("2018" = physeq_rare_2018, "2019" = physeq_rare_2019, "2020" = physeq_rare_2020)
compartments <- c("aboveground","belowground"); comp_suffix <- c("a","b")
managements  <- c("Conventional","No-Till","Organic"); m_suffix <- c("C","N","O")

for (yr in names(year_objs)) {
  psr_year <- year_objs[[yr]]
  assign(paste0("psr_", yr), psr_year)
  for (i in seq_along(compartments)) {
    psr_comp <- prune_samples(sample_names(psr_year)[as.character(sample_data(psr_year)$Compartment) == compartments[i]], psr_year)
    psr_comp <- prune_taxa(taxa_sums(psr_comp) > 0, psr_comp)
    assign(paste0("psr_", yr, "_", comp_suffix[i]), psr_comp)
    for (j in seq_along(managements)) {
      psr_mgmt <- prune_samples(sample_names(psr_comp)[as.character(sample_data(psr_comp)$Management) == managements[j]], psr_comp)
      psr_mgmt <- prune_taxa(taxa_sums(psr_mgmt) > 0, psr_mgmt)
      assign(paste0("psr_", yr, "_", m_suffix[j], "_", comp_suffix[i]), psr_mgmt)
    }
  }
  gc()
}

cat("\n===========================================================\n")
cat("  FUNGAL PIPELINE COMPLETE (PER-YEAR)\n")
cat("===========================================================\n")
for (yr in names(year_objs)) {
  ps <- year_objs[[yr]]
  cat(sprintf("\n[%s] %d samples, %d ASVs\n", yr, nsamples(ps), ntaxa(ps)))
  print(table(Year = sample_data(ps)$Year, Compartment = sample_data(ps)$Compartment))
  cat("  Growth stage:\n"); print(table(sample_data(ps)$Growth_Stage_Description))
}

cat("\n--- PROOF: 2019 id source tags ---\n")
t19 <- taxa_names(physeq_rare_2019)
cat("  _main:", sum(grepl("_main$", t19)),
    "| _2019root200:", sum(grepl("_2019root200$", t19)),
    "| _2019sl180:", sum(grepl("_2019sl180$", t19)),
    "| untagged:", sum(!grepl("_main$|_2019root200$|_2019sl180$", t19)), "\n")

cat("\n--- PROOF: refseq present in slices ---\n")
for (nm in c("psr_2018_a","psr_2019_a","psr_2020_a","psr_2018_b","psr_2019_b","psr_2020_b")) {
  ps <- get(nm)
  rs <- tryCatch(length(refseq(ps)), error = function(e) NA)
  cat(sprintf("  %-14s seqs: %s\n", nm, rs))
}

## PROOF 3 — every slice flat at its year's rarefaction depth.
cat("\n--- PROOF: per-slice depth flat at year target ---\n")
for (yr in names(year_objs)) {
  d <- rare_depth_fun[[yr]]
  for (suf in c("a","b")) {
    ps <- get(paste0("psr_", yr, "_", suf))
    if (nsamples(ps) > 0) {
      dep <- sample_sums(ps)
      cat(sprintf("  psr_%s_%s: n=%d depth all==%d: %s\n", yr, suf, nsamples(ps), d, all(dep == d)))
    } else cat(sprintf("  psr_%s_%s: EMPTY\n", yr, suf))
  }
}

cat("\n--- PROOF: ps_unrare_<yr> present for DESeq2 DA ---\n")
for (yr in names(year_objs)) {
  nm <- paste0("ps_unrare_", yr)
  if (exists(nm) && inherits(get(nm), "phyloseq")) {
    ps <- get(nm)
    cat(sprintf("  %-16s n=%d ASVs=%d | samples match rarefied: %s\n",
                nm, nsamples(ps), ntaxa(ps),
                setequal(sample_names(ps), sample_names(year_objs[[yr]]))))
  } else cat(sprintf("  %-16s MISSING\n", nm))
}

cat("\nObjects available:\n")
cat("  physeq_rare_2018 / _2019 / _2020   (per-year rarefied finals)\n")
cat("  ps_unrare_2018 / _2019 / _2020     (cleaned UNRAREFIED, DESeq2 DA input)\n")
cat("  psr_<yr>, psr_<yr>_<a|b>, psr_<yr>_<C|N|O>_<a|b>\n")
cat("  NO combined physeq_rare — main split by year, only 2019 merged sources.\n")

setwd("~/Desktop/FINAL CODE FOR PAPER")

# =============================================================================
#  SECTION 2
#  Rarefaction curves
#
#  Supplemental Figure 1.

# ====== ENSURE per-year UNRAREFIED objects (curves are pre-rarefaction) ======
if (!exists("ps_unrare_2018_bac") && exists("ps_18")) ps_unrare_2018_bac <- ps_18
if (!exists("ps_unrare_2019_bac") && exists("ps_19")) ps_unrare_2019_bac <- ps_19
if (!exists("ps_unrare_2020_bac") && exists("ps_20")) ps_unrare_2020_bac <- ps_20
if (!all(sapply(c("ps_unrare_2018","ps_unrare_2019","ps_unrare_2020"), exists))) {
  if (exists("process_year_fun") && all(sapply(c("main_2018","main_2020","ps_2019"), exists))) {
    cat("Rebuilding fungal ps_unrare_* ...\n")
    body(process_year_fun)[[length(body(process_year_fun))]] <- quote(list(rare = ps_rare, unrare = ps))
    .t <- process_year_fun(main_2018, "2018"); physeq_rare_2018 <- .t$rare; ps_unrare_2018 <- .t$unrare
    .t <- process_year_fun(main_2020, "2020"); physeq_rare_2020 <- .t$rare; ps_unrare_2020 <- .t$unrare
    .t <- process_year_fun(ps_2019,   "2019"); physeq_rare_2019 <- .t$rare; ps_unrare_2019 <- .t$unrare
  } else stop("Fungal ps_unrare_* missing; source the fungal stage first.")
}
.need <- c("ps_unrare_2018","ps_unrare_2019","ps_unrare_2020",
           "ps_unrare_2018_bac","ps_unrare_2019_bac","ps_unrare_2020_bac")
.miss <- .need[!sapply(.need, exists)]
if (length(.miss) > 0) stop("Still missing: ", paste(.miss, collapse=", "))

# rarefaction depths (dashed reference line per kingdom)
rare_depth_fun <- 7500
rare_depth_bac <- 5000

# ====== plot rarefaction for one already-per-year object (no per-panel legend) ======
plot_rarefaction <- function(physeq_obj, year, kingdom_name, color_palette, depth_line = NA) {
  cat("Plotting", kingdom_name, year, "\n")
  phy_year <- prune_taxa(taxa_sums(physeq_obj) > 0, physeq_obj)
  otu <- as(otu_table(phy_year), "matrix"); if (taxa_are_rows(phy_year)) otu <- t(otu)
  metadata <- as.data.frame(sample_data(phy_year))
  curve_colors <- color_palette[as.character(metadata$Compartment)]
  rarecurve(otu, col = curve_colors, label = FALSE, step = 50,
            main = paste(kingdom_name, "Rarefaction -", year),
            xlab = "Number of DNA Reads", ylab = "Number of ASVs")
  if (!is.na(depth_line)) {
    abline(v = depth_line, lty = 2, col = "grey40")
    text(depth_line, par("usr")[4], paste0("rarefied to ", depth_line),
         pos = 4, offset = 0.3, cex = 0.7, col = "grey40", xpd = TRUE)
  }
}

# ====== LAYOUT: top legend strip + 2x3 panel grid ======
layout(matrix(c(1, 1, 1,
                2, 3, 4,
                5, 6, 7), nrow = 3, byrow = TRUE),
       heights = c(0.18, 1, 1))

# --- shared legend across the top ---
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", horiz = TRUE, bty = "n", cex = 1.1, lwd = 4, seg.len = 2.5,
       legend = c("Fungi aboveground", "Fungi belowground",
                  "Bacteria aboveground", "Bacteria belowground"),
       col    = c(fungi_colors["aboveground"], fungi_colors["belowground"],
                  bacteria_colors["aboveground"], bacteria_colors["belowground"]))

# --- panels ---
par(mar = c(4, 4, 3, 1))
plot_rarefaction(ps_unrare_2018,     "2018", "Fungi",    fungi_colors,    rare_depth_fun)
plot_rarefaction(ps_unrare_2019,     "2019", "Fungi",    fungi_colors,    rare_depth_fun)
plot_rarefaction(ps_unrare_2020,     "2020", "Fungi",    fungi_colors,    rare_depth_fun)
plot_rarefaction(ps_unrare_2018_bac, "2018", "Bacteria", bacteria_colors, rare_depth_bac)
plot_rarefaction(ps_unrare_2019_bac, "2019", "Bacteria", bacteria_colors, rare_depth_bac)
plot_rarefaction(ps_unrare_2020_bac, "2020", "Bacteria", bacteria_colors, rare_depth_bac)

layout(1); par(mfrow = c(1, 1))
cat("\nRarefaction plots complete.\n")

# =============================================================================
#  SECTION 3
#  Community composition tables
#
#  Supplemental Table 3

# ====== ENSURE per-year RAREFIED finals exist (composition uses rarefied) ======
if (!all(sapply(c("physeq_rare_2018","physeq_rare_2019","physeq_rare_2020"), exists))) {
  if (exists("process_year_fun") && all(sapply(c("main_2018","main_2020","ps_2019"), exists))) {
    cat("Rebuilding fungal per-year finals ...\n")
    body(process_year_fun)[[length(body(process_year_fun))]] <- quote(list(rare = ps_rare, unrare = ps))
    .t <- process_year_fun(main_2018, "2018"); physeq_rare_2018 <- .t$rare
    .t <- process_year_fun(main_2020, "2020"); physeq_rare_2020 <- .t$rare
    .t <- process_year_fun(ps_2019,   "2019"); physeq_rare_2019 <- .t$rare
  } else stop("Fungal per-year finals missing; source the fungal stage first.")
}
.need <- c("physeq_rare_2018","physeq_rare_2019","physeq_rare_2020",
           "physeq_rare_bac_2018","physeq_rare_bac_2019","physeq_rare_bac_2020")
.miss <- .need[!sapply(.need, exists)]
if (length(.miss) > 0) stop("Missing: ", paste(.miss, collapse=", "))

# -------------------------------------------------------------------
# Helper: long-format relative abundance from ONE object, optional ASV filter
# aggregate to taxon with rowsum, then per-sample proportions (numeric-safe)
# -------------------------------------------------------------------
build_relabund <- function(ps, rank, keep_asvs = NULL) {
  otu <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) otu <- t(otu)
  storage.mode(otu) <- "double"
  tax <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
  if (!is.null(keep_asvs)) {
    keep_asvs <- intersect(keep_asvs, rownames(otu))
    otu <- otu[keep_asvs, , drop = FALSE]; tax <- tax[keep_asvs, , drop = FALSE]
  }
  tax_vec <- tax[rownames(otu), rank]
  tax_vec[is.na(tax_vec) | tax_vec == ""] <- "Unclassified"
  
  agg  <- rowsum(otu, group = tax_vec)
  prop <- sweep(agg, 2, colSums(agg), "/")
  prop[!is.finite(prop)] <- 0
  
  long <- as.data.table(as.table(prop))
  setnames(long, c("Taxon", "Sample", "Abundance"))
  long[, Sample := as.character(Sample)]
  
  meta <- data.table(as.data.frame(sample_data(ps)), keep.rownames = "Sample")
  meta[, Sample := as.character(Sample)]
  merge(long, meta, by = "Sample")
}

build_relabund_years <- function(ps_list, rank, keep_fun = NULL) {
  rbindlist(lapply(ps_list, function(ps) {
    ka <- if (is.null(keep_fun)) NULL else keep_fun(ps)
    build_relabund(ps, rank, keep_asvs = ka)
  }), use.names = TRUE, fill = TRUE)
}

# -------------------------------------------------------------------
# Helper: summarize by Year × Compartment × Management (mean ± SD)
# -------------------------------------------------------------------
summarize_relabund <- function(long_dt, top_n = NULL, min_mean = NULL) {
  summ <- long_dt[, .(mean_pct = round(mean(Abundance)*100, 2),
                      sd_pct   = round(sd(Abundance)*100, 2)),
                  by = .(Year, Compartment, Management, Taxon)]
  overall <- long_dt[, .(overall_mean_pct = round(mean(Abundance)*100, 2)),
                     by = Taxon][order(-overall_mean_pct)]
  if (!is.null(min_mean)) { keep <- overall[overall_mean_pct >= min_mean, Taxon]
  summ <- summ[Taxon %in% keep]; overall <- overall[Taxon %in% keep] }
  if (!is.null(top_n))   { keep <- overall[1:min(top_n, .N), Taxon]
  summ <- summ[Taxon %in% keep]; overall <- overall[Taxon %in% keep] }
  summ[, value := paste0(formatC(mean_pct, format = "f", digits = 2),
                         " ± ",
                         formatC(sd_pct,   format = "f", digits = 2))]
  summ[, condition := paste(Year, Compartment, Management, sep = "_")]
  wide <- dcast(summ, Taxon ~ condition, value.var = "value", fill = "0.00 ± 0.00")
  wide <- merge(wide, overall, by = "Taxon"); setorder(wide, -overall_mean_pct)
  list(wide = wide, overall = overall)
}

fun_objs <- list(physeq_rare_2018, physeq_rare_2019, physeq_rare_2020)
bac_objs <- list(physeq_rare_bac_2018, physeq_rare_bac_2019, physeq_rare_bac_2020)

# -------------------------------------------------------------------
# FUNGAL TABLES
# -------------------------------------------------------------------
cat("Building fungal phylum table...\n")
fun_phylum_long <- build_relabund_years(fun_objs, "Phylum")
fun_phylum_tbl  <- summarize_relabund(fun_phylum_long, top_n = 10)

cat("Building fungal genus table...\n")
fun_genus_long <- build_relabund_years(fun_objs, "Genus")
fun_genus_long <- fun_genus_long[!grepl("^uncultured|^unidentified|^metagenome|^Unclassified",
                                        Taxon, ignore.case = TRUE)]
fun_genus_tbl <- summarize_relabund(fun_genus_long, min_mean = 1)

# -------------------------------------------------------------------
# BACTERIAL: per-object taxonomy filter, then relabund
# -------------------------------------------------------------------
bac_keep_fun <- function(ps) {
  tax_bac <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
  keep1 <- !is.na(tax_bac$Kingdom) & tax_bac$Kingdom == "Bacteria"
  keep2 <- is.na(tax_bac$Order)  | tax_bac$Order  != "Chloroplast"
  keep3 <- is.na(tax_bac$Family) | tax_bac$Family != "Mitochondria"
  euk_patterns <- c("Prosopis","Phytophthora","Pythium","Peronospora","Albugo","Plasmopara","Phytomonas")
  euk_logical <- Reduce(`|`, lapply(euk_patterns, function(p)
    grepl(p, tax_bac$Genus, ignore.case = TRUE) | grepl(p, tax_bac$Species, ignore.case = TRUE)))
  euk_logical[is.na(euk_logical)] <- FALSE
  rownames(tax_bac)[keep1 & keep2 & keep3 & !euk_logical]
}

cat("\nBuilding bacterial phylum table...\n")
bac_phylum_long <- build_relabund_years(bac_objs, "Phylum", keep_fun = bac_keep_fun)
bac_phylum_tbl  <- summarize_relabund(bac_phylum_long, top_n = 10)

cat("Building bacterial genus table...\n")
bac_genus_long <- build_relabund_years(bac_objs, "Genus", keep_fun = bac_keep_fun)
bac_genus_long <- bac_genus_long[!grepl("^uncultured|^unidentified|^metagenome|^Unclassified",
                                        Taxon, ignore.case = TRUE)]
bac_genus_tbl <- summarize_relabund(bac_genus_long, min_mean = 1)

# -------------------------------------------------------------------
# rankings
# -------------------------------------------------------------------
cat("\n=== TOP 10 FUNGAL PHYLA ===\n");                  print(fun_phylum_tbl$overall)
cat("\n=== FUNGAL GENERA >=1% ===\n");                   print(fun_genus_tbl$overall)
cat("\n=== TOP 10 BACTERIAL PHYLA (post-filter) ===\n"); print(bac_phylum_tbl$overall)
cat("\n=== BACTERIAL GENERA >=1% (post-filter) ===\n");  print(bac_genus_tbl$overall)

# -------------------------------------------------------------------
# COMBINE into single supplemental table
# -------------------------------------------------------------------
tag_table <- function(wide_tbl, kingdom, rank) {
  dt <- copy(wide_tbl); dt[, Kingdom := kingdom]; dt[, Rank := rank]
  setcolorder(dt, c("Kingdom","Rank","Taxon","overall_mean_pct")); dt
}
fun_phy <- tag_table(fun_phylum_tbl$wide, "Fungi",    "Phylum")
fun_gen <- tag_table(fun_genus_tbl$wide,  "Fungi",    "Genus")
bac_phy <- tag_table(bac_phylum_tbl$wide, "Bacteria", "Phylum")
bac_gen <- tag_table(bac_genus_tbl$wide,  "Bacteria", "Genus")

combined_supp_tbl <- rbindlist(list(fun_phy, fun_gen, bac_phy, bac_gen),
                               use.names = TRUE, fill = TRUE)
combined_supp_tbl[, Kingdom := factor(Kingdom, levels = c("Fungi","Bacteria"))]
combined_supp_tbl[, Rank    := factor(Rank,    levels = c("Phylum","Genus"))]
setorder(combined_supp_tbl, Kingdom, Rank, -overall_mean_pct)

cat("\n=== COMBINED TABLE — first 15 rows ===\n"); print(head(combined_supp_tbl, 15))
cat("\nDimensions:", nrow(combined_supp_tbl), "rows x", ncol(combined_supp_tbl), "cols\n")
cat("\nSection sizes:\n"); print(combined_supp_tbl[, .N, by = .(Kingdom, Rank)])

fwrite(combined_supp_tbl, "supp_table_community_composition.csv")
cat("\nExported: supp_table_community_composition.csv\n")

library(phyloseq); library(data.table); library(openxlsx)

# ---- ensure per-year rarefied finals ----
if (!all(sapply(c("physeq_rare_2018","physeq_rare_2019","physeq_rare_2020"), exists))) {
  if (exists("process_year_fun") && all(sapply(c("main_2018","main_2020","ps_2019"), exists))) {
    body(process_year_fun)[[length(body(process_year_fun))]] <- quote(list(rare = ps_rare, unrare = ps))
    .t <- process_year_fun(main_2018,"2018"); physeq_rare_2018 <- .t$rare
    .t <- process_year_fun(main_2020,"2020"); physeq_rare_2020 <- .t$rare
    .t <- process_year_fun(ps_2019,  "2019"); physeq_rare_2019 <- .t$rare
  } else stop("Fungal finals missing; source the fungal stage first.")
}
stopifnot(all(sapply(c("physeq_rare_bac_2018","physeq_rare_bac_2019","physeq_rare_bac_2020"), exists)))

CUTOFF <- 2
COMPS  <- c("aboveground","belowground"); COMP_TITLE <- c(aboveground="Aboveground", belowground="Belowground")
MGMTS  <- c("Conventional","No-Till","Organic"); MGMT_SHORT <- c(Conventional="Conv","No-Till"="NT",Organic="Org")
CROP   <- c("2018"="Soybean","2019"="Wheat","2020"="Maize")
BC <- "#BFBFBF"
NOT_NAMED <- function(x) grepl("^Unclassified$|^uncultured|^unidentified|^metagenome|Incertae", x, ignore.case=TRUE)

# per-sample genus proportions (each sample sums to 1; sample belongs to one compartment)
genus_long <- function(ps, keep_asvs = NULL) {
  otu <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) otu <- t(otu)
  storage.mode(otu) <- "double"
  tax <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
  if (!is.null(keep_asvs)) { ka <- intersect(keep_asvs, rownames(otu)); otu <- otu[ka,,drop=FALSE]; tax <- tax[ka,,drop=FALSE] }
  g <- tax[rownames(otu), "Genus"]; g[is.na(g) | g==""] <- "Unclassified"
  agg <- rowsum(otu, group = g); prop <- sweep(agg, 2, colSums(agg), "/"); prop[!is.finite(prop)] <- 0
  long <- as.data.table(as.table(prop)); setnames(long, c("Genus","Sample","Abundance"))
  long[, Sample := as.character(Sample)]
  meta <- data.table(as.data.frame(sample_data(ps)), keep.rownames="Sample"); meta[, Sample := as.character(Sample)]
  merge(long, meta, by="Sample")
}
bac_keep <- function(ps) {
  tx <- as.data.frame(tax_table(ps), stringsAsFactors=FALSE)
  k1 <- !is.na(tx$Kingdom) & tx$Kingdom=="Bacteria"
  k2 <- is.na(tx$Order)  | tx$Order  != "Chloroplast"
  k3 <- is.na(tx$Family) | tx$Family != "Mitochondria"
  euk <- Reduce(`|`, lapply(c("Prosopis","Phytophthora","Pythium","Peronospora","Albugo","Plasmopara","Phytomonas"),
                            function(p) grepl(p,tx$Genus,ignore.case=TRUE)|grepl(p,tx$Species,ignore.case=TRUE)))
  euk[is.na(euk)] <- FALSE
  rownames(tx)[k1 & k2 & k3 & !euk]
}
# mean % per Year x Compartment x Management (no prune_samples; filter the long table)
cell_means <- function(ps_list, keep_fun=NULL) {
  long <- rbindlist(lapply(ps_list, function(ps)
    genus_long(ps, keep_asvs = if (is.null(keep_fun)) NULL else keep_fun(ps))),
    use.names=TRUE, fill=TRUE)
  long[, Year := as.character(Year)]
  long[, Compartment := as.character(Compartment)]
  long <- long[Compartment %in% COMPS]
  long[, .(pct = round(mean(Abundance)*100, 2)), by=.(Genus, Year, Compartment, Management)]
}
shade_col <- function(pct) {
  if (is.na(pct) || pct<=0) return(NULL)
  fr <- min(pct/40, 1); r<-round(255+fr*(0-255)); g<-round(255+fr*(91-255)); b<-round(255+fr*(64-255))
  sprintf("#%02X%02X%02X", r, g, b)
}

build_kingdom_sheet <- function(wb, sheet, cm, kingdom) {
  addWorksheet(wb, sheet)
  hdr_yr   <- createStyle(fgFill="#2F4858", fontColour="#FFFFFF", textDecoration="bold", halign="left",  border="TopBottomLeftRight", borderColour=BC)
  hdr_comp <- createStyle(fgFill="#5B7A8C", fontColour="#FFFFFF", textDecoration="bold", halign="left",  border="TopBottomLeftRight", borderColour=BC)
  hdr_mgmt <- createStyle(fgFill="#9DB4C0", fontColour="#000000", textDecoration="bold", halign="center",border="TopBottomLeftRight", borderColour=BC)
  gstyle   <- createStyle(textDecoration="italic", halign="left",  border="TopBottomLeftRight", borderColour=BC)
  ustyle   <- createStyle(halign="left",  textDecoration="bold", fontColour="#7a3030", border="TopBottomLeftRight", borderColour=BC, fgFill="#F2E4E4")
  unum     <- createStyle(numFmt='0.0', halign="center", textDecoration="bold", fontColour="#7a3030", border="TopBottomLeftRight", borderColour=BC, fgFill="#F2E4E4")
  ostyle   <- createStyle(halign="left",  textDecoration="bold", border="TopBottomLeftRight", borderColour=BC, fgFill="#EEEEEE")
  onum     <- createStyle(numFmt='0.0', halign="center", textDecoration="bold", border="TopBottomLeftRight", borderColour=BC, fgFill="#EEEEEE")
  tnum     <- createStyle(numFmt='0.0', halign="center", textDecoration="bold", border="TopBottomLeftRight", borderColour=BC, fgFill="#D9D9D9")
  tstyle   <- createStyle(halign="left",  textDecoration="bold", border="TopBottomLeftRight", borderColour=BC, fgFill="#D9D9D9")
  num0     <- createStyle(numFmt='0.0;;""', halign="center", border="TopBottomLeftRight", borderColour=BC)
  title_st <- createStyle(textDecoration="bold", fontSize=12)
  
  writeData(wb, sheet, paste0(kingdom, " genera \u2014 mean relative abundance (%), aboveground and belowground separately"),
            startRow=1, startCol=1); addStyle(wb, sheet, title_st, rows=1, cols=1)
  
  r <- 3
  for (y in c("2018","2019","2020")) {
    writeData(wb, sheet, paste0(y, "  (", CROP[[y]], ")"), startRow=r, startCol=1)
    mergeCells(wb, sheet, cols=1:4, rows=r); addStyle(wb, sheet, hdr_yr, rows=r, cols=1:4, gridExpand=TRUE)
    r <- r+1
    
    for (cmp in COMPS) {
      sub <- cm[Year==y & Compartment==cmp]
      if (nrow(sub)==0) next
      keep <- sub[pct>=CUTOFF & !NOT_NAMED(Genus), unique(Genus)]
      rank <- sub[Genus %in% keep, .(s=sum(pct)), by=Genus][order(-s), Genus]
      
      writeData(wb, sheet, COMP_TITLE[[cmp]], startRow=r, startCol=1)
      mergeCells(wb, sheet, cols=1:4, rows=r); addStyle(wb, sheet, hdr_comp, rows=r, cols=1:4, gridExpand=TRUE)
      r <- r+1
      
      writeData(wb, sheet, "Genus", startRow=r, startCol=1); addStyle(wb, sheet, hdr_mgmt, rows=r, cols=1)
      cidx<-list(); col<-2
      for (m in MGMTS) { writeData(wb, sheet, MGMT_SHORT[[m]], startRow=r, startCol=col)
        addStyle(wb, sheet, hdr_mgmt, rows=r, cols=col); cidx[[m]]<-col; col<-col+1 }
      r<-r+1
      
      shown   <- setNames(numeric(length(MGMTS)), MGMTS)
      unclass <- setNames(numeric(length(MGMTS)), MGMTS)
      for (m in MGMTS) { uv <- sub[NOT_NAMED(Genus) & Management==m, sum(pct)]; unclass[m] <- if (length(uv)==0) 0 else uv }
      
      for (gn in rank) {
        writeData(wb, sheet, gn, startRow=r, startCol=1); addStyle(wb, sheet, gstyle, rows=r, cols=1)
        for (m in MGMTS) {
          v <- sub[Genus==gn & Management==m, pct]; v <- if (length(v)==0) 0 else v
          shown[m] <- shown[m] + v
          cc <- cidx[[m]]
          writeData(wb, sheet, if (v>0) v else NA, startRow=r, startCol=cc)
          st <- num0; hexc <- shade_col(v)
          if (!is.null(hexc)) st <- createStyle(numFmt='0.0;;""', halign="center",
                                                border="TopBottomLeftRight", borderColour=BC, fgFill=hexc)
          addStyle(wb, sheet, st, rows=r, cols=cc)
        }
        r<-r+1
      }
      
      writeData(wb, sheet, "Other named genera (<2%)", startRow=r, startCol=1); addStyle(wb, sheet, ostyle, rows=r, cols=1)
      for (m in MGMTS) { o <- round(100 - shown[m] - unclass[m],1); if (o<0) o<-0
      writeData(wb, sheet, o, startRow=r, startCol=cidx[[m]]); addStyle(wb, sheet, onum, rows=r, cols=cidx[[m]]) }
      r<-r+1
      writeData(wb, sheet, "Unclassified (no genus)", startRow=r, startCol=1); addStyle(wb, sheet, ustyle, rows=r, cols=1)
      for (m in MGMTS) { writeData(wb, sheet, round(unclass[m],1), startRow=r, startCol=cidx[[m]]); addStyle(wb, sheet, unum, rows=r, cols=cidx[[m]]) }
      r<-r+1
      writeData(wb, sheet, "Total", startRow=r, startCol=1); addStyle(wb, sheet, tstyle, rows=r, cols=1)
      for (m in MGMTS) { writeData(wb, sheet, 100, startRow=r, startCol=cidx[[m]]); addStyle(wb, sheet, tnum, rows=r, cols=cidx[[m]]) }
      r<-r+2
    }
    r<-r+1
  }
  
  setColWidths(wb, sheet, cols=1, widths=26)
  setColWidths(wb, sheet, cols=2:4, widths=9)
  writeData(wb, sheet, "Values = mean relative abundance (%) within each compartment. Aboveground and belowground are separate tables with their own genera (\u22652% in that compartment that year). 'Other named genera' = named genera <2%; 'Unclassified' = reads not resolved to genus. Each column sums to 100.",
            startRow=r, startCol=1)
  addStyle(wb, sheet, createStyle(fontSize=8, textDecoration="italic", fontColour="#808080"), rows=r, cols=1)
}

cm_fun <- cell_means(list(physeq_rare_2018, physeq_rare_2019, physeq_rare_2020))
cm_bac <- cell_means(list(physeq_rare_bac_2018, physeq_rare_bac_2019, physeq_rare_bac_2020), keep_fun=bac_keep)

wb <- createWorkbook()
build_kingdom_sheet(wb, "Fungi",    cm_fun, "Fungal")
build_kingdom_sheet(wb, "Bacteria", cm_bac, "Bacterial")
saveWorkbook(wb, "genus_relabund_byyear.xlsx", overwrite=TRUE)
cat("Saved genus_relabund_byyear.xlsx (sheets: Fungi, Bacteria)\n")

# =============================================================================
#  SECTION 4
#  Alpha diversity (Hill Numbers)
#
#  Supplemental Figure 2 and the Kruskal-Wallis table.

stopifnot(exists("physeq_rare_2018"), exists("physeq_rare_2019"), exists("physeq_rare_2020"),
          exists("physeq_rare_bac_2018"), exists("physeq_rare_bac_2019"), exists("physeq_rare_bac_2020"))

q_vals <- c(0, 1, 2)

build_hill_df <- function(ps_list) {
  bind_rows(lapply(ps_list, function(ps) {
    otu <- as(otu_table(ps), "matrix"); if (taxa_are_rows(ps)) otu <- t(otu)
    meta <- as(sample_data(ps), "data.frame"); meta$SampleID <- rownames(meta)
    bind_rows(lapply(q_vals, function(q)
      data.frame(SampleID = rownames(otu), Hill = hill_taxa(otu, q), q = q))) %>%
      left_join(meta, by = "SampleID")
  })) %>%
    mutate(
      Year = factor(Year),
      q = factor(q, labels = c("Richness (q=0)", "Shannon (q=1)", "Simpson (q=2)"))
    )
}

# FUNGI
hill_df_fungi <- build_hill_df(list(physeq_rare_2018, physeq_rare_2019, physeq_rare_2020))
hill_df_above_fungi <- hill_df_fungi %>% filter(Compartment == "aboveground")
hill_df_below_fungi <- hill_df_fungi %>% filter(Compartment == "belowground")

# BACTERIA
hill_df_bac <- build_hill_df(list(physeq_rare_bac_2018, physeq_rare_bac_2019, physeq_rare_bac_2020))
hill_df_above_bac <- hill_df_bac %>% filter(Compartment == "aboveground")
hill_df_below_bac <- hill_df_bac %>% filter(Compartment == "belowground")


# DATA CHECK — confirm frames match the per-year objects

cat("\n=== ALPHA DIVERSITY DATA CHECK ===\n")
for (king in c("fungi","bac")) {
  objs <- if (king=="bac") list("2018"=physeq_rare_bac_2018,"2019"=physeq_rare_bac_2019,"2020"=physeq_rare_bac_2020)
  else            list("2018"=physeq_rare_2018,"2019"=physeq_rare_2019,"2020"=physeq_rare_2020)
  df   <- if (king=="bac") hill_df_bac else hill_df_fungi
  n_obj <- sum(sapply(objs, nsamples)); n_df <- n_distinct(df$SampleID)
  cat(sprintf("%-7s total: objects=%d frame=%d match=%s\n", king, n_obj, n_df, n_obj==n_df))
  for (y in c("2018","2019","2020")) {
    no <- nsamples(objs[[y]]); nd <- df %>% filter(Year==y) %>% distinct(SampleID) %>% nrow()
    cat(sprintf("    %s: object=%d frame=%d match=%s\n", y, no, nd, no==nd))
  }
  cat(sprintf("    NA Hill=%d NA Compartment=%d NA Year=%d | above+below==total: %s\n",
              sum(is.na(df$Hill)), sum(is.na(df$Compartment)), sum(is.na(df$Year)),
              n_distinct(df$SampleID[df$Compartment=="aboveground"]) +
                n_distinct(df$SampleID[df$Compartment=="belowground"]) == n_distinct(df$SampleID)))
}


# HELPER FUNCTION: Year panel

make_year_plot <- function(data, compartment_label, panel_tag = NULL) {
  kw_labels <- data %>%
    group_by(q) %>%
    summarise(p = kruskal.test(Hill ~ Year)$p.value, .groups = "drop") %>%
    mutate(label = case_when(
      p < 0.0001 ~ "Kruskal-Wallis, p < 0.0001",
      TRUE ~ paste0("Kruskal-Wallis, p = ", signif(p, 2))))
  
  p <- ggplot(data, aes(x = Year, y = Hill, color = Year, fill = Year)) +
    geom_hline(yintercept = 0, color = "grey90", size = 0.3) +
    geom_violin(alpha = 0.3, scale = "width", width = 0.9, trim = FALSE, linewidth = 0.5) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.6, width = 0.1, alpha = 0, color = "black") +
    geom_point(position = position_jitter(width = 0.08), size = 1.2, alpha = 0.5, shape = 16) +
    geom_text(data = kw_labels, aes(x = -Inf, y = -Inf, label = label), inherit.aes = FALSE,
              vjust = -0.5, hjust = -0.05, size = 2.8, color = "grey10", fontface = "bold") +
    facet_wrap(~ q, scales = "free_y", nrow = 1) +
    scale_color_manual(values = year_colors, labels = c("2018 (Soy)", "2019 (Wheat)", "2020 (Maize)")) +
    scale_fill_manual(values = year_colors, labels = c("2018 (Soy)", "2019 (Wheat)", "2020 (Maize)")) +
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.15)), labels = comma) +
    labs(y = "Hill Number", x = NULL, title = compartment_label, color = "Year (Crop)", fill = "Year (Crop)") +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey95", size = 0.2),
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
      strip.text = element_text(face = "bold", size = 11, color = "grey20"),
      strip.background = element_rect(fill = "grey96", color = "grey80"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11, face = "bold", color = "black"),
      axis.text.y = element_text(size = 9, color = "grey30"),
      axis.title.y = element_text(size = 12, face = "bold", color = "grey20", margin = ggplot2::margin(r = 10)),
      axis.ticks = element_line(color = "grey80", size = 0.3),
      legend.position = "top", legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9), legend.key.size = unit(0.4, "cm"),
      panel.spacing = unit(0.8, "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "grey20", margin = ggplot2::margin(b = 2)),
      plot.margin = ggplot2::margin(t = 5, r = 10, b = 5, l = 10),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    coord_cartesian(clip = "off") +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
  
  year_comparisons <- list(c("2018", "2019"), c("2018", "2020"), c("2019", "2020"))
  p <- p + stat_compare_means(comparisons = year_comparisons, method = "wilcox.test",
                              p.adjust.method = "BH", label = "p.signif", hide.ns = TRUE,
                              size = 2.8, step.increase = 0.12, color = "grey30")
  
  if (!is.null(panel_tag)) {
    p <- p + labs(tag = paste0("(", panel_tag, ")")) +
      theme(plot.tag = element_text(face = "bold", size = 14), plot.tag.position = c(0, 1))
  }
  p
}


# HELPER FUNCTION: Management / Growth Stage panels

make_hill_plot_violin <- function(data, x_var, color_palette, panel_tag = NULL) {
  comparisons <- combn(levels(factor(data[[x_var]])), 2, simplify = FALSE)
  kw_labels <- data %>%
    group_by(q, Year) %>%
    summarise(p = kruskal.test(Hill ~ .data[[x_var]])$p.value, .groups = "drop") %>%
    mutate(label = case_when(
      p < 0.0001 ~ "Kruskal-Wallis, p < 0.0001",
      TRUE ~ paste0("Kruskal-Wallis, p = ", signif(p, 2))))
  
  p <- ggplot(data, aes(x = .data[[x_var]], y = Hill, color = .data[[x_var]], fill = .data[[x_var]])) +
    geom_hline(yintercept = 0, color = "grey90", size = 0.3) +
    geom_violin(alpha = 0.3, scale = "width", width = 0.9, linewidth = 0.5, trim = FALSE) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.6, width = 0.15, alpha = 0, color = "black") +
    geom_point(position = position_jitter(width = 0.08), size = 1.2, alpha = 0.5, shape = 16) +
    geom_text(data = kw_labels, aes(x = -Inf, y = -Inf, label = label), inherit.aes = FALSE,
              vjust = -0.5, hjust = -0.05, size = 2.8, color = "grey10", fontface = "bold") +
    facet_grid(q ~ Year, scales = "free_y", switch = "y") +
    scale_color_manual(values = color_palette) +
    scale_fill_manual(values = color_palette) +
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.15)), labels = comma) +
    labs(y = "Hill Number", x = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey95", size = 0.2),
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
      strip.text.x = element_text(face = "bold", size = 11, color = "grey20"),
      strip.text.y = element_text(face = "bold", size = 10, angle = 0, color = "grey20"),
      strip.background = element_rect(fill = "grey96", color = "grey80"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold", color = "black"),
      axis.text.y = element_text(size = 9, color = "grey30"),
      axis.title.y = element_text(size = 12, face = "bold", color = "grey20", margin = ggplot2::margin(r = 10)),
      axis.ticks = element_line(color = "grey80", size = 0.3),
      legend.position = "top", legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9), legend.key.size = unit(0.4, "cm"),
      panel.spacing.x = unit(0.5, "lines"), panel.spacing.y = unit(0.4, "lines"),
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 10, l = 10),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    coord_cartesian(clip = "off") +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
  
  p <- p + stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                              p.adjust.method = "BH", label = "p.signif", hide.ns = TRUE,
                              size = 3.2, step.increase = 0.12, tip.length = 0.02, color = "grey30")
  
  if (!is.null(panel_tag)) {
    p <- p + labs(tag = paste0("(", panel_tag, ")")) +
      theme(plot.tag = element_text(face = "bold", size = 14), plot.tag.position = c(0, 1))
  }
  p
}


# BUILD ALL 12 PANELS

p_year_above_fungi <- make_year_plot(hill_df_above_fungi, "Aboveground", "A")
p_year_below_fungi <- make_year_plot(hill_df_below_fungi, "Belowground", "B")
p_mgmt_above_fungi <- make_hill_plot_violin(hill_df_above_fungi, "Management", management_colors, "C")
p_mgmt_below_fungi <- make_hill_plot_violin(hill_df_below_fungi, "Management", management_colors, "D")
p_gs_above_fungi   <- make_hill_plot_violin(hill_df_above_fungi, "Growth_Stage_Description", growth_stage_colors, "E")
p_gs_below_fungi   <- make_hill_plot_violin(hill_df_below_fungi, "Growth_Stage_Description", growth_stage_colors, "F")

p_year_above_bac <- make_year_plot(hill_df_above_bac, "Aboveground", "G")
p_year_below_bac <- make_year_plot(hill_df_below_bac, "Belowground", "H")
p_mgmt_above_bac <- make_hill_plot_violin(hill_df_above_bac, "Management", management_colors, "I")
p_mgmt_below_bac <- make_hill_plot_violin(hill_df_below_bac, "Management", management_colors, "J")
p_gs_above_bac   <- make_hill_plot_violin(hill_df_above_bac, "Growth_Stage_Description", growth_stage_colors, "K")
p_gs_below_bac   <- make_hill_plot_violin(hill_df_below_bac, "Growth_Stage_Description", growth_stage_colors, "L")


# SECTION TITLES

fungi_title    <- textGrob("Fungi", gp = gpar(fontsize = 18, fontface = "bold", col = "grey20"),
                           just = "left", x = unit(0.02, "npc"))
bacteria_title <- textGrob("Bacteria", gp = gpar(fontsize = 18, fontface = "bold", col = "grey20"),
                           just = "left", x = unit(0.02, "npc"))


# COMBINE (stacked)

combined_fig <- wrap_elements(full = fungi_title) /
  (p_year_above_fungi | p_year_below_fungi) /
  (p_mgmt_above_fungi | p_mgmt_below_fungi) /
  (p_gs_above_fungi   | p_gs_below_fungi) /
  wrap_elements(full = bacteria_title) /
  (p_year_above_bac | p_year_below_bac) /
  (p_mgmt_above_bac | p_mgmt_below_bac) /
  (p_gs_above_bac   | p_gs_below_bac) +
  plot_layout(heights = c(0.3, 1, 2.5, 2.5, 0.3, 1, 2.5, 2.5))

print(combined_fig)


# ALTERNATIVE LAYOUT: side by side (Fungi left, Bacteria right)

fungi_combined <- wrap_elements(full = fungi_title) /
  (p_year_above_fungi | p_year_below_fungi) /
  (p_mgmt_above_fungi | p_mgmt_below_fungi) /
  (p_gs_above_fungi   | p_gs_below_fungi) +
  plot_layout(heights = c(0.3, 1, 2.5, 2.5))

bacteria_combined <- wrap_elements(full = bacteria_title) /
  (p_year_above_bac | p_year_below_bac) /
  (p_mgmt_above_bac | p_mgmt_below_bac) /
  (p_gs_above_bac   | p_gs_below_bac) +
  plot_layout(heights = c(0.3, 1, 2.5, 2.5))

combined_fig_sidebyside <- fungi_combined | bacteria_combined
print(combined_fig_sidebyside)


# CONSOLIDATE KRUSKAL-WALLIS RESULTS — one tidy table for the Results section

kw_table <- function(data, factor_var, kingdom, compartment, scope = "overall") {
  if (scope == "overall") {
    data %>%
      group_by(q) %>%
      summarise(
        chi_sq = kruskal.test(Hill ~ .data[[factor_var]])$statistic,
        df     = kruskal.test(Hill ~ .data[[factor_var]])$parameter,
        p      = kruskal.test(Hill ~ .data[[factor_var]])$p.value,
        n      = sum(!is.na(Hill)), .groups = "drop") %>%
      mutate(Kingdom = kingdom, Compartment = compartment, Factor = factor_var, Year = "All")
  } else {
    data %>%
      group_by(q, Year) %>%
      summarise(
        chi_sq = kruskal.test(Hill ~ .data[[factor_var]])$statistic,
        df     = kruskal.test(Hill ~ .data[[factor_var]])$parameter,
        p      = kruskal.test(Hill ~ .data[[factor_var]])$p.value,
        n      = sum(!is.na(Hill)), .groups = "drop") %>%
      mutate(Kingdom = kingdom, Compartment = compartment, Factor = factor_var, Year = as.character(Year))
  }
}

year_tests <- bind_rows(
  kw_table(hill_df_above_fungi, "Year", "Fungi",    "Aboveground", "overall"),
  kw_table(hill_df_below_fungi, "Year", "Fungi",    "Belowground", "overall"),
  kw_table(hill_df_above_bac,   "Year", "Bacteria", "Aboveground", "overall"),
  kw_table(hill_df_below_bac,   "Year", "Bacteria", "Belowground", "overall"))

mgmt_tests <- bind_rows(
  kw_table(hill_df_above_fungi, "Management", "Fungi",    "Aboveground", "by_year"),
  kw_table(hill_df_below_fungi, "Management", "Fungi",    "Belowground", "by_year"),
  kw_table(hill_df_above_bac,   "Management", "Bacteria", "Aboveground", "by_year"),
  kw_table(hill_df_below_bac,   "Management", "Bacteria", "Belowground", "by_year"))

gs_tests <- bind_rows(
  kw_table(hill_df_above_fungi, "Growth_Stage_Description", "Fungi",    "Aboveground", "by_year"),
  kw_table(hill_df_below_fungi, "Growth_Stage_Description", "Fungi",    "Belowground", "by_year"),
  kw_table(hill_df_above_bac,   "Growth_Stage_Description", "Bacteria", "Aboveground", "by_year"),
  kw_table(hill_df_below_bac,   "Growth_Stage_Description", "Bacteria", "Belowground", "by_year"))

all_kw <- bind_rows(year_tests, mgmt_tests, gs_tests) %>%
  mutate(
    Hill_number = as.character(q),
    p_signif = case_when(p < 0.0001 ~ "****", p < 0.001 ~ "***",
                         p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns"),
    p_formatted = ifelse(p < 0.0001, "<0.0001", as.character(signif(p, 3)))) %>%
  select(Kingdom, Compartment, Factor, Year, Hill_number, chi_sq, df, n, p, p_formatted, p_signif) %>%
  arrange(Kingdom, Compartment, Factor, Year, Hill_number)

cat("\n=== YEAR EFFECTS (overall K-W across years) ===\n");  print(all_kw %>% filter(Factor == "Year"))
cat("\n=== MANAGEMENT EFFECTS (within year) ===\n");          print(all_kw %>% filter(Factor == "Management"), n = 100)
cat("\n=== GROWTH STAGE EFFECTS (within year) ===\n");        print(all_kw %>% filter(Factor == "Growth_Stage_Description"), n = 100)

fwrite(all_kw, "alpha_kruskal_wallis_results.csv")
cat("\nExported: alpha_kruskal_wallis_results.csv\n")
cat("\n=== Alpha Diversity Complete ===\n")

# =============================================================================
# =============================================================================
#  SECTION 5
#  Figure 1 and indicator species analysis
#
#  Figure 1, Figure 3, Supplemental Table 6.
#  Also builds seq_bac and seq_fungi, which the PERMANOVA below needs, which
#  is why this section comes first.
#  Writes isa_ALL_long.csv, read by section 15.

#   Figure 1 : PCoA (sequence-keyed cross-year) + Euler, panels A-H
#   Figure 3 : ISA management heatmaps, fungi + bacteria, panels A-D
#
# ALL text bold black throughout.
#
# ISA RUNS ONCE PER STRATUM, UNCAPPED, then gets subset for plotting:
#   isa_cache   = every significant indicator (top_n = Inf)
#   heatmaps    = first ISA_TOP_N rows of each stratum (abundance-ranked)
#   isa_ALL_long.csv  = uncapped, for the 4-way overlap supplemental
#   isa_full_long.csv = capped, matches what the heatmaps plot


relabel_by_seq <- function(ps) {
  otu  <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) otu <- t(otu)
  seqs <- as.character(refseq(ps))[rownames(otu)]
  otu_collapsed <- rowsum(otu, group = seqs)
  phyloseq(otu_table(otu_collapsed, taxa_are_rows = TRUE), sample_data(ps))
}

seq_bac   <- merge_phyloseq(relabel_by_seq(physeq_rare_bac_2018),
                            relabel_by_seq(physeq_rare_bac_2019),
                            relabel_by_seq(physeq_rare_bac_2020))
seq_fungi <- merge_phyloseq(relabel_by_seq(physeq_rare_2018),
                            relabel_by_seq(physeq_rare_2019),
                            relabel_by_seq(physeq_rare_2020))

cat("seq_bac:", ntaxa(seq_bac), "seqs |", nsamples(seq_bac), "samples\n")
cat("seq_fungi:", ntaxa(seq_fungi), "seqs |", nsamples(seq_fungi), "samples\n")


stopifnot(exists("seq_bac"), exists("seq_fungi"))



# SECTION 1 — (sequence-keyed PCoA + Euler, panels A-H)

run_pcoa_plot <- function(ps, compartment, color_palette) {
  keep <- as.character(sample_data(ps)$Compartment) == compartment
  ps   <- prune_samples(sample_names(ps)[keep], ps)
  ps   <- prune_taxa(taxa_sums(ps) > 0, ps)
  
  otu  <- as(otu_table(ps), "matrix"); if (taxa_are_rows(ps)) otu <- t(otu)
  bc   <- vegdist(otu, method = "bray")
  pcoa <- cmdscale(bc, eig = TRUE, k = 2)
  pos  <- pcoa$eig[pcoa$eig > 0]
  
  plot_df <- data.frame(
    PCo1  = pcoa$points[, 1],
    PCo2  = pcoa$points[, 2],
    Group = factor(as.character(sample_data(ps)$Year), levels = names(color_palette))
  )
  ggplot(plot_df, aes(x = PCo1, y = PCo2, color = Group)) +
    geom_point(alpha = 0.85, size = 1.5) +
    stat_ellipse(level = 0.95, linetype = "dashed", linewidth = 0.7) +
    scale_color_manual(values = color_palette) +
    labs(x = sprintf("PCo1 (%.1f%%)", 100*pos[1]/sum(pos)),
         y = sprintf("PCo2 (%.1f%%)", 100*pos[2]/sum(pos))) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", panel.grid = element_blank(),
          axis.title = element_text(face = "bold", colour = "black", size = 13),
          axis.text  = element_text(face = "bold", colour = "black", size = 11),
          axis.ticks = element_line(colour = "black", linewidth = 0.6))
}

make_euler <- function(only_a, only_b, only_c, ab, ac, bc, abc, yr_colors) {
  fit <- euler(c(
    "Soybean\n(2018)" = only_a, "Wheat\n(2019)" = only_b, "Maize\n(2020)" = only_c,
    "Soybean\n(2018)&Wheat\n(2019)" = ab, "Soybean\n(2018)&Maize\n(2020)" = ac,
    "Wheat\n(2019)&Maize\n(2020)" = bc,
    "Soybean\n(2018)&Wheat\n(2019)&Maize\n(2020)" = abc))
  fill_cols <- unname(c(yr_colors["2018"], yr_colors["2019"], yr_colors["2020"]))
  as.ggplot(plot(fit,
                 fills      = list(fill = fill_cols, alpha = 0.8),
                 labels     = list(fontsize = 13, fontface = "bold", col = "black"),
                 quantities = list(fontsize = 11, fontface = "bold", col = "black")))
}

euler_counts_seq <- function(seq_obj, compartment) {
  keep <- as.character(sample_data(seq_obj)$Compartment) == compartment
  ps   <- prune_samples(sample_names(seq_obj)[keep], seq_obj)
  ps   <- prune_taxa(taxa_sums(ps) > 0, ps)
  
  otu  <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) otu <- t(otu)
  yr   <- as.character(sample_data(ps)$Year)
  
  pres <- function(y) {
    cols <- colnames(otu)[yr == y]
    if (length(cols) == 0) return(character(0))
    rownames(otu)[rowSums(otu[, cols, drop = FALSE]) > 0]
  }
  A <- pres("2018"); B <- pres("2019"); C <- pres("2020")
  
  list(
    only_a = length(setdiff(A, union(B, C))),
    only_b = length(setdiff(B, union(A, C))),
    only_c = length(setdiff(C, union(A, B))),
    ab  = length(setdiff(intersect(A, B), C)),
    ac  = length(setdiff(intersect(A, C), B)),
    bc  = length(setdiff(intersect(B, C), A)),
    abc = length(Reduce(intersect, list(A, B, C)))
  )
}

mk_euler_from_seq <- function(seq_obj, compartment, yr_colors) {
  e <- euler_counts_seq(seq_obj, compartment)
  cat(sprintf("  Euler %s: a=%d b=%d c=%d ab=%d ac=%d bc=%d abc=%d\n",
              compartment, e$only_a, e$only_b, e$only_c, e$ab, e$ac, e$bc, e$abc))
  make_euler(e$only_a, e$only_b, e$only_c, e$ab, e$ac, e$bc, e$abc, yr_colors)
}

cat("\n=== PCoA panels (sequence-keyed) ===\n")
pcoa_fa <- run_pcoa_plot(seq_fungi, "aboveground", year_colors)
pcoa_ba <- run_pcoa_plot(seq_bac,   "aboveground", year_colors)
pcoa_fb <- run_pcoa_plot(seq_fungi, "belowground", year_colors)
pcoa_bb <- run_pcoa_plot(seq_bac,   "belowground", year_colors)

cat("\n=== Euler counts (from sequences) ===\n")
cat(" Fungi:\n")
venn_fa <- mk_euler_from_seq(seq_fungi, "aboveground", year_colors)
venn_fb <- mk_euler_from_seq(seq_fungi, "belowground", year_colors)
cat(" Bacteria:\n")
venn_ba <- mk_euler_from_seq(seq_bac,   "aboveground", year_colors)
venn_bb <- mk_euler_from_seq(seq_bac,   "belowground", year_colors)

legend_plot <- pcoa_fa +
  scale_color_manual(values = year_colors,
                     labels = c("Soybean (2018)", "Wheat (2019)", "Maize (2020)")) +
  guides(color = guide_legend(title = NULL, override.aes = list(size = 8),
                              keywidth = unit(2.0, "cm"), keyheight = unit(1.1, "cm"),
                              direction = "horizontal")) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 26, face = "bold", colour = "black"),
        legend.spacing.x = unit(1.2, "cm"), legend.key = element_blank())

shared_legend   <- get_legend(legend_plot)
centered_legend <- arrangeGrob(nullGrob(), shared_legend, nullGrob(), ncol = 3, widths = c(1, 2, 1))

fungi_label  <- textGrob("Fungi",       gp = gpar(fontface = "bold", fontsize = 18, col = "black"))
bac_label    <- textGrob("Bacteria",    gp = gpar(fontface = "bold", fontsize = 18, col = "black"))
pcoa_f_label <- textGrob("PCoA",        gp = gpar(fontface = "bold.italic", fontsize = 14, col = "black"))
venn_f_label <- textGrob("ASV Overlap", gp = gpar(fontface = "bold.italic", fontsize = 14, col = "black"))
pcoa_b_label <- textGrob("PCoA",        gp = gpar(fontface = "bold.italic", fontsize = 14, col = "black"))
venn_b_label <- textGrob("ASV Overlap", gp = gpar(fontface = "bold.italic", fontsize = 14, col = "black"))
above_label  <- textGrob("Aboveground", rot = 90, gp = gpar(fontface = "bold", fontsize = 16, col = "black"))
below_label  <- textGrob("Belowground", rot = 90, gp = gpar(fontface = "bold", fontsize = 16, col = "black"))

label_it <- function(p, lab) {
  arrangeGrob(p, top = textGrob(lab, x = 0.05, just = "left",
                                gp = gpar(fontface = "bold", fontsize = 20, col = "black")))
}

top_block <- arrangeGrob(
  centered_legend,
  nullGrob(), fungi_label, nullGrob(), bac_label, nullGrob(),
  nullGrob(), pcoa_f_label, venn_f_label, pcoa_b_label, venn_b_label,
  above_label,
  label_it(ggplotGrob(pcoa_fa), "A"), label_it(venn_fa, "B"),
  label_it(ggplotGrob(pcoa_ba), "C"), label_it(venn_ba, "D"),
  below_label,
  label_it(ggplotGrob(pcoa_fb), "E"), label_it(venn_fb, "F"),
  label_it(ggplotGrob(pcoa_bb), "G"), label_it(venn_bb, "H"),
  layout_matrix = rbind(c(1,1,1,1,1), c(2,3,4,5,6), c(7,8,9,10,11),
                        c(12,13,14,15,16), c(17,18,19,20,21)),
  widths  = c(0.08, 1.2, 1, 1.2, 1),
  heights = c(0.22, 0.07, 0.05, 1, 1)
)


# SECTION 2 — ISA machinery
# run_isa() returns EVERY significant indicator, abundance-ranked.
# Capping to top N happens later via isa_head(), so the expensive
# multipatt call never depends on the plotting cutoff.


uninformative_labels <- c(
  "uncultured_bacterium", "uncultured_fungus", "metagenome",
  "uncultured", "uncultured_organism", "ambiguous_taxa",
  "unidentified", "unknown", "")

pick_tax_name <- function(tax_row, kingdom) {
  tax_cols <- c("Genus", "Family", "Order", "Class", "Phylum")
  x <- na.omit(unlist(tax_row[tax_cols]))
  x <- x[!x %in% uninformative_labels]
  x <- x[x != kingdom]
  x <- x[!grepl("^uncultured|^metagenome|^unidentified|^unknown", x, ignore.case = TRUE)]
  if (length(x) == 0) "Unclassified" else x[1]
}

run_isa <- function(ps_obj, kingdom, seed = 123) {
  otu <- as(otu_table(ps_obj), "matrix")
  if (!taxa_are_rows(ps_obj)) otu <- t(otu)
  storage.mode(otu) <- "numeric"; otu[is.na(otu)] <- 0
  
  meta        <- data.frame(sample_data(ps_obj), stringsAsFactors = FALSE)
  mgmt_levels <- c("Conventional", "No-Till", "Organic")
  
  keep <- sapply(rownames(otu), function(taxon) {
    any(sapply(mgmt_levels, function(g) {
      samps <- rownames(meta)[meta$Management == g]
      samps <- samps[samps %in% colnames(otu)]
      sum(otu[taxon, samps] > 0) >= 2
    }))
  })
  otu <- otu[keep, , drop = FALSE]
  otu <- otu[rowSums(otu) > 0, , drop = FALSE]
  if (nrow(otu) == 0) return(NULL)
  cat("    Filtered to", nrow(otu), "taxa\n")
  
  mgmt_factor <- factor(meta[colnames(otu), "Management"], levels = mgmt_levels)
  set.seed(seed)
  isa_res <- tryCatch(multipatt(t(otu), mgmt_factor, control = how(nperm = 999)),
                      error = function(e) NULL)
  if (is.null(isa_res)) return(NULL)
  
  sig <- isa_res$sign
  sig <- sig[!is.na(sig$p.value), , drop = FALSE]
  sig$p.value.adj <- p.adjust(sig$p.value, "fdr")
  sig_fdr <- sig[sig$p.value.adj <= 0.05, , drop = FALSE]
  if (nrow(sig_fdr) > 0) {
    sig <- sig_fdr; threshold <- "FDR < 0.05"
  } else {
    cat("    No FDR-significant indicators — falling back to unadjusted p < 0.05\n")
    sig <- sig[sig$p.value < 0.05, , drop = FALSE]; threshold <- "unadjusted p < 0.05"
  }
  if (nrow(sig) == 0) { cat("    No significant indicators\n"); return(NULL) }
  cat("   ", nrow(sig), "significant indicators (", threshold, ")\n")
  
  # abundance rank sets the row order; no truncation here
  sig$abund <- rowSums(otu[rownames(sig), , drop = FALSE])
  sig <- sig[order(-sig$abund), ]
  
  mgmt_cols <- c("s.Conventional", "s.No-Till", "s.Organic")
  mgmt_cols <- mgmt_cols[mgmt_cols %in% colnames(sig)]
  sig$indicator_group <- apply(sig[, mgmt_cols, drop = FALSE], 1, function(row) {
    hit <- mgmt_cols[which(row == 1)]
    if (length(hit) == 1) gsub("s\\.", "", hit) else "Multiple"
  })
  
  col_sums <- colSums(otu); col_sums[col_sums == 0] <- 1
  otu_rel  <- sweep(otu, 2, col_sums, "/")
  mat <- do.call(cbind, lapply(mgmt_levels, function(g) {
    samps <- colnames(otu)[meta[colnames(otu), "Management"] == g]
    samps <- samps[!is.na(samps)]
    if (length(samps) == 0) return(rep(NA, nrow(sig)))
    rowMeans(otu_rel[rownames(sig), samps, drop = FALSE], na.rm = TRUE)
  }))
  colnames(mat) <- mgmt_levels; rownames(mat) <- rownames(sig)
  
  tax <- as.data.frame(tax_table(ps_obj))
  base_names <- sapply(rownames(sig), function(asv) {
    if (asv %in% rownames(tax)) pick_tax_name(tax[asv, ], kingdom) else asv
  })
  
  list(matrix          = mat,
       labels_plain    = unname(base_names),
       labels_with_id  = unname(paste0(base_names, " (", rownames(sig), ")")),
       asv_id          = rownames(sig),
       indicator_group = sig$indicator_group,
       abund           = sig$abund,
       threshold       = threshold,
       n_total         = nrow(sig)) 
}

# take the top n rows of an ISA result; rows are already abundance-ranked
isa_head <- function(res, n) {
  if (is.null(res)) return(NULL)
  k <- min(n, nrow(res$matrix))
  res$matrix          <- res$matrix[seq_len(k), , drop = FALSE]
  res$labels_plain    <- res$labels_plain[seq_len(k)]
  res$labels_with_id  <- res$labels_with_id[seq_len(k)]
  res$asv_id          <- res$asv_id[seq_len(k)]
  res$indicator_group <- res$indicator_group[seq_len(k)]
  res$abund           <- res$abund[seq_len(k)]
  res
}



# ISA slice lists

fungi_a <- list("2018" = psr_2018_a,     "2019" = psr_2019_a,     "2020" = psr_2020_a)
fungi_b <- list("2018" = psr_2018_b,     "2019" = psr_2019_b,     "2020" = psr_2020_b)
bac_a   <- list("2018" = psr_2018_a_bac, "2019" = psr_2019_a_bac, "2020" = psr_2020_a_bac)
bac_b   <- list("2018" = psr_2018_b_bac, "2019" = psr_2019_b_bac, "2020" = psr_2020_b_bac)



# BUILD THE CACHE — multipatt runs here and nowhere else, uncapped
# 12 strata: 2 kingdoms x 2 compartments x 3 years
# Set FORCE_ISA <- TRUE to recompute after changing the data.

ISA_TOP_N <- 15          # heatmap cutoff only; does not affect the cache
if (!exists("FORCE_ISA")) FORCE_ISA <- FALSE

isa_strata <- list(
  list(km = "Fungi",    tis = "aboveground", lst = fungi_a),
  list(km = "Fungi",    tis = "belowground", lst = fungi_b),
  list(km = "Bacteria", tis = "aboveground", lst = bac_a),
  list(km = "Bacteria", tis = "belowground", lst = bac_b)
)
isa_seeds <- c("2018" = 123, "2019" = 456, "2020" = 789)

if (!exists("isa_cache") || FORCE_ISA) {
  cat("\n########## RUNNING ISA (12 strata, uncapped, cached) ##########\n")
  isa_cache <- list()
  for (s in isa_strata) for (yr in c("2018", "2019", "2020")) {
    key <- paste(s$km, s$tis, yr, sep = "|")
    cat("ISA:", key, "\n")
    isa_cache[[key]] <- run_isa(s$lst[[yr]], kingdom = s$km, seed = isa_seeds[yr])
  }
} else {
  cat("\n########## ISA cache found — skipping recomputation ##########\n")
  cat("Set FORCE_ISA <- TRUE and rerun to refresh.\n")
}

get_isa <- function(kingdom, tissue, year) isa_cache[[paste(kingdom, tissue, year, sep = "|")]]

make_year_heatmap <- function(res, year, show_asv) {
  mat      <- res$matrix
  row_labs <- if (show_asv) res$labels_with_id else res$labels_plain
  col_fun  <- colorRamp2(c(0, max(mat, na.rm = TRUE)), c("white", year_colors[year]))
  crop      <- switch(year, "2018" = "Soybean", "2019" = "Wheat", "2020" = "Maize")
  col_title <- paste0(crop, "\n(", year, ")")
  if (res$threshold != "FDR < 0.05") col_title <- paste0(col_title, "\n(unadjusted p<0.05)")
  Heatmap(
    mat, name = paste0("Rel. Abund. ", year), col = col_fun,
    cluster_rows = FALSE, cluster_columns = FALSE,
    column_title = col_title,
    column_title_gp = gpar(fontsize = 14, fontface = "bold", col = "black"),
    row_labels = row_labs, row_names_side = "left",
    row_names_gp = gpar(fontsize = 10, fontface = "bold", col = "black"),
    column_names_gp  = gpar(fontsize = 11, fontface = "bold", col = "black"),
    column_names_rot = 45,
    show_heatmap_legend = TRUE,
    heatmap_legend_param = list(
      title = "Rel. Abund.", title_gp = gpar(fontsize = 11, fontface = "bold", col = "black"),
      labels_gp = gpar(fontsize = 9, fontface = "bold", col = "black"),
      legend_direction = "horizontal", legend_width = unit(2.5, "cm")),
    na_col = "grey95", border = TRUE, border_gp = gpar(col = "black", lwd = 1.5),
    height = unit(min(nrow(mat) * 5.5, 130), "mm")
  )
}

heatmap_grob <- function(res, year, show_asv, w = 5.5, h = 10) {
  if (is.null(res)) return(nullGrob())
  ht <- make_year_heatmap(res, year, show_asv)
  grid.grabExpr(draw(ht, heatmap_legend_side = "bottom", padding = unit(c(2,2,8,2), "mm")),
                width = unit(w, "cm"), height = unit(h, "cm"))
}

add_panel_letter <- function(g, letter) {
  arrangeGrob(g, top = textGrob(letter, x = 0.02, hjust = 0,
                                gp = gpar(fontface = "bold", fontsize = 20, col = "black")))
}

build_isa_block <- function(kingdom, letter_above, letter_below, show_asv = FALSE) {
  years <- c("2018", "2019", "2020")
  
  mk <- function(tissue, comp_letter) {
    lapply(seq_along(years), function(i) {
      res <- isa_head(get_isa(kingdom, tissue, years[i]), ISA_TOP_N)
      g   <- heatmap_grob(res, years[i], show_asv)
      if (i == 1) g <- add_panel_letter(g, comp_letter)
      g
    })
  }
  ga <- mk("aboveground", letter_above)
  gb <- mk("belowground", letter_below)
  
  kingdom_lab <- textGrob(kingdom,       gp = gpar(fontface = "bold", fontsize = 20, col = "black"))
  above_lab   <- textGrob("Aboveground", gp = gpar(fontface = "bold", fontsize = 16, col = "black"))
  below_lab   <- textGrob("Belowground", gp = gpar(fontface = "bold", fontsize = 16, col = "black"))
  spacer      <- nullGrob()
  
  arrangeGrob(
    kingdom_lab, above_lab, below_lab, spacer,
    ga[[1]], ga[[2]], ga[[3]], gb[[1]], gb[[2]], gb[[3]],
    layout_matrix = rbind(c(1,1,1, 1,1,1,1), c(2,2,2, 4,3,3,3), c(5,6,7, 4,8,9,10)),
    widths  = c(1, 1, 1, 0.05, 1, 1, 1),
    heights = c(0.06, 0.06, 1)
  )
}

# SECTION 2a — STANDALONE ISA STRIPS, WITH Zotu ids (verification)

cat("\n########## STANDALONE ISA (WITH Zotu ids, verification) ##########\n")

fungi_block_ids    <- build_isa_block("Fungi",    "A", "B", show_asv = TRUE)
grid.newpage(); grid.draw(fungi_block_ids)

bacteria_block_ids <- build_isa_block("Bacteria", "C", "D", show_asv = TRUE)
grid.newpage(); grid.draw(bacteria_block_ids)



# SECTION 3 — FIGURE VERSION, ISA WITHOUT Zotu ids, panels A-D

cat("\n########## FIGURE VERSION (ISA WITHOUT Zotu ids) ##########\n")

fungi_block    <- build_isa_block("Fungi",    "A", "B", show_asv = FALSE)
bacteria_block <- build_isa_block("Bacteria", "C", "D", show_asv = FALSE)

isa_figure <- arrangeGrob(fungi_block, bacteria_block, ncol = 1, heights = c(1, 1))

# SECTION 2b — EXPORT ISA TABLES
#   isa_full_long.csv     — top ISA_TOP_N per stratum, matches the heatmaps
#   isa_group_summary.csv — indicator-group counts + % (from the FULL set,
#                           which is where the manuscript percentages come from)
#   isa_ALL_long.csv      — every significant indicator, for the 4-way
#                           overlap supplemental figure

isa_result_to_df <- function(res, kingdom, year, tissue) {
  if (is.null(res)) {
    return(data.frame(Kingdom=kingdom, Year=year, Tissue=tissue,
                      ASV=NA, ASV_label=NA, Indicator_group=NA,
                      Abundance=NA, Threshold="none (no indicators)",
                      Conventional=NA, `No-Till`=NA, Organic=NA,
                      check.names=FALSE, stringsAsFactors=FALSE))
  }
  m <- res$matrix
  data.frame(
    Kingdom = kingdom, Year = year, Tissue = tissue,
    ASV             = res$asv_id,          
    ASV_label       = res$labels_with_id,
    Indicator_group = res$indicator_group,
    Abundance       = res$abund,
    Threshold       = res$threshold,
    Conventional = m[, "Conventional"],
    `No-Till`    = m[, "No-Till"],
    Organic      = m[, "Organic"],
    check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL
  )
}

collect_isa <- function(cap = Inf) {
  out <- list()
  for (s in isa_strata) for (yr in c("2018","2019","2020")) {
    res <- get_isa(s$km, s$tis, yr)
    if (is.finite(cap)) res <- isa_head(res, cap)
    out[[length(out)+1]] <- isa_result_to_df(res, s$km, yr, s$tis)
  }
  do.call(rbind, out)
}

isa_all    <- collect_isa(cap = Inf)        
isa_capped <- collect_isa(cap = ISA_TOP_N)  

# group summary computed on the FULL set
valid <- isa_all[!is.na(isa_all$Indicator_group), , drop = FALSE]
summ  <- aggregate(ASV ~ Kingdom + Year + Tissue + Threshold + Indicator_group,
                   data = valid, FUN = length)
names(summ)[names(summ) == "ASV"] <- "n"
tot <- aggregate(n ~ Kingdom + Year + Tissue, data = summ, FUN = sum)
names(tot)[names(tot) == "n"] <- "stratum_total"
summ <- merge(summ, tot, by = c("Kingdom","Year","Tissue"))
summ$pct <- round(100 * summ$n / summ$stratum_total, 1)
summ <- summ[order(summ$Kingdom, summ$Year, summ$Tissue, -summ$n), ]

write.csv(isa_all,    "isa_ALL_long.csv",      row.names = FALSE)
write.csv(isa_capped, "isa_full_long.csv",     row.names = FALSE)
write.csv(summ,       "isa_group_summary.csv", row.names = FALSE)

cat("\nWrote isa_ALL_long.csv (", nrow(isa_all), " rows, uncapped)\n", sep = "")
cat("Wrote isa_full_long.csv (", nrow(isa_capped), " rows, top ", ISA_TOP_N, ")\n", sep = "")
cat("Wrote isa_group_summary.csv (", nrow(summ), " rows)\n", sep = "")
print(summ, row.names = FALSE)


# OUTPUT — two separate figures

pdf("Fig1_PCoA_Euler.pdf", width = 20, height = 12)
grid.newpage(); grid.draw(top_block)
dev.off()

pdf("Fig3_ISA.pdf", width = 20, height = 16)
grid.newpage(); grid.draw(isa_figure)
dev.off()

# on-screen check
grid.newpage(); grid.draw(top_block)
grid.newpage(); grid.draw(isa_figure)

# =============================================================================
# =============================================================================
#  SECTION 6
#  PERMANOVA
#
#  Supplemental Tables 4 and 5.

stopifnot(exists("seq_bac"), exists("seq_fungi"),
          exists("physeq_rare_2018"), exists("physeq_rare_bac_2018"))

run_one <- function(ps, formula_str) {
  otu <- as(otu_table(ps), "matrix"); if (taxa_are_rows(ps)) otu <- t(otu)
  meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
  for (col in c("Year","Management","Growth_Stage_Description","Compartment"))
    if (col %in% colnames(meta)) meta[[col]] <- factor(meta[[col]])
  meta <- meta[rownames(otu), ]
  bc <- vegdist(otu, method = "bray")
  set.seed(42)
  res <- adonis2(as.formula(paste("bc ~", formula_str)), data = meta, permutations = 999, by = "terms")
  pr <- as.data.frame(res)[formula_str, ]
  list(R2 = pr$R2, F_stat = pr$F, p = pr$`Pr(>F)`, dist = bc, meta = meta, n = nrow(otu))
}
run_bd <- function(bc_dist, meta, factor_col) {
  groups <- meta[[factor_col]]
  if (length(unique(groups)) < 2) return(list(F_stat = NA, p = NA))
  bd <- betadisper(bc_dist, groups); bdt <- permutest(bd, permutations = 999)
  list(F_stat = bdt$tab$F[1], p = bdt$tab$`Pr(>F)`[1])
}
# split a phyloseq by an explicit logical (no subset_samples NSE)
slice_ps <- function(ps, keep_logical) {
  ps2 <- prune_samples(sample_names(ps)[keep_logical], ps)
  prune_taxa(taxa_sums(ps2) > 0, ps2)
}

rows <- list()
add <- function(yr, king, comp, test, strat_by, strat_level, n, pf, pr2, pp, df, dp) {
  rows[[length(rows) + 1]] <<- data.frame(
    Year = as.character(yr), Kingdom = king, Compartment = comp,
    Test = test, Stratify_by = strat_by, Stratify_level = strat_level, n_samples = n,
    PERMANOVA_F = round(pf, 2), PERMANOVA_R2 = round(pr2, 3), PERMANOVA_P_value = round(pp, 3),
    DISPERSION_F = round(df, 2), DISPERSION_P_value = round(dp, 3), stringsAsFactors = FALSE)
}

crop_map      <- c("2018" = "Soybean", "2019" = "Wheat", "2020" = "Maize")
years         <- c("2018","2019","2020")
compartments  <- c("aboveground","belowground")
managements   <- c("Conventional","No-Till","Organic")
growth_stages <- c("Vegetative","Inflorescence","Reproductive")

# SECTION 1: CROP SPECIES (Year effect) — SEQUENCE-KEYED

cat("\n=== Section 1: Crop Species (sequence-keyed) ===\n")
seq_list <- list(list(ps = seq_fungi, name = "Fungi"), list(ps = seq_bac, name = "Bacteria"))
for (k in seq_list) {
  ps_tmp <- prune_taxa(taxa_sums(k$ps) > 0, k$ps)
  cat("Crop Species:", k$name, "| Pooled\n")
  res <- run_one(ps_tmp, "Year"); bd <- run_bd(res$dist, res$meta, "Year")
  add("All", k$name, "Pooled", "Crop_Species", NA, "Pooled",
      res$n, res$F_stat, res$R2, res$p, bd$F_stat, bd$p)
  for (comp in compartments) {
    keep <- as.character(sample_data(k$ps)$Compartment) == comp
    ps_c <- slice_ps(k$ps, keep)
    cat("Crop Species:", k$name, "|", comp, "\n")
    res <- run_one(ps_c, "Year"); bd <- run_bd(res$dist, res$meta, "Year")
    add("All", k$name, comp, "Crop_Species", "Compartment", comp,
        res$n, res$F_stat, res$R2, res$p, bd$F_stat, bd$p)
  }
}

# SECTION 2: WITHIN YEAR × KINGDOM × COMPARTMENT
# each kingdom loops its OWN per-year objects

cat("\n=== Section 2: Within-Year ===\n")
kingdom_year_obj <- list(
  Fungi    = function(y) get(paste0("physeq_rare_", y)),
  Bacteria = function(y) get(paste0("physeq_rare_bac_", y)))

for (yr in years) {
  for (king in names(kingdom_year_obj)) {
    ps_year <- kingdom_year_obj[[king]](yr)
    for (comp in compartments) {
      keep_comp <- as.character(sample_data(ps_year)$Compartment) == comp
      ps_sub <- slice_ps(ps_year, keep_comp)
      n_sub <- nsamples(ps_sub)
      cat(crop_map[yr], yr, king, comp, "n =", n_sub, "\n")
      if (n_sub < 10) { cat("  Skipping\n"); next }
      
      otu_tmp <- as(otu_table(ps_sub), "matrix"); if (taxa_are_rows(ps_sub)) otu_tmp <- t(otu_tmp)
      meta_tmp <- data.frame(sample_data(ps_sub), stringsAsFactors = FALSE)
      meta_tmp$Management <- factor(meta_tmp$Management)
      meta_tmp$Growth_Stage_Description <- factor(meta_tmp$Growth_Stage_Description)
      meta_tmp <- meta_tmp[rownames(otu_tmp), ]
      bc_sub <- vegdist(otu_tmp, method = "bray")
      yr_label <- yr
      
      # Pooled Management
      res <- run_one(ps_sub, "Management"); bd <- run_bd(bc_sub, meta_tmp, "Management")
      add(yr_label, king, comp, "Management", NA, "Pooled",
          res$n, res$F_stat, res$R2, res$p, bd$F_stat, bd$p); yr_label <- NA
      
      # Pooled Growth Stage
      res <- run_one(ps_sub, "Growth_Stage_Description"); bd <- run_bd(bc_sub, meta_tmp, "Growth_Stage_Description")
      add(NA, NA, NA, "Growth_Stage", NA, "Pooled",
          res$n, res$F_stat, res$R2, res$p, bd$F_stat, bd$p)
      
      # Management stratified by Growth Stage
      for (gs in growth_stages) {
        keep_gs <- as.character(sample_data(ps_sub)$Growth_Stage_Description) == gs
        ps_gs <- slice_ps(ps_sub, keep_gs)
        if (nsamples(ps_gs) < 6) next
        if (length(unique(sample_data(ps_gs)$Management)) < 2) next
        res <- run_one(ps_gs, "Management")
        otu_gs <- as(otu_table(ps_gs), "matrix"); if (taxa_are_rows(ps_gs)) otu_gs <- t(otu_gs)
        meta_gs <- data.frame(sample_data(ps_gs)); meta_gs$Management <- factor(meta_gs$Management)
        meta_gs <- meta_gs[rownames(otu_gs), ]
        bd <- run_bd(vegdist(otu_gs, method = "bray"), meta_gs, "Management")
        if (gs == growth_stages[1]) add(NA, NA, NA, "Management", "Growth_Stage", gs,
                                        res$n, res$F_stat, res$R2, res$p, bd$F_stat, bd$p)
        else add(NA, NA, NA, NA, NA, gs, res$n, res$F_stat, res$R2, res$p, bd$F_stat, bd$p)
      }
      
      # Growth Stage stratified by Management
      for (mg in managements) {
        keep_mg <- as.character(sample_data(ps_sub)$Management) == mg
        ps_mg <- slice_ps(ps_sub, keep_mg)
        if (nsamples(ps_mg) < 6) next
        if (length(unique(sample_data(ps_mg)$Growth_Stage_Description)) < 2) next
        res <- run_one(ps_mg, "Growth_Stage_Description")
        otu_mg <- as(otu_table(ps_mg), "matrix"); if (taxa_are_rows(ps_mg)) otu_mg <- t(otu_mg)
        meta_mg <- data.frame(sample_data(ps_mg)); meta_mg$Growth_Stage_Description <- factor(meta_mg$Growth_Stage_Description)
        meta_mg <- meta_mg[rownames(otu_mg), ]
        bd <- run_bd(vegdist(otu_mg, method = "bray"), meta_mg, "Growth_Stage_Description")
        if (mg == managements[1]) add(NA, NA, NA, "Growth_Stage", "Management", mg,
                                      res$n, res$F_stat, res$R2, res$p, bd$F_stat, bd$p)
        else add(NA, NA, NA, NA, NA, mg, res$n, res$F_stat, res$R2, res$p, bd$F_stat, bd$p)
      }
    }
  }
}

final <- bind_rows(rows)
cat("\n=== Complete Table ===\n"); print(as.data.frame(final), row.names = FALSE)
write.csv(final, "PERMANOVA_COMPLETE_TABLE_NEW.csv", row.names = FALSE)
cat("\nSaved. Total rows:", nrow(final), "\n")

# =============================================================================
# =============================================================================

#  SECTION 7
#  Pairwise PERMANOVA heatmap
#
#  Figure 2.

# PAIRWISE PERMANOVA (+ betadisper flag) -> "which level" heatmap
# Two facets: Management + Growth stage

set.seed(517)
N_PERM  <- 999
MIN_GRP <- 3
MIN_TOT <- 6

pairwise_permanova <- function(ps, group_var) {
  if (is.null(ps) || nsamples(ps) < MIN_TOT) return(NULL)
  g <- as.character(sample_data(ps)[[group_var]])
  groups <- sort(unique(g[!is.na(g)]))
  if (length(groups) < 2) return(NULL)
  pairs <- combn(groups, 2, simplify = FALSE)
  purrr::map_dfr(pairs, function(pr) {
    sub  <- prune_samples(g %in% pr, ps); sub <- prune_taxa(taxa_sums(sub) > 0, sub)
    meta <- data.frame(sample_data(sub)); meta[[group_var]] <- droplevels(factor(meta[[group_var]]))
    tb <- table(meta[[group_var]])
    if (length(tb) < 2 || min(tb) < MIN_GRP || nsamples(sub) < MIN_TOT)
      return(tibble(contrast = paste(pr, collapse = " vs "), R2 = NA_real_, p = NA_real_, n = nsamples(sub)))
    d   <- phyloseq::distance(sub, method = "bray")
    fit <- adonis2(as.formula(paste("d ~", group_var)), data = meta, permutations = N_PERM)
    tibble(contrast = paste(pr, collapse = " vs "), R2 = fit$R2[1], p = fit$`Pr(>F)`[1], n = nsamples(sub))
  })
}

pairwise_betadisper <- function(ps, group_var) {
  if (is.null(ps) || nsamples(ps) < MIN_TOT) return(NULL)
  g <- as.character(sample_data(ps)[[group_var]])
  groups <- sort(unique(g[!is.na(g)]))
  if (length(groups) < 2) return(NULL)
  pairs <- combn(groups, 2, simplify = FALSE)
  purrr::map_dfr(pairs, function(pr) {
    sub  <- prune_samples(g %in% pr, ps); sub <- prune_taxa(taxa_sums(sub) > 0, sub)
    meta <- data.frame(sample_data(sub)); grp <- droplevels(factor(meta[[group_var]]))
    tb <- table(grp)
    if (length(tb) < 2 || min(tb) < MIN_GRP || nsamples(sub) < MIN_TOT)
      return(tibble(contrast = paste(pr, collapse = " vs "), disp_p = NA_real_))
    d  <- phyloseq::distance(sub, method = "bray")
    bd <- betadisper(d, grp); pt <- permutest(bd, permutations = N_PERM)
    tibble(contrast = paste(pr, collapse = " vs "), disp_p = pt$tab$`Pr(>F)`[1])
  })
}

spec <- expand_grid(
  kingdom     = c("Fungi", "Bacteria"),
  compartment = c("aboveground", "belowground"),
  year        = c("2018", "2019", "2020")
) %>%
  mutate(comp_suf = if_else(compartment == "aboveground", "a", "b"),
         obj_name = if_else(kingdom == "Fungi",
                            paste0("psr_", year, "_", comp_suf),
                            paste0("psr_", year, "_", comp_suf, "_bac")))

run_all_strata <- function(fun, group_var, empty_row) {
  spec %>%
    mutate(res = purrr::map(obj_name, function(nm) {
      if (!exists(nm)) return(NULL); fun(get(nm), group_var = group_var)
    })) %>%
    mutate(res = purrr::map2(res, obj_name, ~ if (is.null(.x)) empty_row else .x)) %>%
    unnest(res)
}

empty_perm <- tibble(contrast = NA_character_, R2 = NA_real_, p = NA_real_, n = 0L)
empty_disp <- tibble(contrast = NA_character_, disp_p = NA_real_)

pw_mgmt  <- run_all_strata(pairwise_permanova, "Management",               empty_perm) %>% mutate(facet = "Management")
pw_stage <- run_all_strata(pairwise_permanova, "Growth_Stage_Description", empty_perm) %>% mutate(facet = "Growth stage")
pw <- bind_rows(pw_mgmt, pw_stage) %>%
  mutate(p_adj = p.adjust(p, method = "BH"), sig = !is.na(p_adj) & p_adj < 0.05)

disp_mgmt  <- run_all_strata(pairwise_betadisper, "Management",               empty_disp) %>% mutate(facet = "Management")
disp_stage <- run_all_strata(pairwise_betadisper, "Growth_Stage_Description", empty_disp) %>% mutate(facet = "Growth stage")
disp <- bind_rows(disp_mgmt, disp_stage) %>%
  mutate(disp_p_adj = p.adjust(disp_p, method = "BH"), disp_sig = !is.na(disp_p_adj) & disp_p_adj < 0.05)

pw <- pw %>%
  left_join(disp %>% select(kingdom, compartment, year, contrast, facet, disp_sig),
            by = c("kingdom", "compartment", "year", "contrast", "facet"))

mgmt_levels  <- c("Conventional vs No-Till", "Conventional vs Organic", "No-Till vs Organic")
stage_levels <- c("Inflorescence vs Vegetative", "Reproductive vs Vegetative", "Inflorescence vs Reproductive")

gmax <- pw %>% filter(sig) %>% pull(R2) %>% max(na.rm = TRUE)

pw_plot <- pw %>%
  filter(!is.na(contrast)) %>%
  mutate(
    kc        = factor(paste(kingdom, compartment), levels = kc_order),
    year_lab  = factor(crop_labels[year],
                       levels = crop_labels[c("2020", "2019", "2018")]),
    contrast  = factor(contrast, levels = c(mgmt_levels, stage_levels)),
    facet     = factor(facet, levels = c("Management", "Growth stage")),
    R2sig     = if_else(sig, R2, NA_real_),
    star = case_when(is.na(p_adj) ~ "", p_adj < 0.001 ~ "***",
                     p_adj < 0.01 ~ "**", p_adj < 0.05 ~ "*", TRUE ~ ""),
    label = case_when(is.na(R2) ~ "\u2014",
                      sig       ~ paste0(sprintf("%.3f", R2), star),
                      TRUE      ~ "n.s."),
    frac = R2sig / gmax,
    txt_col   = if_else(!is.na(frac) & frac > 0.62, "white", "black"),
    flag_disp = sig & !is.na(disp_sig) & disp_sig
  )

make_block <- function(kc_name, tag_start, show_x = FALSE) {
  d <- dplyr::filter(pw_plot, kc == kc_name)
  
  tags <- tibble(
    facet = factor(c("Management", "Growth stage"),
                   levels = c("Management", "Growth stage")),
    tag   = LETTERS[tag_start:(tag_start + 1)]
  )
  
  ggplot(d, aes(contrast, year_lab)) +
    geom_tile(aes(fill = R2sig), colour = "white", linewidth = 0.9) +
    geom_tile(data = dplyr::filter(d, flag_disp),
              fill = NA, colour = "black", linewidth = 0.5, linetype = "22") +
    geom_text(aes(label = label, colour = txt_col), size = 3.0,
              fontface = "bold", show.legend = FALSE) +
    geom_text(data = tags, aes(label = tag), x = -Inf, y = Inf,
              hjust = 1.5, vjust = -0.3, size = 5, fontface = "bold",
              colour = "black", inherit.aes = FALSE) +
    facet_grid2(. ~ facet, scales = "free_x", space = "free_x") +
    scale_fill_gradient(low = "white", high = unname(kc_colors[kc_name]),
                        limits = c(0, gmax), na.value = "#e8e6df",
                        name = expression(bold(R^2)),
                        guide = guide_colourbar(barwidth = unit(0.35, "cm"),
                                                barheight = unit(2.2, "cm"),
                                                ticks.colour = "black",
                                                frame.colour = "black")) +
    scale_colour_identity() +
    coord_cartesian(clip = "off") +
    labs(title = kc_name, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title        = element_text(face = "bold", colour = "black",
                                       size = 11, hjust = 0,
                                       margin = ggplot2::margin(b = 4)),
      panel.grid        = element_blank(),
      axis.text.y       = element_text(face = "bold", colour = "black"),
      axis.text.x       = if (show_x)
        element_text(angle = 25, hjust = 1, face = "bold", colour = "black")
      else element_blank(),
      strip.text.x      = if (tag_start == 1)
        element_text(face = "bold", colour = "black", size = 11)
      else element_blank(),
      legend.title      = element_text(face = "bold", colour = "black", size = 9),
      legend.text       = element_text(colour = "black", size = 8),
      panel.spacing.x   = unit(8, "pt"),
      plot.margin       = ggplot2::margin(4, 4, 4, 22)
    )
}

blocks <- list(
  make_block(kc_order[1], 1),                 # A, B  Fungi aboveground
  make_block(kc_order[2], 3),                 # C, D  Fungi belowground
  make_block(kc_order[3], 5),                 # E, F  Bacteria aboveground
  make_block(kc_order[4], 7, show_x = TRUE)   # G, H  Bacteria belowground
)

p_pairwise <- plot_grid(plotlist = blocks, ncol = 1,
                        align = "v", axis = "lr",
                        rel_heights = c(1, 1, 1, 1.45))

p_pairwise
# =============================================================================
# =============================================================================

#  SECTION 8
#  Distance-based redundancy analysis
#
#  Supplemental Table 9.
# db-RDA SUPPLEMENTARY TABLE S2 — per-year clean slices
# Terms = Management + Growth_Stage (matches the PERMANOVA table).

extract_dbrda_stats <- function(ps_obj, year, compartment, kingdom) {
  Y <- as(otu_table(ps_obj), "matrix"); if (taxa_are_rows(ps_obj)) Y <- t(Y)
  Y <- Y[, colSums(Y) > 0, drop = FALSE]
  if (nrow(Y) < 10) return(NULL)
  
  meta <- as(sample_data(ps_obj), "data.frame"); meta <- meta[rownames(Y), , drop = FALSE]
  # normalize stage spelling + factor
  if ("Growth_Stage_Description" %in% names(meta)) {
    gs <- gsub("Infloresence", "Inflorescence", as.character(meta$Growth_Stage_Description), ignore.case = TRUE)
    meta$Growth_Stage_Description <- factor(gs, levels = c("Vegetative","Inflorescence","Reproductive"))
  }
  if ("Management" %in% names(meta)) meta$Management <- factor(as.character(meta$Management))
  
  # keep a term only if it has >=2 levels each with >=2 samples
  usable <- function(col) {
    if (!col %in% names(meta)) return(FALSE)
    tb <- table(meta[[col]]); sum(tb >= 2) >= 2
  }
  terms <- c()
  if (usable("Management")) terms <- c(terms, "Management")
  if (usable("Growth_Stage_Description")) terms <- c(terms, "Growth_Stage_Description")
  if (length(terms) == 0) { cat("  skip", kingdom, year, compartment, "- no usable terms\n"); return(NULL) }
  
  # drop unused levels so capscale doesn't see empty ones
  for (t in terms) meta[[t]] <- droplevels(meta[[t]])
  
  D   <- vegdist(Y, method = "bray")
  fit <- tryCatch(capscale(as.formula(paste("D ~", paste(terms, collapse = " + "))), data = meta),
                  error = function(e) { cat("  capscale failed:", kingdom, year, compartment, "-", conditionMessage(e), "\n"); NULL })
  if (is.null(fit)) return(NULL)
  
  set.seed(42); ctrl <- how(nperm = 999)
  aov_overall <- anova.cca(fit, permutations = ctrl)
  aov_terms   <- anova.cca(fit, by = "term", permutations = ctrl)
  r2_adj    <- RsquareAdj(fit)$adj.r.squared
  p_overall <- aov_overall["Model", "Pr(>F)"]
  eig       <- eigenvals(fit, model = "constrained")
  cap1 <- round(100 * eig[1] / sum(eig), 1)
  cap2 <- if (length(eig) > 1) round(100 * eig[2] / sum(eig), 1) else NA
  term_rows <- rownames(aov_terms)[rownames(aov_terms) != "Residual"]
  crop <- switch(year, "2018"="Soybean", "2019"="Wheat", "2020"="Maize")
  
  data.frame(
    Kingdom = kingdom, Year = year, Crop = crop, Compartment = compartment,
    Term = c("Overall Model", term_rows),
    df   = c(aov_overall["Model","Df"], aov_terms[term_rows,"Df"]),
    F_statistic = round(c(aov_overall["Model","F"], aov_terms[term_rows,"F"]), 2),
    R2_adj = c(round(r2_adj,3), rep(NA, length(term_rows))),
    p_value = round(c(p_overall, aov_terms[term_rows,"Pr(>F)"]), 3),
    CAP1_pct = cap1, CAP2_pct = cap2, stringsAsFactors = FALSE)
}

dbrda_results <- list(); counter <- 1
for (yr in c("2018","2019","2020")) {
  for (suf in c("a","b")) {
    comp <- ifelse(suf == "a", "Aboveground", "Belowground")
    dbrda_results[[counter]] <- extract_dbrda_stats(get(paste0("psr_", yr, "_", suf)),        yr, comp, "Fungi");    counter <- counter + 1
    dbrda_results[[counter]] <- extract_dbrda_stats(get(paste0("psr_", yr, "_", suf, "_bac")), yr, comp, "Bacteria"); counter <- counter + 1
  }
}
dbrda_table <- do.call(rbind, dbrda_results)
dbrda_table$Significance <- ifelse(dbrda_table$p_value < 0.001, "***",
                                   ifelse(dbrda_table$p_value < 0.01,  "**",
                                          ifelse(dbrda_table$p_value < 0.05,  "*", "ns")))
write.csv(dbrda_table, "Supplementary_Table_S2_dbRDA.csv", row.names = FALSE)
cat("\nSaved Supplementary_Table_S2_dbRDA.csv | rows:", nrow(dbrda_table),
    "| NA cells:", sum(is.na(dbrda_table)), "\n")
print(dbrda_table, row.names = FALSE)
# =============================================================================
# =============================================================================
#  SECTION 9
#  Compartment overlap and cultivar confounding
#  Supplemental Figure 5 (aboveground against belowground overlap) and
#  Supplemental Figure 6 (conventional against no-till within Pioneer).

mgmt_cols <- c("Conventional" = "#bd461d", "No-Till" = "#0d3660")

ab_bg_sets <- function(ps) {
  md  <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
  comp <- as.character(md$Compartment)
  seqs <- if (!is.null(refseq(ps, errorIfNULL = FALSE))) as.character(refseq(ps)) else taxa_names(ps)
  names(seqs) <- taxa_names(ps)
  otu <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) otu <- t(otu)
  seqs_in <- function(wc) {
    samp <- rownames(md)[comp == wc]
    keep <- rowSums(otu[, samp, drop = FALSE]) > 0
    unique(seqs[names(keep)[keep]])
  }
  list(ag = seqs_in("aboveground"), bg = seqs_in("belowground"))
}

mk_panel <- function(ps, title, fill_ab, fill_bg) {
  s <- ab_bg_sets(ps)
  shared  <- length(intersect(s$ag, s$bg))
  ag_only <- length(setdiff(s$ag, s$bg))
  bg_only <- length(setdiff(s$bg, s$ag))
  total   <- length(union(s$ag, s$bg))
  pct     <- round(100 * shared / total, 1)
  fit <- euler(c("Aboveground" = ag_only, "Belowground" = bg_only,
                 "Aboveground&Belowground" = shared))
  p <- plot(fit,
            fills  = list(fill = c(fill_ab, fill_bg), alpha = 0.55),
            edges  = list(col = "grey40", lwd = 1),
            labels = list(labels = c("Aboveground","Belowground"), fontsize = 10, font = 2),
            quantities = list(fontsize = 9))
  # wrap the euler grob with a title block above it, with padding
  arrangeGrob(p,
              top = textGrob(paste0(title, "\n", pct, "% shared"),
                             gp = gpar(fontface = "bold", fontsize = 11),
                             just = "centre"))
}

pA <- mk_panel(psr_2018,     "A  Fungi \u2014 Soybean (2018)",   col_ab, col_bg)
pB <- mk_panel(psr_2019,     "B  Fungi \u2014 Wheat (2019)",     col_ab, col_bg)
pC <- mk_panel(psr_2020,     "C  Fungi \u2014 Maize (2020)",     col_ab, col_bg)
pD <- mk_panel(psr_2018_bac, "D  Bacteria \u2014 Soybean (2018)", col_ab, col_bg)
pE <- mk_panel(psr_2019_bac, "E  Bacteria \u2014 Wheat (2019)",   col_ab, col_bg)
pF <- mk_panel(psr_2020_bac, "F  Bacteria \u2014 Maize (2020)",   col_ab, col_bg)

grid.newpage()
grid.arrange(
  textGrob("Aboveground vs. Belowground Sequence Overlap",
           gp = gpar(fontface = "bold", fontsize = 15)),
  arrangeGrob(pA, pB, pC, ncol = 3),
  arrangeGrob(pD, pE, pF, ncol = 3),
  heights = c(0.10, 1, 1)
)

# CULTIVAR ANALYSIS
# Cultivar is fully confounded with management (Pioneer = Conv/NoTill, Viking = Organic).
# (2) demonstrates the confounding via varpart (which CANNOT partition them — the failure IS the result),
# (3) tests Conv-vs-NoTill within Pioneer (cultivar held constant),
# (4) draws the Pioneer-only Conv-vs-NoTill PCoA figure.
# 2019 excluded throughout: single cultivar (Pioneer 25R40).

library(phyloseq); library(vegan); library(ggplot2); library(grid); library(gridExtra)

# 0. Cultivar backfill for bacteria

cult_for <- function(yr, mgmt) {
  pio <- c("2018"="Pioneer P22T69R","2019"="Pioneer 25R40","2020"="Pioneer P0306Q")
  vik <- c("2018"="Viking O.2188AT12N","2020"="Viking O.84-95UP")
  ifelse(as.character(mgmt)=="Organic", vik[as.character(yr)], pio[as.character(yr)])
}

# CHECK 1: rule must reproduce the (real) fungal Cultivar exactly
verify_against_fungi <- function(ps_fun, yr) {
  sd <- data.frame(sample_data(ps_fun))
  if (!"Cultivar" %in% colnames(sd)) stop(yr, ": fungal object has no Cultivar to verify against")
  derived <- cult_for(yr, sd$Management); actual <- as.character(sd$Cultivar)
  if (!all(derived == actual)) {
    bad <- which(derived != actual)[1]
    stop(sprintf("%s: backfill rule mismatch (mgmt=%s -> '%s' but data='%s').",
                 yr, sd$Management[bad], derived[bad], actual[bad]))
  }
  cat(sprintf("  CHECK 1 PASS (%s): rule reproduces fungal Cultivar exactly (%d samples)\n", yr, nrow(sd)))
}

# Backfill: skip if Cultivar already present, else patch
patch_cultivar <- function(ps, yr) {
  sd <- sample_data(ps); mg <- as.character(sd$Management)
  if ("Cultivar" %in% colnames(sd) && sum(!is.na(sd$Cultivar)) > 0) {
    cat(sprintf("  SKIP (%s bacteria): Cultivar already present\n", yr)); return(ps)
  }
  if (any(is.na(mg) | !(mg %in% c("Conventional","No-Till","Organic"))))
    stop(yr, " bacteria: missing/unexpected Management; cannot backfill safely.")
  new <- cult_for(yr, mg)
  if (any(is.na(new))) stop(yr, " bacteria: backfill produced NA cultivar.")
  sd$Cultivar <- new; sample_data(ps) <- sd
  cat(sprintf("  PATCHED (%s bacteria): %d samples backfilled\n", yr, nrow(sd)))
  ps
}

# CHECK 4: confounding must be clean
confound_ok <- function(ps, yr) {
  sd <- data.frame(sample_data(ps))
  pio_in_org <- sum(grepl("Pioneer", sd$Cultivar) & sd$Management=="Organic")
  vik_in_cnt <- sum(!grepl("Pioneer", sd$Cultivar) & sd$Management %in% c("Conventional","No-Till"))
  if (pio_in_org != 0 || vik_in_cnt != 0)
    stop(sprintf("%s: confounding not clean (Pioneer-in-Org=%d, Viking-in-Conv/NT=%d).", yr, pio_in_org, vik_in_cnt))
  cat(sprintf("  CHECK 4 PASS (%s bacteria): clean confounding\n", yr))
  print(table(sd$Cultivar, sd$Management))
}

cat("=== Pre-flight checks ===\n")
verify_against_fungi(physeq_rare_2018, "2018")
verify_against_fungi(physeq_rare_2020, "2020")
physeq_rare_bac_2018 <- patch_cultivar(physeq_rare_bac_2018, "2018")
physeq_rare_bac_2020 <- patch_cultivar(physeq_rare_bac_2020, "2020")
cat("\n=== Post-patch confounding ===\n")
confound_ok(physeq_rare_bac_2018, "2018")
confound_ok(physeq_rare_bac_2020, "2020")

# 1. VARPART — demonstrates cultivar & management are INSEPARABLE
#    (expected to fail with collinearity; the failure is the finding)
cat("\n=== Variance partitioning: Cultivar vs Management ===\n")
cat("    (perfect confounding -> varpart cannot partition; collinearity expected)\n")
run_varpart <- function(ps, yr, king) {
  for (comp in c("aboveground","belowground")) {
    md <- data.frame(sample_data(ps))
    keep <- as.character(md$Compartment) == comp
    if (sum(keep) < 10) next
    ps_c <- prune_samples(sample_names(ps)[keep], ps); ps_c <- prune_taxa(taxa_sums(ps_c) > 0, ps_c)
    md <- data.frame(sample_data(ps_c)); md$Cultivar <- factor(as.character(md$Cultivar))
    md$Management <- factor(as.character(md$Management))
    otu <- as(otu_table(ps_c),"matrix"); if (taxa_are_rows(ps_c)) otu <- t(otu)
    md <- md[rownames(otu), ]; D <- vegdist(otu, "bray")
    res <- tryCatch({
      vp <- varpart(D, ~ Cultivar, ~ Management, data = md)
      f <- vp$part$indfract$Adj.R.squared
      sprintf("    %s %s %s: Cultivar_alone=%.3f Shared=%.3f Mgmt_alone=%.3f", yr, king, comp, f[1], f[2], f[3])
    }, warning = function(w) sprintf("    %s %s %s: CONFOUNDED - %s", yr, king, comp, conditionMessage(w)),
    error   = function(e) sprintf("    %s %s %s: CONFOUNDED - cannot partition (%s)", yr, king, comp, conditionMessage(e)))
    cat(res, "\n")
  }
}
run_varpart(physeq_rare_2018,     "2018","Fungi")
run_varpart(physeq_rare_2020,     "2020","Fungi")
run_varpart(physeq_rare_bac_2018, "2018","Bacteria")
run_varpart(physeq_rare_bac_2020, "2020","Bacteria")
cat("    -> cultivar and management share identical information; effects inseparable by design.\n")

# 2. PERMANOVA — Conv vs No-Till WITHIN Pioneer (cultivar held constant)

cat("\n=== Conventional vs No-Till within Pioneer cultivar ===\n")
run_conv_vs_nt_pioneer <- function(ps, yr, king) {
  out <- list()
  for (comp in c("aboveground","belowground")) {
    md <- data.frame(sample_data(ps))
    keep <- as.character(md$Compartment)==comp & grepl("Pioneer", as.character(md$Cultivar)) &
      as.character(md$Management) %in% c("Conventional","No-Till")
    if (sum(keep) < 10) { cat(sprintf("  %s %s %s: n=%d <10, skip\n", yr,king,comp,sum(keep))); next }
    ps_c <- prune_samples(sample_names(ps)[keep], ps); ps_c <- prune_taxa(taxa_sums(ps_c)>0, ps_c)
    md <- data.frame(sample_data(ps_c)); md$Management <- factor(as.character(md$Management),
                                                                 levels=c("Conventional","No-Till"))
    otu <- as(otu_table(ps_c),"matrix"); if (taxa_are_rows(ps_c)) otu <- t(otu)
    md <- md[rownames(otu), ]; bc <- vegdist(otu,"bray"); set.seed(42)
    res <- adonis2(bc ~ Management, data=md, permutations=999, by="terms")
    pr <- as.data.frame(res)["Management", ]
    bd <- betadisper(bc, md$Management); bdt <- permutest(bd, permutations=999)
    out[[length(out)+1]] <- data.frame(Year=yr, Kingdom=king, Compartment=comp,
                                       Cultivar="Pioneer (fixed)", n=nrow(otu), R2=round(pr$R2,3), F=round(pr$F,2),
                                       p=round(pr$`Pr(>F)`,3), disp_p=round(bdt$tab$`Pr(>F)`[1],3), stringsAsFactors=FALSE)
  }
  out
}
cvnt <- do.call(rbind, c(
  run_conv_vs_nt_pioneer(physeq_rare_2018,     "2018","Fungi"),
  run_conv_vs_nt_pioneer(physeq_rare_2020,     "2020","Fungi"),
  run_conv_vs_nt_pioneer(physeq_rare_bac_2018, "2018","Bacteria"),
  run_conv_vs_nt_pioneer(physeq_rare_bac_2020, "2020","Bacteria")
))
print(cvnt, row.names=FALSE)
write.csv(cvnt, "PERMANOVA_ConvVsNT_withinPioneer.csv", row.names=FALSE)

# 3. PCoA FIGURE — Pioneer-only Conv vs No-Till, stats computed live per panel

pioneer_pcoa <- function(ps, yr, comp, panel_lab, title) {
  md <- data.frame(sample_data(ps))
  keep <- as.character(md$Compartment)==comp & grepl("Pioneer", as.character(md$Cultivar)) &
    as.character(md$Management) %in% c("Conventional","No-Till")
  ps_c <- prune_samples(sample_names(ps)[keep], ps); ps_c <- prune_taxa(taxa_sums(ps_c)>0, ps_c)
  otu <- as(otu_table(ps_c),"matrix"); if (taxa_are_rows(ps_c)) otu <- t(otu)
  meta <- data.frame(sample_data(ps_c)); meta$Management <- factor(as.character(meta$Management),
                                                                   levels=c("Conventional","No-Till"))
  meta <- meta[rownames(otu), ]; bc <- vegdist(otu,"bray")
  pcoa <- cmdscale(bc, k=2, eig=TRUE); ev <- pcoa$eig; pct <- round(100*ev[1:2]/sum(ev[ev>0]),1)
  df <- data.frame(PCo1=pcoa$points[,1], PCo2=pcoa$points[,2], Management=meta$Management)
  set.seed(42); res <- adonis2(bc ~ Management, data=meta, permutations=999, by="terms")
  r2 <- round(as.data.frame(res)["Management","R2"],3); pv <- as.data.frame(res)["Management","Pr(>F)"]
  pv_lab <- if (pv<0.001) "p < 0.001" else paste0("p = ", formatC(pv, format="f", digits=3))
  sub <- paste0("R\u00B2 = ", formatC(r2, format="f", digits=3)," | ", pv_lab," | n = ", nrow(df))
  ggplot(df, aes(PCo1, PCo2, color=Management)) +
    stat_ellipse(aes(group=Management), linetype="dashed", linewidth=0.4, show.legend=FALSE) +
    geom_point(size=1.3, alpha=0.85) + scale_color_manual(values=mgmt_cols) +
    labs(title=title, subtitle=sub, x=paste0("PCo1 (",pct[1],"%)"), y=paste0("PCo2 (",pct[2],"%)"), tag=panel_lab) +
    theme_bw(base_size=9) +
    theme(legend.position="none", plot.title=element_text(face="bold", size=9, hjust=0.5),
          plot.subtitle=element_text(size=6, hjust=0.5, colour="grey25"),
          plot.tag=element_text(face="bold", size=11), panel.grid.minor=element_blank())
}

pA <- pioneer_pcoa(physeq_rare_2018,     "2018","aboveground","A","Fungi | Aboveground")
pB <- pioneer_pcoa(physeq_rare_2018,     "2018","belowground","B","Fungi | Belowground")
pC <- pioneer_pcoa(physeq_rare_bac_2018, "2018","aboveground","C","Bacteria | Aboveground")
pD <- pioneer_pcoa(physeq_rare_bac_2018, "2018","belowground","D","Bacteria | Belowground")
pE <- pioneer_pcoa(physeq_rare_2020,     "2020","aboveground","E","Fungi | Aboveground")
pF <- pioneer_pcoa(physeq_rare_2020,     "2020","belowground","F","Fungi | Belowground")
pG <- pioneer_pcoa(physeq_rare_bac_2020, "2020","aboveground","G","Bacteria | Aboveground")
pH <- pioneer_pcoa(physeq_rare_bac_2020, "2020","belowground","H","Bacteria | Belowground")

leg_src <- pA + theme(legend.position="bottom") +
  guides(color=guide_legend(title="Management", override.aes=list(size=4))) +
  theme(legend.title=element_text(face="bold", size=11), legend.text=element_text(size=10))
g <- ggplotGrob(leg_src); legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name)=="guide-box")]]

row18 <- textGrob("2018\n(Soybean)", rot=90, gp=gpar(fontface="bold", fontsize=11))
row20 <- textGrob("2020\n(Maize)",   rot=90, gp=gpar(fontface="bold", fontsize=11))
title <- textGrob("Conventional vs. No-Till within Pioneer cultivar (cultivar held constant)",
                  gp=gpar(fontface="bold", fontsize=14))

grid.newpage()
grid.arrange(
  title,
  arrangeGrob(
    arrangeGrob(row18, pA, pB, pC, pD, ncol=5, widths=c(0.18,1,1,1,1)),
    arrangeGrob(row20, pE, pF, pG, pH, ncol=5, widths=c(0.18,1,1,1,1)),
    ncol=1),
  legend,
  heights=c(0.08, 1, 0.10)
)

cat("\n=== DONE ===\n")



####WITH BOLDING FOR SIG
library(phyloseq); library(vegan); library(ggplot2); library(grid); library(gridExtra)

pioneer_pcoa <- function(ps, yr, comp, panel_lab, title) {
  md <- data.frame(sample_data(ps))
  keep <- as.character(md$Compartment) == comp &
    grepl("Pioneer", as.character(md$Cultivar)) &
    as.character(md$Management) %in% c("Conventional", "No-Till")
  ps_c <- prune_samples(sample_names(ps)[keep], ps)
  ps_c <- prune_taxa(taxa_sums(ps_c) > 0, ps_c)
  otu <- as(otu_table(ps_c), "matrix"); if (taxa_are_rows(ps_c)) otu <- t(otu)
  meta <- data.frame(sample_data(ps_c))
  meta$Management <- factor(as.character(meta$Management), levels = c("Conventional", "No-Till"))
  meta <- meta[rownames(otu), ]
  bc <- vegdist(otu, "bray")
  pcoa <- cmdscale(bc, k = 2, eig = TRUE)
  ev <- pcoa$eig; pct <- round(100 * ev[1:2] / sum(ev[ev > 0]), 1)
  df <- data.frame(PCo1 = pcoa$points[, 1], PCo2 = pcoa$points[, 2], Management = meta$Management)
  
  set.seed(42)
  res <- adonis2(bc ~ Management, data = meta, permutations = 999, by = "terms")
  r2  <- round(as.data.frame(res)["Management", "R2"], 3)
  pv  <- as.data.frame(res)["Management", "Pr(>F)"]
  
  # Significance label
  sig_lab <- if (pv < 0.001) "***" else if (pv < 0.01) "**" else if (pv < 0.05) "*" else "ns"
  pv_lab  <- if (pv < 0.001) "p < 0.001" else paste0("p = ", formatC(pv, format = "f", digits = 3))
  
  # Subtitle text (asterisks indicate significance)
  sub <- paste0("R\u00B2 = ", formatC(r2, format = "f", digits = 3),
                " | ", pv_lab, " ", sig_lab,
                " | n = ", nrow(df))
  
  # Bold subtitle if significant, plain if not
  sub_face <- if (pv < 0.05) "bold" else "plain"
  sub_col  <- if (pv < 0.05) "black" else "grey40"
  
  ggplot(df, aes(PCo1, PCo2, color = Management)) +
    stat_ellipse(aes(group = Management), linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
    geom_point(size = 1.3, alpha = 0.85) +
    scale_color_manual(values = mgmt_cols) +
    labs(title = title, subtitle = sub,
         x = paste0("PCo1 (", pct[1], "%)"),
         y = paste0("PCo2 (", pct[2], "%)"),
         tag = panel_lab) +
    theme_bw(base_size = 9) +
    theme(legend.position = "none",
          plot.title    = element_text(face = "bold", size = 9, hjust = 0.5),
          plot.subtitle = element_text(face = sub_face, size = 6, hjust = 0.5, colour = sub_col),
          plot.tag      = element_text(face = "bold", size = 11),
          panel.grid.minor = element_blank())
}

pA <- pioneer_pcoa(physeq_rare_2018,     "2018", "aboveground", "A", "Fungi | Aboveground")
pB <- pioneer_pcoa(physeq_rare_2018,     "2018", "belowground", "B", "Fungi | Belowground")
pC <- pioneer_pcoa(physeq_rare_bac_2018, "2018", "aboveground", "C", "Bacteria | Aboveground")
pD <- pioneer_pcoa(physeq_rare_bac_2018, "2018", "belowground", "D", "Bacteria | Belowground")
pE <- pioneer_pcoa(physeq_rare_2020,     "2020", "aboveground", "E", "Fungi | Aboveground")
pF <- pioneer_pcoa(physeq_rare_2020,     "2020", "belowground", "F", "Fungi | Belowground")
pG <- pioneer_pcoa(physeq_rare_bac_2020, "2020", "aboveground", "G", "Bacteria | Aboveground")
pH <- pioneer_pcoa(physeq_rare_bac_2020, "2020", "belowground", "H", "Bacteria | Belowground")

leg_src <- pA + theme(legend.position = "bottom") +
  guides(color = guide_legend(title = "Management", override.aes = list(size = 4))) +
  theme(legend.title = element_text(face = "bold", size = 11),
        legend.text  = element_text(size = 10))
g <- ggplotGrob(leg_src)
legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]

row18 <- textGrob("2018\n(Soybean)", rot = 90, gp = gpar(fontface = "bold", fontsize = 11))
row20 <- textGrob("2020\n(Maize)",   rot = 90, gp = gpar(fontface = "bold", fontsize = 11))
title <- textGrob("Conventional vs. No-Till within Pioneer cultivar (cultivar held constant)",
                  gp = gpar(fontface = "bold", fontsize = 14))

grid.newpage()
grid.arrange(
  title,
  arrangeGrob(
    arrangeGrob(row18, pA, pB, pC, pD, ncol = 5, widths = c(0.18, 1, 1, 1, 1)),
    arrangeGrob(row20, pE, pF, pG, pH, ncol = 5, widths = c(0.18, 1, 1, 1, 1)),
    ncol = 1),
  legend,
  heights = c(0.08, 1, 0.10))

# =============================================================================
# =============================================================================

#  SECTION 10
#  Niche-level clustering
#
#  Supplemental Figures 3 and 4.

norm_origin <- function(x){
  v <- tolower(as.character(x))
  out <- rep(NA_character_, length(v))
  out[grepl("stem.*leaf|leaf.*stem", v)] <- "Stem+Leaf"    # combined FIRST
  out[is.na(out) & grepl("soil", v)] <- "Soil"
  out[is.na(out) & grepl("stem", v)] <- "Stem"
  out[is.na(out) & grepl("leaf", v)] <- "Leaf"
  out[is.na(out) & grepl("root", v)] <- "Root"
  factor(out, levels=origin_levels)
}

pcoa_small <- function(ps, color_by, pal, title=NULL, ocol=ORIGIN_COL){
  ph <- function() ggdraw()   # blank cell
  if (is.null(ps) || nsamples(ps) < 3) return(ph())
  ps <- tryCatch({
    p <- prune_taxa(taxa_sums(ps) > 0, ps)
    prune_samples(sample_sums(p) > 0, p)
  }, error=function(e) NULL)
  if (is.null(ps) || nsamples(ps) < 3 || ntaxa(ps) < 2) return(ph())
  
  ord <- tryCatch(ordinate(ps,"PCoA","bray"), error=function(e) NULL)
  if (is.null(ord)) return(ph())
  
  eig <- ord$values$Relative_eig
  d <- as.data.frame(ord$vectors[,1:2]); colnames(d) <- c("PCo1","PCo2")
  meta <- as(sample_data(ps),"data.frame")
  grp <- switch(color_by,
                "Management" = as.character(meta$Management),
                "Stage"      = gsub("Infloresence","Inflorescence",
                                    as.character(meta$Growth_Stage_Description), ignore.case=TRUE),
                "Origin"     = as.character(norm_origin(meta[[ocol]])))
  d$Group <- factor(grp, levels=names(pal))
  
  ggplot(d, aes(PCo1,PCo2,color=Group)) +
    geom_point(size=1.2, alpha=.75) +
    stat_ellipse(linewidth=.45, alpha=.5, na.rm=TRUE) +
    scale_color_manual(values=pal, drop=FALSE, name=NULL) +
    labs(title=title, x=sprintf("Ax1 %.0f%%",100*eig[1]), y=sprintf("Ax2 %.0f%%",100*eig[2])) +
    theme_classic(base_size=8) +
    theme(legend.position="none",
          plot.title=element_text(face="bold",size=8.5,hjust=.5),
          axis.title=element_text(face="bold",size=6.5,color="black"),
          axis.text =element_text(face="bold",size=6,  color="black"),
          axis.line=element_line(color="black",linewidth=.4),
          axis.ticks=element_line(color="black",linewidth=.4))
}

mk_legend <- function(pal,title){
  df <- data.frame(x=seq_along(pal), G=factor(names(pal),levels=names(pal)))
  get_legend(ggplot(df,aes(x,x,color=G))+geom_point(size=4)+
               scale_color_manual(values=pal,name=title)+
               theme(legend.position="bottom",
                     legend.title=element_text(face="bold",size=10),
                     legend.text =element_text(size=9)))
}

build_year_col <- function(yr, kingdom){
  suf  <- if (kingdom=="Bacteria") "_bac" else ""
  ps   <- merge_phyloseq(get(paste0("psr_",yr,"_a",suf)), get(paste0("psr_",yr,"_b",suf)))
  crop <- switch(yr,"2018"="Soybean","2019"="Wheat","2020"="Maize")
  
  origin_p <- pcoa_small(ps,"Origin",origin_colors,"All samples \u00b7 niche")
  
  niche_vec <- norm_origin(as(sample_data(ps),"data.frame")[[ORIGIN_COL]])
# keep only niches actually present this year (>=3 samples)
  present <- origin_levels[sapply(origin_levels,
                                  function(nm) sum(as.character(niche_vec)==nm, na.rm=TRUE) >= 3)]
  
  cells <- list()
  for (nm in present){
    keep <- which(as.character(niche_vec)==nm)
    ps_n <- prune_samples(sample_names(ps)[keep], ps)
    cells[[length(cells)+1]] <- pcoa_small(ps_n,"Management",management_colors, paste0(nm," \u00b7 Mgmt"))
    cells[[length(cells)+1]] <- pcoa_small(ps_n,"Stage",     stage_colors,      paste0(nm," \u00b7 Stage"))
  }
  grid_cells <- plot_grid(plotlist=cells, ncol=2, align="hv")
  n_rows <- length(present)
  hdr <- ggdraw()+draw_label(paste0(crop," (",yr,")"),fontface="bold",size=14)
  plot_grid(hdr, origin_p, grid_cells, ncol=1,
            rel_heights=c(.04, .36, .42*n_rows))
}

make_niche_figure <- function(kingdom){
  cols  <- lapply(years, build_year_col, kingdom=kingdom)
  body  <- plot_grid(plotlist=cols, ncol=length(years))
  legs  <- plot_grid(mk_legend(origin_colors,"Niche"),
                     mk_legend(management_colors,"Management"),
                     mk_legend(stage_colors,"Growth stage"), ncol=3)
  title <- ggdraw()+draw_label(
    paste0(kingdom," \u2014 community clustering by niche, management, and growth stage"),
    fontface="bold",size=16)
  plot_grid(title, body, legs, ncol=1, rel_heights=c(.03,1,.05))
}

# ---- build ----
fig_niche_fungi    <- make_niche_figure("Fungi")
fig_niche_bacteria <- make_niche_figure("Bacteria")
print(fig_niche_fungi)
print(fig_niche_bacteria)

# =============================================================================
# =============================================================================
#  SECTION 11
#  Management and growth stage PCoA and UpSet
#
#  Supplemental Figures 7-10.

# OUTPUTS
#   fig_mgmt_quad    combined: management PCoA on top, UpSet underneath
#   fig_stage_quad   combined: growth stage PCoA on top, UpSet underneath
#   fig_mgmt_pcoa    split -> Supplemental Figure 7
#   fig_mgmt_upset   split -> Supplemental Figure 8
#   fig_stage_pcoa   split -> Supplemental Figure 9
#   fig_stage_upset  split -> Supplemental Figure 10
#   mgmt_upset_counts.csv / stage_upset_counts.csv

# DA/UpSet uses UNRAREFIED whole-year objects and ERRORS if missing
# PCoA uses RAREFIED

years <- c("2018", "2019", "2020")

# ---- resolver: unrarefied required for DA, no silent fallback --------------
resolve_year_obj <- function(yr, kingdom, require_unrare = TRUE) {
  if (kingdom == "bac") {
    cands <- c(paste0("ps_unrare_", yr, "_bac"), paste0("ps_", substr(yr, 3, 4)))
    rare  <- paste0("physeq_rare_bac_", yr)
  } else {
    cands <- c(paste0("ps_unrare_", yr))
    rare  <- paste0("physeq_rare_", yr)
  }
  for (nm in cands) if (exists(nm)) return(get(nm))
  if (require_unrare)
    stop("UNRAREFIED object missing for ", kingdom, " ", yr,
         " (looked for: ", paste(cands, collapse = ", "),
         "). Rebuild ps_unrare_", yr, if (kingdom == "bac") "_bac" else "",
         " before running — do NOT fall back to ", rare, ".")
  if (exists(rare)) { message("WARNING: rarefied fallback for ", kingdom, " ", yr); return(get(rare)) }
  stop("No object for ", kingdom, " ", yr)
}

preflight_unrare <- function() {
  ok <- TRUE
  for (yr in years) {
    f <- exists(paste0("ps_unrare_", yr))
    b <- exists(paste0("ps_unrare_", yr, "_bac")) || exists(paste0("ps_", substr(yr, 3, 4)))
    cat(sprintf("  %s  fungi_unrare=%s  bac_unrare=%s\n", yr, f, b))
    ok <- ok && f && b
  }
  if (!ok) stop("Missing unrarefied objects above — rebuild before running.")
  cat("preflight OK: all unrarefied objects present.\n")
}
preflight_unrare()


# ---- shared helpers --------------------------------------------------------
blend_colors <- function(cv) {
  if (length(cv) == 1) return(cv)
  rc <- t(sapply(cv, function(x) col2rgb(x)[, 1])); a <- pmin(colMeans(rc) * 0.8, 255)
  rgb(a[1], a[2], a[3], maxColorValue = 255)
}

process_contrast <- function(contrast, label, dds_obj, ps_obj) {
  as.data.frame(results(dds_obj, contrast = contrast)) %>%
    rownames_to_column("ASV") %>% filter(!is.na(padj)) %>%
    mutate(Comparison = label,
           Enriched = case_when(padj < 0.05 & log2FoldChange >  1 ~ contrast[2],
                                padj < 0.05 & log2FoldChange < -1 ~ contrast[3], TRUE ~ "NS")) %>%
    left_join(as.data.frame(tax_table(ps_obj)) %>% rownames_to_column("ASV"), by = "ASV")
}

blank_plot <- function(msg = "No enriched ASVs") ggplot() +
  annotate("text", x = .5, y = .5, label = msg, size = 4, color = "gray50", fontface = "italic") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + theme_void()

get_present_labels <- function(df, sets, ints) {
  out <- logical(length(ints))
  for (i in seq_along(ints)) {
    cols <- ints[[i]]; oth <- setdiff(sets, cols)
    m <- Reduce(`&`, lapply(cols, function(c) df[[c]]))
    if (length(oth)) m <- m & Reduce(`&`, lapply(oth, function(c) !df[[c]]))
    out[i] <- any(m, na.rm = TRUE)
  }
  sapply(ints, function(x) paste(x, collapse = "&"))[out]
}

make_upset_with_sets <- function(df, sets, ints, queries_fn, set_cols) {
  present <- get_present_labels(df, sets, ints)
  p_up <- ComplexUpset::upset(
    as.data.frame(df), intersect = sets, intersections = ints,
    sort_intersections_by = NULL, n_intersections = length(ints),
    min_size = 0, set_sizes = FALSE, queries = queries_fn(present),
    base_annotations = list('Intersection size' = intersection_size(text = list(size = 3))),
    themes = upset_modify_themes(list(
      'intersections_matrix' = theme(panel.grid = element_blank(), panel.background = element_blank(),
                                     text = element_text(color = "black", size = 9, face = "bold")),
      'Intersection size' = theme(panel.grid.major.y = element_line(color = "gray90", linewidth = .3),
                                  panel.grid.major.x = element_blank(), panel.background = element_blank(),
                                  text = element_text(color = "black", size = 9, face = "bold"),
                                  axis.title.y = element_text(size = 9, face = "bold"),
                                  axis.text = element_text(color = "black", face = "bold", size = 9),
                                  axis.ticks = element_line(color = "black")))))
  sc <- data.frame(set = factor(sets, levels = rev(sets)),
                   count = sapply(sets, function(s) sum(df[[s]], na.rm = TRUE)))
  p_sets <- ggplot(sc, aes(count, set, fill = set)) + geom_col(width = .6) +
    geom_text(aes(label = count), hjust = -.2, size = 3, fontface = "bold", color = "black") +
    scale_fill_manual(values = set_cols, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, .3))) + labs(x = "Set Size") +
    theme_classic(base_size = 9) +
    theme(axis.title.y = element_blank(),
          axis.text.y  = element_text(size = 9, face = "bold", color = "black"),
          axis.title.x = element_text(size = 9, face = "bold", color = "black"),
          axis.text.x  = element_text(size = 9, face = "bold", color = "black"),
          axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black"),
          axis.line.y = element_blank(), axis.ticks.y = element_blank())
  plot_grid(p_sets, p_up, ncol = 2, rel_widths = c(1, 4))
}

pcoa_core <- function(ps, group_vec, pal) {
  ord <- ordinate(ps, method = "PCoA", distance = "bray"); eig <- ord$values$Relative_eig
  d <- as.data.frame(as.matrix(ord$vectors[, 1:2])); colnames(d) <- c("PCo1", "PCo2")
  d$Group <- factor(group_vec, levels = names(pal))
  p <- ggplot(d, aes(PCo1, PCo2, color = Group)) + geom_point(size = 2, alpha = .75) +
    stat_ellipse(aes(group = Group), linewidth = .6, alpha = .5) +
    scale_color_manual(values = pal, drop = FALSE) + theme_classic(base_size = 10) +
    labs(x = sprintf("Axis.1 (%.1f%%)", 100 * eig[1]), y = sprintf("Axis.2 (%.1f%%)", 100 * eig[2])) +
    theme(legend.position = "none",
          axis.title = element_text(size = 11, face = "bold", color = "black"),
          axis.text  = element_text(size = 10, face = "bold", color = "black"),
          axis.line  = element_line(color = "black", linewidth = .7),
          axis.ticks = element_line(color = "black", linewidth = .7))
  tryCatch(ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE, alpha = .4, size = 4),
           error = function(e) p,
           warning = function(w) tryCatch(
             suppressWarnings(ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE, alpha = .4, size = 4)),
             error = function(e2) p))
}


# ---- MANAGEMENT ------------------------------------------------------------
ensure_factors_mgmt <- function(ps) {
  md <- as.data.frame(sample_data(ps))
  md$Management <- factor(gsub("-", "_", as.character(md$Management)),
                          levels = c("Conventional", "No_Till", "Organic"))
  md$Tissue <- gsub("-", "", tolower(as.character(md$Compartment)))
  sample_data(ps) <- sample_data(md); prune_samples(!is.na(sample_data(ps)$Management), ps)
}

make_pcoa_mgmt <- function(ps, pal) pcoa_core(ps, as.character(sample_data(ps)$Management), pal)

fixed_sets_mgmt <- c("Conventional", "No-Till", "Organic")
fixed_ints_mgmt <- list("Conventional", "No-Till", "Organic",
                        c("Conventional", "No-Till"), c("Conventional", "Organic"), c("No-Till", "Organic"))
mgmt_icol <- c(
  "Conventional&No-Till" = blend_colors(management_colors[c("Conventional", "No-Till")]),
  "Conventional&Organic" = blend_colors(management_colors[c("Conventional", "Organic")]),
  "No-Till&Organic"      = blend_colors(management_colors[c("No-Till", "Organic")]))

queries_mgmt <- function(pl) {
  q <- list()
  add <- function(n, i, f) if (n %in% pl) q[[length(q) + 1]] <<- upset_query(intersect = i, fill = f)
  add("Conventional", "Conventional", management_colors["Conventional"])
  add("No-Till", "No-Till", management_colors["No-Till"])
  add("Organic", "Organic", management_colors["Organic"])
  add("Conventional&No-Till", c("Conventional", "No-Till"), mgmt_icol["Conventional&No-Till"])
  add("Conventional&Organic", c("Conventional", "Organic"), mgmt_icol["Conventional&Organic"])
  add("No-Till&Organic", c("No-Till", "Organic"), mgmt_icol["No-Till&Organic"])
  q
}

# shared DA engine — figures AND count tables both read this
mgmt_membership <- function(ps_year, tissue) {
  ps <- ensure_factors_mgmt(ps_year)
  ps_sub <- prune_samples(as.character(sample_data(ps)$Tissue) == tissue, ps)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub); ps_sub <- prune_samples(sample_sums(ps_sub) > 0, ps_sub)
  if (nsamples(ps_sub) < 3) return(NULL)
  tryCatch({
    dds <- phyloseq_to_deseq2(ps_sub, ~ Management) %>% DESeq(sfType = "poscounts", quiet = TRUE)
    has <- function(x) x %in% levels(droplevels(sample_data(ps_sub)$Management))
    comps <- list(); labs <- character(0)
    if (all(has(c("Organic","Conventional")))) { comps <- c(comps, list(c("Management","Organic","Conventional"))); labs <- c(labs,"OvC") }
    if (all(has(c("No_Till","Conventional")))) { comps <- c(comps, list(c("Management","No_Till","Conventional"))); labs <- c(labs,"NvC") }
    if (all(has(c("Organic","No_Till"))))      { comps <- c(comps, list(c("Management","Organic","No_Till")));      labs <- c(labs,"OvN") }
    if (length(comps) == 0) return(NULL)
    en <- bind_rows(map2(comps, labs, ~process_contrast(.x, .y, dds, ps_sub))) %>% filter(Enriched != "NS")
    if (nrow(en) == 0) return(NULL)
    en %>% group_by(ASV) %>% summarise(
      Conventional = any(Enriched == "Conventional"), `No-Till` = any(Enriched == "No_Till"),
      Organic = any(Enriched == "Organic"), .groups = "drop")
  }, error = function(e) NULL)
}

run_upset_mgmt2 <- function(ps_year, tissue) {
  m <- mgmt_membership(ps_year, tissue); if (is.null(m)) return(blank_plot())
  make_upset_with_sets(m, fixed_sets_mgmt, fixed_ints_mgmt, queries_mgmt, management_colors)
}


# ---- GROWTH STAGE ----------------------------------------------------------
stage_map <- c("Vegetative" = "Veg", "Inflorescence" = "Inf", "Reproductive" = "Rep")

ensure_factors_stage <- function(ps) {
  md <- as.data.frame(sample_data(ps))
  raw <- gsub("Infloresence", "Inflorescence", as.character(md$Growth_Stage_Description), ignore.case = TRUE)
  md$Stage  <- factor(stage_map[raw], levels = c("Veg", "Inf", "Rep"))
  md$Tissue <- gsub("-", "", tolower(as.character(md$Compartment)))
  sample_data(ps) <- sample_data(md); prune_samples(!is.na(sample_data(ps)$Stage), ps)
}

make_pcoa_stage <- function(ps, pal) {
  grp <- gsub("Infloresence", "Inflorescence",
              as.character(sample_data(ps)$Growth_Stage_Description), ignore.case = TRUE)
  pcoa_core(ps, grp, pal)
}

fixed_sets_stage <- c("Vegetative", "Inflorescence", "Reproductive")
fixed_ints_stage <- list("Vegetative", "Inflorescence", "Reproductive",
                         c("Vegetative","Inflorescence"), c("Vegetative","Reproductive"), c("Inflorescence","Reproductive"))
stage_icol <- c(
  "Vegetative&Inflorescence"   = blend_colors(stage_colors[c("Vegetative","Inflorescence")]),
  "Vegetative&Reproductive"    = blend_colors(stage_colors[c("Vegetative","Reproductive")]),
  "Inflorescence&Reproductive" = blend_colors(stage_colors[c("Inflorescence","Reproductive")]))

queries_stage <- function(pl) {
  q <- list()
  add <- function(n, i, f) if (n %in% pl) q[[length(q) + 1]] <<- upset_query(intersect = i, fill = f)
  add("Vegetative", "Vegetative", stage_colors["Vegetative"])
  add("Inflorescence", "Inflorescence", stage_colors["Inflorescence"])
  add("Reproductive", "Reproductive", stage_colors["Reproductive"])
  add("Vegetative&Inflorescence", c("Vegetative","Inflorescence"), stage_icol["Vegetative&Inflorescence"])
  add("Vegetative&Reproductive", c("Vegetative","Reproductive"), stage_icol["Vegetative&Reproductive"])
  add("Inflorescence&Reproductive", c("Inflorescence","Reproductive"), stage_icol["Inflorescence&Reproductive"])
  q
}

# per-cell n>=3 stage guard handles the 2019 aboveground gap
stage_membership <- function(ps_year, tissue) {
  ps <- ensure_factors_stage(ps_year)
  ps_sub <- prune_samples(as.character(sample_data(ps)$Tissue) == tissue, ps)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub); ps_sub <- prune_samples(sample_sums(ps_sub) > 0, ps_sub)
  if (nsamples(ps_sub) < 3) return(NULL)
  sn <- table(droplevels(sample_data(ps_sub)$Stage)); ok <- names(sn)[sn >= 3]
  if (length(ok) < 2) return(NULL)
  tryCatch({
    dds <- phyloseq_to_deseq2(ps_sub, ~ Stage) %>% DESeq(sfType = "poscounts", quiet = TRUE)
    pk <- function(a, b) all(c(a, b) %in% ok)
    comps <- list(); labs <- character(0)
    if (pk("Inf","Veg")) { comps <- c(comps, list(c("Stage","Inf","Veg"))); labs <- c(labs,"IvV") }
    if (pk("Rep","Inf")) { comps <- c(comps, list(c("Stage","Rep","Inf"))); labs <- c(labs,"RvI") }
    if (pk("Rep","Veg")) { comps <- c(comps, list(c("Stage","Rep","Veg"))); labs <- c(labs,"RvV") }
    if (length(comps) == 0) return(NULL)
    en <- bind_rows(map2(comps, labs, ~process_contrast(.x, .y, dds, ps_sub))) %>% filter(Enriched != "NS")
    if (nrow(en) == 0) return(NULL)
    en %>% group_by(ASV) %>% summarise(
      Vegetative = any(Enriched == "Veg"), Inflorescence = any(Enriched == "Inf"),
      Reproductive = any(Enriched == "Rep"), .groups = "drop")
  }, error = function(e) NULL)
}

run_upset_stage2 <- function(ps_year, tissue) {
  m <- stage_membership(ps_year, tissue); if (is.null(m)) return(blank_plot("No enriched ASVs"))
  make_upset_with_sets(m, fixed_sets_stage, fixed_ints_stage, queries_stage, stage_colors)
}

# THE SINGLE DA PASS — everything below reads M and S, nothing recomputes.

build_metric_lists <- function(metric) {
  pf <- list(); pb <- list(); uf <- list(); ub <- list()
  pal <- if (metric == "mgmt") management_colors else stage_colors
  mkp <- if (metric == "mgmt") make_pcoa_mgmt   else make_pcoa_stage
  mku <- if (metric == "mgmt") run_upset_mgmt2  else run_upset_stage2
  for (yr in years) {
    fy <- resolve_year_obj(yr, "fun")
    by <- resolve_year_obj(yr, "bac")
    for (suf in c("a", "b")) {
      tis <- ifelse(suf == "a", "aboveground", "belowground"); key <- paste0(yr, "_", tis)
      message(metric, "  ", key)
      pf[[key]] <- mkp(get(paste0("psr_", yr, "_", suf)), pal)
      pb[[key]] <- mkp(get(paste0("psr_", yr, "_", suf, "_bac")), pal)
      uf[[key]] <- mku(fy, tis); ub[[key]] <- mku(by, tis)
    }
  }
  list(pcoa_fungi = pf, pcoa_bac = pb, upset_fungi = uf, upset_bac = ub)
}

M <- build_metric_lists("mgmt")
S <- build_metric_lists("stage")


# ---- shared layout pieces --------------------------------------------------
crop_of   <- function(y) switch(y, "2018"="Soybean", "2019"="Wheat", "2020"="Maize")
mk_header <- function(t, sz = 26) ggdraw() + draw_label(t, fontface = "bold", size = sz)
mk_ylab   <- function(t, sz = 18) ggdraw() + draw_label(t, fontface = "bold", size = sz, angle = 90)

make_legend <- function(color_palette, legend_title, pt = 10, ts = 26, ls = 22, kw = 1.4) {
  ld <- data.frame(x = seq_along(color_palette),
                   Group = factor(names(color_palette), levels = names(color_palette)))
  lp <- ggplot(ld, aes(x, x, color = Group)) + geom_point(size = pt) +
    scale_color_manual(values = color_palette, name = legend_title) +
    guides(color = guide_legend(
      title.theme = element_text(size = ts, face = "bold"),
      label.theme = element_text(size = ls, face = "bold"),
      override.aes = list(size = pt), title.position = "left",
      direction = "horizontal", keywidth = unit(kw, "cm"), keyheight = unit(kw, "cm"))) +
    theme(legend.position = "top")
  get_legend(lp)
}

# one kingdom: 3 year-rows x 2 compartment-cols. off = panel-letter offset.
build_grid <- function(plist, off, hdr_sz = 20, lbl_sz = 22, ylab_sz = 18) {
  ch <- plot_grid(ggdraw(), mk_header("Aboveground", hdr_sz), mk_header("Belowground", hdr_sz),
                  ncol = 3, rel_widths = c(.12, 1, 1))
  ls <- LETTERS[(off + 1):(off + 6)]; rows <- list()
  for (i in seq_along(years)) {
    yr   <- years[i]
    body <- plot_grid(plist[[paste0(yr, "_aboveground")]],
                      plist[[paste0(yr, "_belowground")]],
                      ncol = 2, labels = ls[(2*i - 1):(2*i)],
                      label_size = lbl_sz, label_fontface = "bold")
    rows[[yr]] <- plot_grid(mk_ylab(paste0(crop_of(yr), "\n(", yr, ")"), ylab_sz), body,
                            ncol = 2, rel_widths = c(.12, 2))
  }
  plot_grid(ch, rows[["2018"]], rows[["2019"]], rows[["2020"]],
            ncol = 1, rel_heights = c(.10, 1, 1, 1))
}


# ---- COMBINED: PCoA on top, UpSet underneath. Letters A-X. -----------------
assemble_quad <- function(bundle, color_palette, legend_title) {
  leg      <- make_legend(color_palette, legend_title)
  king_row <- plot_grid(mk_header("Fungi", 28), mk_header("Bacteria", 28), ncol = 2)
  pcoa_blk <- plot_grid(mk_header("PCoA (Bray-Curtis)", 30), king_row,
                        plot_grid(build_grid(bundle$pcoa_fungi, 0),
                                  build_grid(bundle$pcoa_bac,   6), ncol = 2),
                        ncol = 1, rel_heights = c(.05, .05, 1))
  upst_blk <- plot_grid(mk_header("Differential abundance (UpSet)", 30), king_row,
                        plot_grid(build_grid(bundle$upset_fungi, 12),
                                  build_grid(bundle$upset_bac,   18), ncol = 2),
                        ncol = 1, rel_heights = c(.05, .05, 1))
  plot_grid(leg, pcoa_blk, upst_blk, ncol = 1, rel_heights = c(.04, 1, 1))
}


# ---- SPLIT: one content type per figure. Letters restart at A. -------------
assemble_single <- function(fungi_list, bac_list, color_palette, legend_title) {
  leg <- make_legend(color_palette, legend_title)
  plot_grid(leg,
            plot_grid(mk_header("Fungi", 28), mk_header("Bacteria", 28), ncol = 2),
            plot_grid(build_grid(fungi_list, 0), build_grid(bac_list, 6), ncol = 2),
            ncol = 1, rel_heights = c(.05, .05, 1))
}


# ---- BUILD EVERY FIGURE FROM M AND S ---------------------------------------
fig_mgmt_quad   <- assemble_quad(M, management_colors, "Management")
fig_stage_quad  <- assemble_quad(S, stage_colors,      "Growth Stage")

fig_mgmt_pcoa   <- assemble_single(M$pcoa_fungi,  M$pcoa_bac,  management_colors, "Management")
fig_mgmt_upset  <- assemble_single(M$upset_fungi, M$upset_bac, management_colors, "Management")
fig_stage_pcoa  <- assemble_single(S$pcoa_fungi,  S$pcoa_bac,  stage_colors,      "Growth Stage")
fig_stage_upset <- assemble_single(S$upset_fungi, S$upset_bac, stage_colors,      "Growth Stage")

print(fig_mgmt_quad)     # combined management
print(fig_stage_quad)    # combined growth stage
print(fig_mgmt_pcoa)     # Supplemental Figure 7
print(fig_mgmt_upset)    # Supplemental Figure 8
print(fig_stage_pcoa)    # Supplemental Figure 9
print(fig_stage_upset)   # Supplemental Figure 10

# UPSET COUNT TABLES
cell_counts <- function(m, sets) {
  combos <- list(sets[1], sets[2], sets[3],
                 c(sets[1], sets[2]), c(sets[1], sets[3]), c(sets[2], sets[3]),
                 c(sets[1], sets[2], sets[3]))
  nm <- c(paste(sets[1], "only"), paste(sets[2], "only"), paste(sets[3], "only"),
          paste(sets[1], "+", sets[2]), paste(sets[1], "+", sets[3]),
          paste(sets[2], "+", sets[3]), "All three")
  setNames(sapply(combos, function(isect) {
    others <- setdiff(sets, isect)
    inall  <- Reduce(`&`, lapply(isect, function(c) m[[c]]))
    inno   <- if (length(others)) Reduce(`&`, lapply(others, function(c) !m[[c]])) else TRUE
    sum(inall & inno, na.rm = TRUE)
  }), nm)
}

empty7 <- function(sets) setNames(rep(0L, 7),
                                  c(paste(sets[1], "only"), paste(sets[2], "only"), paste(sets[3], "only"),
                                    paste(sets[1], "+", sets[2]), paste(sets[1], "+", sets[3]),
                                    paste(sets[2], "+", sets[3]), "All three"))

export_counts <- function(kind) {
  sets <- if (kind == "mgmt") fixed_sets_mgmt else fixed_sets_stage
  memf <- if (kind == "mgmt") mgmt_membership else stage_membership
  rows <- list()
  for (yr in years) {
    fy <- resolve_year_obj(yr, "fun"); by <- resolve_year_obj(yr, "bac")
    for (suf in c("a", "b")) {
      tis <- ifelse(suf == "a", "aboveground", "belowground")
      for (km in c("Fungi", "Bacteria")) {
        ps  <- if (km == "Fungi") fy else by
        m   <- memf(ps, tis)
        cnt <- if (is.null(m)) empty7(sets) else cell_counts(m, sets)
        rows[[length(rows) + 1]] <- data.frame(Kingdom = km, Year = yr, Tissue = tis,
                                               as.list(cnt), TOTAL = sum(cnt), check.names = FALSE)
      }
    }
  }
  do.call(rbind, rows)
}

mgmt_counts_table  <- export_counts("mgmt")
stage_counts_table <- export_counts("stage")

cat("\n=========== MGMT UPSET COUNTS (UNRAREFIED) ===========\n")
print(mgmt_counts_table, row.names = FALSE)
cat("\n=========== STAGE UPSET COUNTS (UNRAREFIED) ===========\n")
print(stage_counts_table, row.names = FALSE)

write.csv(mgmt_counts_table,  "mgmt_upset_counts.csv",  row.names = FALSE)
write.csv(stage_counts_table, "stage_upset_counts.csv", row.names = FALSE)
cat("\nwrote mgmt_upset_counts.csv and stage_upset_counts.csv\n")

# =============================================================================
# =============================================================================
#  SECTION 12
#  Cross-kingdom co-occurrence networks
#
#  Table 1 and Supplemental Table 12.
#  Writes hub_taxa_YrCompMgmt.csv, read by section 15.

# Year × Compartment × Management. Bacteria >=50%, Fungi >=30%.
# NO CAP, FDR 0.05. Hubs: degree + betweenness, top 10% on BOTH.
# Hub table + edges carry DNA Sequence (cross-year identity = sequence, not id).

# ---- ensure per-year unrarefied objects exist ------------------------

if (!exists("ps_unrare_2018_bac") && exists("ps_18")) ps_unrare_2018_bac <- ps_18
if (!exists("ps_unrare_2019_bac") && exists("ps_19")) ps_unrare_2019_bac <- ps_19
if (!exists("ps_unrare_2020_bac") && exists("ps_20")) ps_unrare_2020_bac <- ps_20


.fun_missing <- !all(sapply(c("ps_unrare_2018","ps_unrare_2019","ps_unrare_2020"), exists))
if (.fun_missing) {
  if (exists("process_year_fun") && all(sapply(c("main_2018","main_2020","ps_2019"), exists))) {
    cat("Rebuilding fungal ps_unrare_* (process_year_fun -> list(rare, unrare))...\n")
    body(process_year_fun)[[length(body(process_year_fun))]] <-
      quote(list(rare = ps_rare, unrare = ps))
    .t <- process_year_fun(main_2018, "2018"); physeq_rare_2018 <- .t$rare; ps_unrare_2018 <- .t$unrare
    .t <- process_year_fun(main_2020, "2020"); physeq_rare_2020 <- .t$rare; ps_unrare_2020 <- .t$unrare
    .t <- process_year_fun(ps_2019,   "2019"); physeq_rare_2019 <- .t$rare; ps_unrare_2019 <- .t$unrare
  } else {
    stop("Fungal ps_unrare_* missing AND can't rebuild (need process_year_fun + main_2018/main_2020/ps_2019).\n",
         " -> run the fungal stage first, or source the edited stage file.")
  }
}

.need <- c("ps_unrare_2018","ps_unrare_2019","ps_unrare_2020",
           "ps_unrare_2018_bac","ps_unrare_2019_bac","ps_unrare_2020_bac")
.miss <- .need[!sapply(.need, exists)]
if (length(.miss) > 0) stop("Still missing: ", paste(.miss, collapse = ", "))

get_bac <- function(yr) get(paste0("ps_unrare_", yr, "_bac"))
get_fun <- function(yr) get(paste0("ps_unrare_", yr))


# FUNCTIONS

build_correlation_network <- function(ps_bac, ps_fun, cor_threshold = 0.3, fdr_threshold = 0.05) {
  sd_bac_sub <- as.data.frame(sample_data(ps_bac))
  sd_fun_sub <- as.data.frame(sample_data(ps_fun))
  bac_lookup <- setNames(rownames(sd_bac_sub), sd_bac_sub$CoOccurence)
  fun_lookup <- setNames(rownames(sd_fun_sub), sd_fun_sub$CoOccurence)
  common <- sort(intersect(names(bac_lookup), names(fun_lookup)))
  cat("  n =", length(common), "samples\n")
  
  otu_bac <- as(otu_table(ps_bac), "matrix"); if (taxa_are_rows(ps_bac)) otu_bac <- t(otu_bac)
  otu_fun <- as(otu_table(ps_fun), "matrix"); if (taxa_are_rows(ps_fun)) otu_fun <- t(otu_fun)
  otu_bac_aligned <- otu_bac[bac_lookup[common], , drop = FALSE]
  otu_fun_aligned <- otu_fun[fun_lookup[common], , drop = FALSE]
  
  clr_transform <- function(x) { x[x == 0] <- 0.5; log(x) - mean(log(x)) }
  otu_bac_clr <- t(apply(otu_bac_aligned, 1, clr_transform))
  otu_fun_clr <- t(apply(otu_fun_aligned, 1, clr_transform))
  n_bac <- ncol(otu_bac_clr); n_fun <- ncol(otu_fun_clr)
  cat("  Calculating", n_bac, "x", n_fun, "=", n_bac * n_fun, "correlations...\n")
  
  cor_mat <- matrix(NA, n_bac, n_fun); p_mat <- matrix(NA, n_bac, n_fun)
  for (i in 1:n_bac) for (j in 1:n_fun) {
    test <- cor.test(otu_bac_clr[, i], otu_fun_clr[, j], method = "spearman", exact = FALSE)
    cor_mat[i, j] <- test$estimate; p_mat[i, j] <- test$p.value
  }
  p_adj <- matrix(p.adjust(p_mat, method = "fdr"), n_bac, n_fun)
  sig_edges <- which(abs(cor_mat) >= cor_threshold & p_adj < fdr_threshold, arr.ind = TRUE)
  if (nrow(sig_edges) == 0) return(NULL)
  
  edges <- data.frame(bac_idx = sig_edges[, 1], fun_idx = sig_edges[, 2],
                      correlation = cor_mat[sig_edges], p_value = p_mat[sig_edges], p_adj = p_adj[sig_edges])
  list(edges = edges, bac_names = colnames(otu_bac_aligned), fun_names = colnames(otu_fun_aligned),
       bac_refseq = tryCatch(refseq(ps_bac), error = function(e) NULL),
       fun_refseq = tryCatch(refseq(ps_fun), error = function(e) NULL))
}

identify_hubs <- function(ig, q = 0.90) {
  deg   <- degree(ig)
  btw_w <- if (!is.null(E(ig)$weight)) 1 / E(ig)$weight else NA
  btw   <- betweenness(ig, weights = btw_w, directed = FALSE)
  hubs  <- names(deg)[deg >= quantile(deg, q) & btw >= quantile(btw, q)]
  list(degree = deg, betweenness = btw, hubs = hubs)
}

# MAIN

BAC_PREV <- 0.5; FUN_PREV <- 0.3; COR_THRESHOLD <- 0.3; FDR_THRESHOLD <- 0.05
years <- c("2018","2019","2020"); compartments <- c("aboveground","belowground")
managements <- c("Conventional","No-Till","Organic")

all_network_results <- list(); all_hub_taxa <- list(); network_summary <- data.frame()

for (yr in years) {
  ps_bac_year <- get_bac(yr); ps_fun_year <- get_fun(yr)
  for (comp in compartments) {
    for (mgmt in managements) {
      cat("\n========================================\n")
      cat("Processing", yr, comp, mgmt, "\n")
      cat("========================================\n")
      
      sd_bac <- as.data.frame(sample_data(ps_bac_year))
      sd_fun <- as.data.frame(sample_data(ps_fun_year))
      bac_subset <- sd_bac[sd_bac$Compartment == comp & sd_bac$Management == mgmt & !is.na(sd_bac$CoOccurence), ]
      fun_subset <- sd_fun[sd_fun$Compartment == comp & sd_fun$Management == mgmt & !is.na(sd_fun$CoOccurence), ]
      common_cooccur <- intersect(bac_subset$CoOccurence, fun_subset$CoOccurence)
      cat("Paired samples:", length(common_cooccur), "\n")
      if (length(common_cooccur) < 20) { cat("Skipping - insufficient samples\n"); next }
      
      bac_samples <- rownames(bac_subset)[bac_subset$CoOccurence %in% common_cooccur]
      fun_samples <- rownames(fun_subset)[fun_subset$CoOccurence %in% common_cooccur]
      ps_bac_sub <- prune_taxa(taxa_sums(prune_samples(bac_samples, ps_bac_year)) > 0,
                               prune_samples(bac_samples, ps_bac_year))
      ps_fun_sub <- prune_taxa(taxa_sums(prune_samples(fun_samples, ps_fun_year)) > 0,
                               prune_samples(fun_samples, ps_fun_year))
      ps_bac_filt <- filter_taxa(ps_bac_sub, function(x) sum(x > 0) > (BAC_PREV * length(x)), TRUE)
      ps_fun_filt <- filter_taxa(ps_fun_sub, function(x) sum(x > 0) > (FUN_PREV * length(x)), TRUE)
      cat("After filtering:", ntaxa(ps_bac_filt), "bacteria,", ntaxa(ps_fun_filt), "fungi\n")
      if (ntaxa(ps_bac_filt) < 10 | ntaxa(ps_fun_filt) < 10) { cat("Skipping - too few taxa\n"); next }
      
      cat("Building network...\n")
      net_result <- build_correlation_network(ps_bac_filt, ps_fun_filt, COR_THRESHOLD, FDR_THRESHOLD)
      if (is.null(net_result)) { cat("No significant correlations\n"); next }
      
      bac_tax <- as.data.frame(tax_table(ps_bac_filt))
      fun_tax <- as.data.frame(tax_table(ps_fun_filt))
      edges_df <- net_result$edges
      edges_df$from_asv <- net_result$bac_names[edges_df$bac_idx]
      edges_df$to_asv   <- net_result$fun_names[edges_df$fun_idx]
      edges_df$from_kingdom <- "Bacteria"; edges_df$to_kingdom <- "Fungi"
      edges_df$from_phylum <- bac_tax[edges_df$from_asv, "Phylum"]
      edges_df$from_family <- bac_tax[edges_df$from_asv, "Family"]
      edges_df$from_genus  <- bac_tax[edges_df$from_asv, "Genus"]
      edges_df$to_phylum <- fun_tax[edges_df$to_asv, "Phylum"]
      edges_df$to_family <- fun_tax[edges_df$to_asv, "Family"]
      edges_df$to_genus  <- fun_tax[edges_df$to_asv, "Genus"]
      brs <- net_result$bac_refseq; frs <- net_result$fun_refseq
      edges_df$from_seq <- if (!is.null(brs)) as.character(brs[edges_df$from_asv]) else NA_character_
      edges_df$to_seq   <- if (!is.null(frs)) as.character(frs[edges_df$to_asv])   else NA_character_
      edges_df$Year <- yr; edges_df$Compartment <- comp; edges_df$Management <- mgmt
      edges_df$weight <- edges_df$correlation
      
      result_name <- paste(yr, comp, mgmt, sep = "_")
      all_network_results[[result_name]] <- edges_df
      
      ig <- simplify(graph_from_data_frame(edges_df[, c("from_asv","to_asv","weight")], directed = FALSE),
                     remove.multiple = TRUE, remove.loops = TRUE)
      E(ig)$weight <- abs(E(ig)$weight)
      modules <- cluster_fast_greedy(ig, weights = E(ig)$weight)
      mod <- modularity(modules); n_mod <- max(membership(modules))
      
      hubinfo <- identify_hubs(ig, q = 0.90)
      deg <- hubinfo$degree; btw <- hubinfo$betweenness; hub_names <- hubinfo$hubs
      bac_asvs <- unique(edges_df$from_asv); fun_asvs <- unique(edges_df$to_asv)
      hub_kingdoms <- ifelse(hub_names %in% bac_asvs, "Bacteria",
                             ifelse(hub_names %in% fun_asvs, "Fungi", "Unknown"))
      
      hub_seq <- mapply(function(asv, kg) {
        rs <- if (kg == "Bacteria") brs else if (kg == "Fungi") frs else NULL
        if (is.null(rs) || !asv %in% names(rs)) return(NA_character_)
        as.character(rs[asv])
      }, hub_names, hub_kingdoms, USE.NAMES = FALSE)
      
      bac_tax_lookup <- edges_df %>% distinct(from_asv, .keep_all = TRUE) %>%
        select(ASV = from_asv, Phylum = from_phylum, Family = from_family, Genus = from_genus)
      fun_tax_lookup <- edges_df %>% distinct(to_asv, .keep_all = TRUE) %>%
        select(ASV = to_asv, Phylum = to_phylum, Family = to_family, Genus = to_genus)
      tax_lookup <- bind_rows(bac_tax_lookup, fun_tax_lookup) %>% distinct(ASV, .keep_all = TRUE)
      
      hub_df <- data.frame(
        ASV = hub_names, Kingdom = hub_kingdoms, Sequence = hub_seq,
        Degree = deg[hub_names], Betweenness = round(btw[hub_names], 2),
        Module = membership(modules)[hub_names],
        Year = yr, Compartment = comp, Management = mgmt, stringsAsFactors = FALSE
      ) %>% left_join(tax_lookup, by = "ASV")
      all_hub_taxa[[result_name]] <- hub_df
      
      cat("  Edges:", nrow(edges_df), "| Modularity:", round(mod, 3),
          "| Modules:", n_mod, "| Avg degree:", round(mean(deg), 2), "\n")
      cat("  Hubs:", nrow(hub_df), "(", sum(hub_df$Kingdom == "Bacteria"), "bac,",
          sum(hub_df$Kingdom == "Fungi"), "fun )\n")
      if (nrow(hub_df) > 0) {
        hub_print <- hub_df %>% arrange(Kingdom, desc(Degree)) %>%
          mutate(label = paste0("    ", Kingdom, " | ", ASV, " | ", Genus, " (", Family, ")",
                                " | deg=", Degree, " | btw=", Betweenness, " | mod=", Module))
        cat(paste(hub_print$label, collapse = "\n"), "\n")
      }
      
      network_summary <- rbind(network_summary, data.frame(
        Year = yr, Compartment = comp, Management = mgmt,
        n_paired_samples = length(common_cooccur),
        n_bacteria = ntaxa(ps_bac_filt), n_fungi = ntaxa(ps_fun_filt),
        n_edges = nrow(edges_df), n_positive = sum(edges_df$correlation > 0),
        n_negative = sum(edges_df$correlation < 0),
        prop_negative = round(sum(edges_df$correlation < 0) / nrow(edges_df), 3),
        mean_abs_correlation = round(mean(abs(edges_df$correlation)), 3),
        modularity = round(mod, 3), n_modules = n_mod, avg_degree = round(mean(deg), 2),
        density = round(edge_density(ig), 4), n_hubs = nrow(hub_df),
        hub_bac = sum(hub_df$Kingdom == "Bacteria"), hub_fun = sum(hub_df$Kingdom == "Fungi"),
        stringsAsFactors = FALSE))
    }
  }
}


# COMPILE + SAVE
hub_all <- bind_rows(all_hub_taxa)
write.csv(network_summary, "cross_kingdom_network_summary.csv", row.names = FALSE)
if (length(all_network_results) > 0) {
  all_edges <- bind_rows(all_network_results)
  write.csv(all_edges, "cross_kingdom_all_edges.csv", row.names = FALSE)
}
write.csv(hub_all, "hub_taxa_YrCompMgmt.csv", row.names = FALSE)

cat("\n========================================\nANALYSIS COMPLETE\n========================================\n")
cat("Networks built:", nrow(network_summary), "\n")
cat("Total edges:", nrow(all_edges), "\n")
cat("Total hubs:", nrow(hub_all), "|", sum(!is.na(hub_all$Sequence)), "with sequence\n")
print(network_summary)

# =============================================================================
# =============================================================================
#  SECTION 13
#  Random forest
#
#  Figure 4, Supplemental Figure 11, Supplemental Tables 10 and 11.
#  Writes Supplementary_Table_RF2_ClassifierTaxa.csv, read by section 15.
#  Redefines `compartments` with a third level, "pooled".
# RANDOM FOREST: Management Prediction from Microbiome
# 27 models: 3 years × 3 kingdoms × 3 compartment sets
# Full pipeline: RFE → hyperparameter tuning → final model

RENDER_DIAGNOSTICS <- TRUE

# CUSTOM RANGER RFE FUNCTIONS
# Avoids "importance matched by multiple arguments" clash between caret's rfFuncs and ranger's importance parameter

rangerFuncs <- list(
  summary = defaultSummary,
  
  fit = function(x, y, first, last, ...) {
    df_tmp   <- as.data.frame(x)
    df_tmp$y <- y
    ranger(
      y ~ .,
      data          = df_tmp,
      num.trees     = 500,
      importance    = "permutation",
      min.node.size = 5,
      num.threads   = 1
    )
  },
  
  pred = function(object, x) {
    predict(object, data = as.data.frame(x))$predictions
  },
  
  rank = function(object, x, y) {
    imp <- object$variable.importance
    data.frame(
      var       = names(imp),
      Overall   = as.numeric(imp),
      row.names = names(imp)
    )
  },
  
  selectSize = pickSizeBest,
  selectVar  = pickVars
)

# HELPER: CLR-transform OTU table from phyloseq

get_clr_matrix <- function(ps) {
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu <- t(otu)
  otu <- matrix(as.numeric(otu), nrow = nrow(otu), ncol = ncol(otu),
                dimnames = dimnames(otu))
  otu[otu == 0] <- 0.5
  as.matrix(compositions::clr(otu))
}


# HELPER: Prefix taxa names

prefix_taxa <- function(mat, prefix) {
  colnames(mat) <- paste0(prefix, "_", colnames(mat))
  mat
}


# HELPER: Keep top N most variable features (by CLR variance)

filter_top_var <- function(mat, top_n = 500) {
  vars    <- apply(mat, 2, var)
  top_idx <- order(vars, decreasing = TRUE)[1:min(top_n, ncol(mat))]
  mat[, top_idx, drop = FALSE]
}


# HELPER: Safe Management vector extraction

get_mgmt <- function(ps, sample_order) {
  sd <- as.data.frame(sample_data(ps))
  make.names(as.character(sd[["Management"]][match(sample_order, rownames(sd))]))
}

# HELPER: Get phyloseq objects for a given year/compartment

get_ps_objects <- function(yr, cs) {
  if (cs == "all") {
    ps_b_name <- paste0("psr_", yr, "_bac")
    ps_f_name <- paste0("psr_", yr)
  } else {
    ps_b_name <- paste0("psr_", yr, "_", cs, "_bac")
    ps_f_name <- paste0("psr_", yr, "_", cs)
  }
  if (!exists(ps_b_name) || !exists(ps_f_name)) return(NULL)
  list(ps_bac = get(ps_b_name), ps_fun = get(ps_f_name))
}

# HELPER: Build feature matrix + Management label

build_feature_df <- function(ps_bac = NULL, ps_fun = NULL, kingdom = "combined",
                             top_var = 500) {
  
  if (kingdom == "bac") {
    feat <- filter_top_var(prefix_taxa(get_clr_matrix(ps_bac), "bac"), top_var)
    df   <- as.data.frame(feat)
    df$Management <- factor(get_mgmt(ps_bac, rownames(feat)), levels = mgmt_levels)
    
  } else if (kingdom == "fun") {
    feat <- filter_top_var(prefix_taxa(get_clr_matrix(ps_fun), "fun"), top_var)
    df   <- as.data.frame(feat)
    df$Management <- factor(get_mgmt(ps_fun, rownames(feat)), levels = mgmt_levels)
    
  } else {
    mat_bac <- prefix_taxa(get_clr_matrix(ps_bac), "bac")
    mat_fun <- prefix_taxa(get_clr_matrix(ps_fun), "fun")
    
    sd_bac <- as.data.frame(as.matrix(sample_data(ps_bac)))
    sd_fun <- as.data.frame(as.matrix(sample_data(ps_fun)))
    sd_bac$sample_name_bac <- rownames(sd_bac)
    sd_fun$sample_name_fun <- rownames(sd_fun)
    
    sd_bac <- sd_bac[!is.na(sd_bac$CoOccurence) & sd_bac$CoOccurence != "", ]
    sd_fun <- sd_fun[!is.na(sd_fun$CoOccurence) & sd_fun$CoOccurence != "", ]
    
    matched <- inner_join(
      sd_bac[, c("sample_name_bac", "CoOccurence", "Management")],
      sd_fun[, c("sample_name_fun", "CoOccurence")],
      by = "CoOccurence"
    )
    matched <- distinct(matched, CoOccurence, .keep_all = TRUE)
    
    if (nrow(matched) == 0) stop("No shared CoOccurence values")
    cat("  Matched", nrow(matched), "sample pairs via CoOccurence\n")
    
    feat_bac <- mat_bac[matched$sample_name_bac, , drop = FALSE]
    feat_fun <- mat_fun[matched$sample_name_fun, , drop = FALSE]
    rownames(feat_bac) <- matched$CoOccurence
    rownames(feat_fun) <- matched$CoOccurence
    
    feat_bac <- filter_top_var(feat_bac, top_n = top_var)
    feat_fun <- filter_top_var(feat_fun, top_n = top_var)
    
    df <- as.data.frame(cbind(feat_bac, feat_fun))
    df$Management <- factor(
      make.names(as.character(matched$Management)),
      levels = mgmt_levels
    )
  }
  
  df <- df[!is.na(df$Management), ]
  return(df)
}

make_tax_labels <- function(ps, prefix) {
  if (is.null(ps)) return(NULL)
  tax     <- as.data.frame(tax_table(ps))
  tax$ASV <- paste0(prefix, "_", rownames(tax))
  
  kingdom_tag <- ifelse(prefix == "bac", "[B] ", ifelse(prefix == "fun", "[F] ", ""))
  
  clean <- function(x) {
    x <- sub("_[0-9.]+$", "", x)
    x <- gsub("_", " ", x)
    x[is.na(x) | x == ""] <- NA
    # Null out plant/eukaryote SINTAX misassignments
    bad <- grepl("Triticum|Bromus|Zea mays|eukaryote|archaeote",
                 x, ignore.case = TRUE)
    x[bad] <- NA
    x
  }
  
  # Skip SINTAX columns entirely — use RDP first, then BLAST, then plain SILVA
  genus <- if ("Genus_RDP"    %in% colnames(tax)) clean(tax$Genus_RDP)
  else if ("Genus_BLAST"  %in% colnames(tax)) clean(tax$Genus_BLAST)
  else clean(tax$Genus)
  fam   <- if ("Family_RDP"   %in% colnames(tax)) clean(tax$Family_RDP)
  else if ("Family_BLAST" %in% colnames(tax)) clean(tax$Family_BLAST)
  else clean(tax$Family)
  ord   <- if ("Order_RDP"    %in% colnames(tax)) clean(tax$Order_RDP)
  else if ("Order_BLAST"  %in% colnames(tax)) clean(tax$Order_BLAST)
  else clean(tax$Order)
  
  best_name <- ifelse(!is.na(genus), genus,
                      ifelse(!is.na(fam), fam,
                             ifelse(!is.na(ord), ord, NA)))
  
  tax$label <- ifelse(
    !is.na(best_name),
    paste0(kingdom_tag, best_name, " (", rownames(tax), ")"),
    paste0(kingdom_tag, rownames(tax))
  )
  
  tax[, c("ASV", "label")]
}


# HELPER: RFE — single threaded

run_rfe <- function(df, seed = 42) {
  set.seed(seed)
  
  n_feat <- ncol(df) - 1
  sizes  <- unique(round(exp(seq(log(10), log(n_feat), length.out = 8))))
  sizes  <- sort(sizes[sizes <= n_feat], decreasing = TRUE)
  
  rfe_ctrl <- rfeControl(
    functions     = rangerFuncs,
    method        = "cv",
    number        = 5,
    allowParallel = FALSE
  )
  
  rfe_result <- rfe(
    x          = df[, -ncol(df), drop = FALSE],
    y          = df$Management,
    sizes      = sizes,
    rfeControl = rfe_ctrl
  )
  
  cat("  RFE optimal features:", rfe_result$optsize, "\n")
  cat("  RFE best accuracy:   ",
      round(max(rfe_result$results$Accuracy) * 100, 1), "%\n")
  
  optimal_vars <- predictors(rfe_result)
  df_rfe       <- df[, c(optimal_vars, "Management"), drop = FALSE]
  
  return(list(
    rfe_result   = rfe_result,
    df_rfe       = df_rfe,
    optimal_vars = optimal_vars
  ))
}


# HELPER: Full RF with hyperparameter tuning + CV

run_rf_cv <- function(df, seed = 42) {
  set.seed(seed)
  
  ctrl <- trainControl(
    method          = "repeatedcv",
    number          = 5,
    repeats         = 5,
    savePredictions = "final",
    classProbs      = TRUE,
    summaryFunction = multiClassSummary,
    allowParallel   = FALSE
  )
  
  n_features    <- ncol(df) - 1
  mtry_vals     <- unique(c(floor(sqrt(n_features)),
                            floor(sqrt(n_features) * 2),
                            floor(n_features / 3)))
  nodesize_vals <- c(1, 5, 10)
  ntree_vals    <- c(500, 1000, 2000)
  
  tune_grid <- expand.grid(
    mtry          = mtry_vals,
    min.node.size = nodesize_vals,
    splitrule     = "gini"
  )
  
  best_model  <- NULL
  best_acc    <- -Inf
  best_ntrees <- 500
  
  for (nt in ntree_vals) {
    set.seed(seed)
    m <- tryCatch(
      train(
        Management ~ .,
        data        = df,
        method      = "ranger",
        trControl   = ctrl,
        tuneGrid    = tune_grid,
        metric      = "Accuracy",
        importance  = "permutation",
        num.trees   = nt,
        num.threads = 1
      ),
      error = function(e) NULL
    )
    if (is.null(m)) next
    acc <- max(m$results$Accuracy)
    if (acc > best_acc) {
      best_acc    <- acc
      best_model  <- m
      best_ntrees <- nt
    }
  }
  
  if (is.null(best_model)) stop("All num.trees configurations failed")
  cat("  Best num.trees:", best_ntrees, "\n")
  cat("  Best mtry:", best_model$bestTune$mtry, "\n")
  cat("  Best min.node.size:", best_model$bestTune$min.node.size, "\n")
  
  # Confusion matrix from held-out CV predictions
  cv_preds <- best_model$pred
  cv_preds <- cv_preds[
    cv_preds$mtry          == best_model$bestTune$mtry &
      cv_preds$min.node.size == best_model$bestTune$min.node.size, ]
  
  actual    <- name_map[as.character(cv_preds$obs)]
  predicted <- name_map[as.character(cv_preds$pred)]
  conf_tbl       <- as.data.frame(table(Actual = actual, Predicted = predicted))
  conf_tbl$prop  <- ave(conf_tbl$Freq, conf_tbl$Actual, FUN = function(x) x / sum(x))
  conf_tbl$n     <- conf_tbl$Freq
  conf_tbl$Freq  <- NULL
  
  # Final ranger for importance
  set.seed(seed)
  model_ranger <- ranger(
    formula       = Management ~ .,
    data          = df,
    num.trees     = best_ntrees,
    mtry          = best_model$bestTune$mtry,
    min.node.size = best_model$bestTune$min.node.size,
    importance    = "permutation",
    num.threads   = 1
  )
  
  imp_df <- data.frame(
    ASV        = names(model_ranger$variable.importance),
    Importance = as.numeric(model_ranger$variable.importance),
    stringsAsFactors = FALSE
  )
  
  return(list(
    model_cv     = best_model,
    model_ranger = model_ranger,
    conf_df      = conf_tbl,
    imp_df       = imp_df,
    cv_accuracy  = best_acc,
    best_ntrees  = best_ntrees
  ))
}


# HELPER: Per-class importance (one-vs-rest ranger)

get_perclass_importance <- function(df, seed = 42, best_mtry = NULL,
                                    best_nodesize = 5, best_ntrees = 1000) {
  if (is.null(best_mtry)) best_mtry <- floor(sqrt(ncol(df) - 1))
  
  imp_list <- lapply(levels(df$Management), function(cls) {
    df_bin <- df
    df_bin$Management <- factor(
      ifelse(df$Management == cls, "target", "other"),
      levels = c("target", "other")
    )
    set.seed(seed)
    m <- ranger(
      formula       = Management ~ .,
      data          = df_bin,
      num.trees     = best_ntrees,
      mtry          = best_mtry,
      min.node.size = best_nodesize,
      importance    = "permutation",
      num.threads   = 1
    )
    data.frame(
      ASV        = names(m$variable.importance),
      Importance = as.numeric(m$variable.importance),
      Management = name_map[cls],
      stringsAsFactors = FALSE
    )
  })
  bind_rows(imp_list)
}

# DIAGNOSTIC PLOT: 3×3 per-class importance grid (from Doc 5)
# One panel per Management × Year for a given kingdom × compartment

make_kingdom_grid <- function(kingdom_set, comp_suf, top_n = 15) {
  
  comp_full  <- compartments[comp_suffix == comp_suf]
  comp_label <- compartment_labels[comp_full]
  grid_title <- paste0(kingdom_labels[kingdom_set], " — ", comp_label)
  plot_list  <- list()
  
  for (mg in mgmt_display) {
    for (yr in years) {
      
      plot_id  <- paste0(mg, "_", yr)
      model_id <- paste0(yr, "_", comp_suf, "_", kingdom_set)
      
      if (!model_id %in% names(all_results)) {
        plot_list[[plot_id]] <- ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                   label = "No result", size = 3, color = "grey50") +
          theme_void()
        next
      }
      
      result   <- all_results[[model_id]]
      ps_objs  <- get_ps_objects(yr, comp_suf)
      ps_b     <- ps_objs$ps_bac
      ps_f     <- ps_objs$ps_fun
      df_final <- all_df_final[[model_id]]
      
      imp_df <- get_perclass_importance(
        df            = df_final,
        seed          = 42,
        best_mtry     = result$model_cv$bestTune$mtry,
        best_nodesize = result$model_cv$bestTune$min.node.size,
        best_ntrees   = result$best_ntrees
      )
      
      imp_mg <- imp_df %>%
        filter(Management == mg) %>%
        arrange(desc(Importance)) %>%
        slice_head(n = top_n)
      
      tax_labels <- rbind(
        make_tax_labels(if (kingdom_set %in% c("bac", "combined")) ps_b else NULL, "bac"),
        make_tax_labels(if (kingdom_set %in% c("fun", "combined")) ps_f else NULL, "fun")
      )
      
      if (!is.null(tax_labels)) {
        imp_mg <- left_join(imp_mg, tax_labels, by = "ASV")
        imp_mg$label <- ifelse(is.na(imp_mg$label), imp_mg$ASV, imp_mg$label)
      } else {
        imp_mg$label <- imp_mg$ASV
      }
      
      imp_mg <- imp_mg %>% mutate(label = factor(label, levels = rev(label)))
      acc    <- round(result$cv_accuracy * 100, 1)
      
      p <- ggplot(imp_mg, aes(x = Importance, y = label)) +
        geom_segment(aes(x = 0, xend = Importance, y = label, yend = label),
                     color = management_colors[mg], linewidth = 0.6, alpha = 0.7) +
        geom_point(color = management_colors[mg], size = 3) +
        labs(
          title    = paste0(mg, " — ", year_labels[yr]),
          subtitle = paste0("CV Accuracy: ", acc, "%"),
          x        = "Permutation Importance",
          y        = NULL
        ) +
        theme_classic(base_size = 9) +
        theme(
          axis.text.y   = element_text(size = 6),
          plot.title    = element_text(face = "bold", size = 9,
                                       color = management_colors[mg]),
          plot.subtitle = element_text(size = 8, color = "grey40")
        )
      
      plot_list[[plot_id]] <- p
    }
  }
  
  inner_grid <- plot_grid(
    plot_list[["Conventional_2018"]], plot_list[["Conventional_2019"]], plot_list[["Conventional_2020"]],
    plot_list[["No-Till_2018"]],      plot_list[["No-Till_2019"]],      plot_list[["No-Till_2020"]],
    plot_list[["Organic_2018"]],      plot_list[["Organic_2019"]],      plot_list[["Organic_2020"]],
    ncol = 3, nrow = 3,
    labels = c("A","B","C","D","E","F","G","H","I"), label_size = 10
  )
  
  title_grob <- ggdraw() +
    draw_label(grid_title, fontface = "bold", size = 13, x = 0.5, hjust = 0.5)
  
  plot_grid(title_grob, inner_grid, ncol = 1, rel_heights = c(0.05, 1))
}


# DIAGNOSTIC PLOT: confusion matrix (from Doc 5)

plot_confusion <- function(conf_df, title = "", accuracy = NULL) {
  conf_df$Actual    <- factor(conf_df$Actual,    levels = mgmt_display)
  conf_df$Predicted <- factor(conf_df$Predicted, levels = mgmt_display)
  acc_label <- if (!is.null(accuracy)) paste0("\nCV Accuracy: ", round(accuracy * 100, 1), "%") else ""
  
  ggplot(conf_df, aes(x = Predicted, y = Actual, fill = prop)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = paste0(round(prop * 100, 1), "%\n(n=", n, ")")),
              size = 3.5, color = "white", fontface = "bold") +
    scale_fill_gradient(low = "#f0f0f0", high = "#1a3a5c",
                        limits = c(0, 1), name = "Proportion\nCorrect") +
    scale_x_discrete(position = "top") +
    labs(title = paste0(title, acc_label),
         x = "Predicted Management", y = "Actual Management") +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x     = element_text(angle = 30, hjust = 0, size = 10),
      axis.text.y     = element_text(size = 10),
      plot.title      = element_text(face = "bold", size = 11),
      legend.position = "right"
    )
}

# MAIN LOOP: 27 models

all_results  <- list()
all_rfe      <- list()
all_df_final <- list()

for (yr in years) {
  for (ci in seq_along(compartments)) {
    
    comp    <- compartments[ci]
    cs      <- comp_suffix[ci]
    ps_objs <- get_ps_objects(yr, cs)
    
    if (is.null(ps_objs)) {
      cat("  Skipping", yr, comp, "— objects not found\n")
      next
    }
    
    ps_b <- ps_objs$ps_bac
    ps_f <- ps_objs$ps_fun
    
    for (kg in kingdoms) {
      
      model_id <- paste0(yr, "_", cs, "_", kg)
      cat("\n=== RF:", yr, "|", comp, "|", kg, "===\n")
      
      df <- tryCatch(
        build_feature_df(
          ps_bac  = if (kg %in% c("bac", "combined")) ps_b else NULL,
          ps_fun  = if (kg %in% c("fun", "combined")) ps_f else NULL,
          kingdom = kg
        ),
        error = function(e) {
          cat("  ERROR building features:", conditionMessage(e), "\n"); NULL
        }
      )
      
      if (is.null(df) || nrow(df) < 10) {
        cat("  Skipping — insufficient samples\n"); next
      }
      
      cat("  Initial: Samples =", nrow(df), "| Features =", ncol(df) - 1, "\n")
      print(table(df$Management))
      
      cat("  Running RFE...\n")
      rfe_out <- tryCatch(
        run_rfe(df, seed = 42),
        error = function(e) {
          cat("  ERROR in RFE:", conditionMessage(e), "\n"); NULL
        }
      )
      
      if (is.null(rfe_out)) next
      all_rfe[[model_id]]      <- rfe_out$rfe_result
      all_df_final[[model_id]] <- rfe_out$df_rfe
      cat("  Post-RFE features:", ncol(rfe_out$df_rfe) - 1, "\n")
      
      cat("  Running full RF with hyperparameter tuning...\n")
      result <- tryCatch(
        run_rf_cv(rfe_out$df_rfe, seed = 42),
        error = function(e) {
          cat("  ERROR in RF:", conditionMessage(e), "\n"); NULL
        }
      )
      
      if (is.null(result)) next
      cat("  Final CV Accuracy:", round(result$cv_accuracy * 100, 1), "%\n")
      all_results[[model_id]] <- result
      
      # Incremental save — so crashes don't lose progress
      saveRDS(result,         file = paste0("rf_result_",   model_id, ".rds"))
      saveRDS(rfe_out,        file = paste0("rf_rfe_",      model_id, ".rds"))
      saveRDS(rfe_out$df_rfe, file = paste0("rf_df_final_", model_id, ".rds"))
    }
  }
}

# Save full collected objects
saveRDS(all_results,  "rf_all_results.rds")
saveRDS(all_rfe,      "rf_all_rfe.rds")
saveRDS(all_df_final, "rf_all_df_final.rds")
cat("\n✓ All models saved\n")

# SUMMARY TABLE

summary_df <- data.frame(
  Model              = names(all_results),
  Year               = sapply(strsplit(names(all_results), "_"), `[`, 1),
  Compartment        = sapply(strsplit(names(all_results), "_"), function(x)
    compartment_labels[compartments[comp_suffix == x[2]]]),
  Kingdom            = sapply(strsplit(names(all_results), "_"), `[`, 3),
  N_features_initial = 500,
  N_features_rfe     = sapply(names(all_results), function(mid)
    ncol(all_df_final[[mid]]) - 1),
  Accuracy           = sapply(all_results, function(r) round(r$cv_accuracy * 100, 1)),
  Best_ntrees        = sapply(all_results, function(r) r$best_ntrees),
  Best_mtry          = sapply(all_results, function(r) r$model_cv$bestTune$mtry),
  Best_nodesize      = sapply(all_results, function(r) r$model_cv$bestTune$min.node.size)
)
rownames(summary_df) <- NULL
print(summary_df)
saveRDS(summary_df, "rf_summary_table.rds")


# FIGURE: PANEL A — Accuracy heatmap

heatmap_df <- summary_df %>%
  mutate(
    Kingdom     = factor(Kingdom,
                         levels = c("bac", "fun", "combined"),
                         labels = c("Bacteria", "Fungi", "Combined")),
    Year_label  = factor(Year,
                         levels = c("2018", "2019", "2020"),
                         labels = c("2018\nSoybean", "2019\nWheat", "2020\nMaize")),
    Compartment = factor(Compartment,
                         levels = c("All Compartments", "Aboveground", "Belowground"))
  )

p_heatmap <- ggplot(heatmap_df, aes(x = Year_label, y = Kingdom, fill = Accuracy)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(
    label = paste0(Accuracy, "%"),
    color = Accuracy > 75
  ),  size = 4, fontface = "bold") +
  scale_color_manual(values = c("black", "white"), guide = "none") +
  scale_fill_gradient(
    low    = "#f5f1e8",
    high   = "#1a3a5c",
    limits = c(0, 100),
    name   = "CV Accuracy (%)"
  ) +
  facet_wrap(~ Compartment, ncol = 3) +
  labs(
    title    = "Random Forest Classification Accuracy",
    subtitle = "Management prediction from microbiome composition — 5-fold repeated CV",
    x = NULL, y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text       = element_text(face = "bold", size = 12, color = "black"),
    strip.background = element_rect(fill = "grey92", color = NA),
    axis.text        = element_text(face = "bold", color = "black", size = 11),
    axis.text.x      = element_text(lineheight = 0.9),
    plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle    = element_text(size = 10, hjust = 0.5, color = "grey40"),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 10),
    panel.spacing    = unit(1, "lines")
  )


# FIGURE: PANEL B — Cross-year consistent features
# Combined pooled model only

consistent_df <- lapply(years, function(yr) {
  model_id <- paste0(yr, "_all_combined")
  if (!model_id %in% names(all_results)) return(NULL)
  
  result   <- all_results[[model_id]]
  df_final <- all_df_final[[model_id]]
  ps_objs  <- get_ps_objects(yr, "all")
  
  imp_df <- get_perclass_importance(
    df            = df_final,
    seed          = 42,
    best_mtry     = result$model_cv$bestTune$mtry,
    best_nodesize = result$model_cv$bestTune$min.node.size,
    best_ntrees   = result$best_ntrees
  )
  
  tax_labels <- rbind(
    make_tax_labels(ps_objs$ps_bac, "bac"),
    make_tax_labels(ps_objs$ps_fun, "fun")
  )
  
  imp_df <- left_join(imp_df, tax_labels, by = "ASV")
  imp_df$label <- ifelse(is.na(imp_df$label), imp_df$ASV, imp_df$label)
  imp_df$Year  <- yr
  imp_df
}) %>%
  bind_rows()

# ASVs present in 2+ years with positive importance
consistent_asvs <- consistent_df %>%
  filter(Importance > 0) %>%
  group_by(ASV) %>%
  summarise(n_years = n_distinct(Year), .groups = "drop") %>%
  filter(n_years >= 2) %>%
  pull(ASV)

plot_consistent <- consistent_df %>%
  filter(ASV %in% consistent_asvs) %>%
  group_by(ASV, label, Management) %>%
  summarise(
    mean_importance = mean(Importance),
    n_years         = n_distinct(Year),
    .groups         = "drop"
  ) %>%
  group_by(Management) %>%
  slice_max(mean_importance, n = 15) %>%
  ungroup() %>%
  mutate(
    label      = factor(label, levels = unique(label[order(mean_importance)])),
    Management = factor(Management, levels = mgmt_display)
  )

p_consistent <- ggplot(plot_consistent,
                       aes(x = mean_importance, y = label, color = Management)) +
  geom_segment(aes(x = 0, xend = mean_importance, y = label, yend = label),
               linewidth = 0.7, alpha = 0.6) +
  geom_point(aes(size = n_years), alpha = 0.9) +
  scale_color_manual(values = management_colors, name = "Management") +
  scale_size_continuous(
    range  = c(2, 5),
    breaks = c(2, 3),
    labels = c("2 years", "3 years"),
    name   = "Consistency"
  ) +
  facet_wrap(~ Management, ncol = 3, scales = "free_x") +
  labs(
    title    = "Cross-Year Consistent Features — Combined Model (Pooled)",
    subtitle = "ASVs with positive importance in ≥2 years; point size = years present",
    x        = "Mean Permutation Importance",
    y        = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.text       = element_text(face = "bold", size = 11, color = "black"),
    strip.background = element_rect(fill = "grey92", color = NA),
    axis.text.y      = element_text(size = 7,  face = "bold", color = "black"),
    axis.text.x      = element_text(size = 9,  face = "bold", color = "black"),
    axis.title.x     = element_text(size = 10, face = "bold", color = "black"),
    plot.title       = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle    = element_text(size = 9,  hjust = 0.5, color = "grey40"),
    legend.position  = "right",
    legend.title     = element_text(face = "bold"),
    panel.spacing    = unit(1, "lines")
  )

# COMBINE AND PRINT
# Export manually from RStudio at 7000 x 5000

fig_rf_summary <- plot_grid(
  p_heatmap,
  p_consistent,
  ncol        = 1,
  rel_heights = c(1, 1.6),
  labels      = c("A", "B"),
  label_size  = 16,
  label_fontface = "bold"
)
print(fig_rf_summary)

# OOB vs CV ACCURACY CHECK

oob_df <- data.frame(
  Model        = names(all_results),
  CV_accuracy  = sapply(all_results, function(r) round(max(r$model_cv$results$Accuracy) * 100, 1)),
  OOB_accuracy = sapply(all_results, function(r) round((1 - r$model_ranger$prediction.error) * 100, 1))
)
oob_df$Delta <- oob_df$CV_accuracy - oob_df$OOB_accuracy
print(oob_df)

# OPTIONAL DIAGNOSTICS (from Doc 5)
# Confusion matrices + 3×3 per-class importance grids.

if (RENDER_DIAGNOSTICS) {
  
  # --- 9 importance grids: each kingdom × compartment ---
  for (kg in kingdoms) {
    for (cs in comp_suffix) {
      cat("\n=== Grid:", kingdom_labels[kg], "|",
          compartment_labels[compartments[comp_suffix == cs]], "===\n")
      g <- tryCatch(
        make_kingdom_grid(kg, cs),
        error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
      )
      if (!is.null(g)) print(g)
    }
  }
  
  # --- Confusion matrices, grouped by compartment ---
  for (cs in comp_suffix) {
    comp_label <- compartment_labels[compartments[comp_suffix == cs]]
    conf_plots <- list()
    for (kg in kingdoms) {
      for (yr in years) {
        model_id <- paste0(yr, "_", cs, "_", kg)
        if (!model_id %in% names(all_results)) next
        conf_plots[[model_id]] <- plot_confusion(
          conf_df  = all_results[[model_id]]$conf_df,
          title    = paste0(kingdom_labels[kg], "\n", year_labels[yr], " — ", comp_label),
          accuracy = all_results[[model_id]]$cv_accuracy
        )
      }
    }
    if (length(conf_plots) > 0) {
      cat("\n=== Confusion Matrices —", comp_label, "===\n")
      print(plot_grid(plotlist = conf_plots, ncol = 3, labels = "AUTO"))
    }
  }
  
  # --- Cross-year consistency tables (console), all kg × cs ---
  compare_importance <- function(kingdom_set = "bac", comp_suf = "all", top_n = 15) {
    ids <- paste0(years, "_", comp_suf, "_", kingdom_set)
    ids <- ids[ids %in% names(all_results)]
    lapply(ids, function(mid) {
      all_results[[mid]]$imp_df %>%
        arrange(desc(Importance)) %>%
        slice_head(n = top_n) %>%
        mutate(Year = sub("_.*", "", mid))
    }) %>%
      bind_rows() %>%
      group_by(ASV) %>%
      mutate(n_years = n_distinct(Year)) %>%
      ungroup()
  }
  
  for (kg in kingdoms) {
    for (cs in comp_suffix) {
      cat("\n=== Cross-year consistency:", kg, "|",
          compartment_labels[compartments[comp_suffix == cs]], "===\n")
      print(compare_importance(kg, cs) %>% filter(n_years > 1) %>% arrange(ASV, Year))
    }
  }
}

message("\n========== ALL COMPLETE ==========")
message("fig_rf_summary — export manually at 7000 x 5000")
message("rf_all_results.rds / rf_all_rfe.rds / rf_all_df_final.rds — saved")
if (RENDER_DIAGNOSTICS) message("Diagnostics rendered: grids + confusion matrices + consistency tables")

# RF VISUALIZATION — SEQUENCE-KEYED, self-contained
# Cross-year panels (C/D, persistence, Table RF2) keyed by DNA sequence.
# Main E = kingdom contribution; Supplemental D = persistence.

setwd("~/Desktop/Paper/NEW_R_Code/R_Code_Files/paper-final-versions/final_figures_for_realsies/rf/rf_for_man")
all_results  <- readRDS("rf_all_results.rds")
all_rfe      <- readRDS("rf_all_rfe.rds")
all_df_final <- readRDS("rf_all_df_final.rds")
summary_df   <- readRDS("rf_summary_table.rds")

get_ps_objects <- function(yr, cs) {
  if (cs == "all") { ps_b <- paste0("psr_", yr, "_bac"); ps_f <- paste0("psr_", yr) }
  else             { ps_b <- paste0("psr_", yr, "_", cs, "_bac"); ps_f <- paste0("psr_", yr, "_", cs) }
  if (!exists(ps_b) || !exists(ps_f)) return(NULL)
  list(ps_bac = get(ps_b), ps_fun = get(ps_f))
}

get_perclass_importance <- function(df, seed = 42, best_mtry = NULL, best_nodesize = 5, best_ntrees = 1000) {
  if (is.null(best_mtry)) best_mtry <- floor(sqrt(ncol(df) - 1))
  bind_rows(lapply(levels(df$Management), function(cls) {
    df_bin <- df
    df_bin$Management <- factor(ifelse(df$Management == cls, "target", "other"), levels = c("target","other"))
    set.seed(seed)
    m <- ranger(Management ~ ., data = df_bin, num.trees = best_ntrees, mtry = best_mtry,
                min.node.size = best_nodesize, importance = "permutation", num.threads = 1)
    data.frame(ASV = names(m$variable.importance), Importance = as.numeric(m$variable.importance),
               Management = name_map[cls], stringsAsFactors = FALSE)
  }))
}

make_tax_labels <- function(ps, prefix) {
  if (is.null(ps)) return(NULL)
  tax <- as.data.frame(tax_table(ps)); tax$ASV <- paste0(prefix, "_", rownames(tax))
  ktag <- ifelse(prefix == "bac", "[B] ", ifelse(prefix == "fun", "[F] ", ""))
  clean <- function(x) {
    x <- sub("_[0-9.]+$", "", x); x <- gsub("_", " ", x); x[is.na(x) | x == ""] <- NA
    bad <- grepl("Triticum|Bromus|Zea mays|eukaryote|archaeote|Mitochondria|Chloroplast|uncultured|unidentified|metagenome",
                 x, ignore.case = TRUE)
    x[bad] <- NA; x
  }
  has_rdp <- "Genus_RDP" %in% colnames(tax)
  get_rank <- function(rank) {
    plain <- if (rank %in% colnames(tax)) clean(tax[[rank]]) else rep(NA, nrow(tax))
    if (!has_rdp) return(plain)
    rdp   <- if (paste0(rank,"_RDP")   %in% colnames(tax)) clean(tax[[paste0(rank,"_RDP")]])   else rep(NA, nrow(tax))
    blast <- if (paste0(rank,"_BLAST") %in% colnames(tax)) clean(tax[[paste0(rank,"_BLAST")]]) else rep(NA, nrow(tax))
    coalesce(plain, rdp, blast)
  }
  best <- coalesce(get_rank("Genus"), get_rank("Family"), get_rank("Order"), get_rank("Class"), get_rank("Phylum"))
  tax$label <- ifelse(!is.na(best), paste0(ktag, best, " (", rownames(tax), ")"), paste0(ktag, rownames(tax)))
  tax[, c("ASV","label")]
}

get_imp <- function(model_id) {
  result <- all_results[[model_id]]; df_final <- all_df_final[[model_id]]
  parts <- strsplit(model_id, "_")[[1]]; yr <- parts[1]; cs <- parts[2]; kg <- parts[3]
  ps_objs <- get_ps_objects(yr, cs); if (is.null(ps_objs)) return(NULL)
  imp <- get_perclass_importance(df_final, seed = 42,
                                 best_mtry = result$model_cv$bestTune$mtry,
                                 best_nodesize = result$model_cv$bestTune$min.node.size,
                                 best_ntrees = result$best_ntrees)
  tax_labels <- if (kg == "bac") make_tax_labels(ps_objs$ps_bac, "bac")
  else if (kg == "fun") make_tax_labels(ps_objs$ps_fun, "fun")
  else rbind(make_tax_labels(ps_objs$ps_bac, "bac"), make_tax_labels(ps_objs$ps_fun, "fun"))
  imp <- left_join(imp, tax_labels, by = "ASV")
  imp$label <- ifelse(is.na(imp$label), imp$ASV, imp$label)
  imp$Year <- yr; imp$Compartment <- compartment_labels[compartments[comp_suffix == cs]]
  imp$Kingdom <- kingdom_labels[kg]; imp$model_id <- model_id
  imp
}

cat("Building master importance table...\n")
imp_all <- bind_rows(lapply(names(all_results), function(mid) {
  cat("  ", mid, "\n"); tryCatch(get_imp(mid), error = function(e) { cat("  ERROR:", e$message, "\n"); NULL })
}))

# ---- attach DNA sequence (the cross-year key) ------------------------
attach_sequence <- function(imp_df) {
  one <- function(asv, yr, cs) {
    ps_objs <- get_ps_objects(yr, cs); if (is.null(ps_objs)) return(NA_character_)
    if (startsWith(asv, "bac_")) { ps <- ps_objs$ps_bac; id <- sub("^bac_", "", asv) }
    else if (startsWith(asv, "fun_")) { ps <- ps_objs$ps_fun; id <- sub("^fun_", "", asv) }
    else return(NA_character_)
    rs <- tryCatch(refseq(ps), error = function(e) NULL)
    if (is.null(rs) || !id %in% names(rs)) return(NA_character_)
    as.character(rs[id])
  }
  cs_of <- setNames(comp_suffix, compartment_labels[compartments])
  imp_df$Sequence <- mapply(one, imp_df$ASV, imp_df$Year, cs_of[as.character(imp_df$Compartment)], USE.NAMES = FALSE)
  imp_df
}
imp_all <- attach_sequence(imp_all)
cat("Done:", nrow(imp_all), "rows |", sum(!is.na(imp_all$Sequence)), "with sequence\n")

rf_theme <- theme_classic(base_size = 10) +
  theme(strip.text = element_text(face="bold", size=9, color="black"),
        strip.background = element_rect(fill="grey92", color=NA),
        axis.text = element_text(face="bold", color="black"),
        axis.text.y = element_text(size=7), axis.text.x = element_text(size=8),
        plot.title = element_text(face="bold", size=13, hjust=0.5),
        plot.subtitle = element_text(size=9, hjust=0.5, color="grey40"),
        legend.position = "right", panel.spacing = unit(0.5, "lines"))

heatmap_df <- summary_df %>% mutate(
  Kingdom = factor(Kingdom, levels=c("bac","fun","combined"), labels=c("Bacteria","Fungi","Combined")),
  Year_label = factor(Year, levels=c("2018","2019","2020"), labels=c("2018\nSoybean","2019\nWheat","2020\nMaize")),
  Compartment = factor(Compartment, levels=c("All Compartments","Aboveground","Belowground")))

oob_df <- data.frame(
  Model = names(all_results),
  CV_accuracy  = sapply(all_results, function(r) round(max(r$model_cv$results$Accuracy)*100,1)),
  OOB_accuracy = sapply(all_results, function(r) round((1-r$model_ranger$prediction.error)*100,1))) %>%
  mutate(Year = sapply(strsplit(Model,"_"), `[`, 1), cs = sapply(strsplit(Model,"_"), `[`, 2),
         Compartment = factor(compartment_labels[compartments[match(cs, comp_suffix)]],
                              levels=c("All Compartments","Aboveground","Belowground")),
         Kingdom = factor(sapply(strsplit(Model,"_"), `[`, 3), levels=c("bac","fun","combined"),
                          labels=c("Bacteria","Fungi","Combined")))

acc_gain <- heatmap_df %>% select(Year, Compartment, Kingdom, Accuracy) %>%
  pivot_wider(names_from = Kingdom, values_from = Accuracy) %>%
  mutate(Best_single = pmax(Bacteria, Fungi, na.rm = TRUE), Gain_combined = Combined - Best_single) %>%
  filter(!is.na(Gain_combined))
acc_ab <- heatmap_df %>% filter(Compartment %in% c("Aboveground","Belowground")) %>%
  select(Year, Kingdom, Compartment, Accuracy) %>%
  pivot_wider(names_from = Compartment, values_from = Accuracy)

# ---- PANEL A: accuracy heatmap ---------------------------------------
p_accuracy <- ggplot(heatmap_df, aes(x = Year_label, y = Kingdom, fill = Accuracy)) +
  geom_tile(color="white", linewidth=1) +
  geom_text(aes(label=paste0(Accuracy,"%"), color=Accuracy>75), size=4, fontface="bold") +
  scale_color_manual(values=c("black","white"), guide="none") +
  scale_fill_gradient(low="#f5f1e8", high="#4c4184", limits=c(0,100), name="CV Accuracy (%)") +
  facet_wrap(~ Compartment, ncol=3) +
  labs(title="Random Forest Classification Accuracy",
       subtitle="Management prediction from microbiome composition, 5-fold repeated CV", x=NULL, y=NULL) +
  theme_classic(base_size=12) +
  theme(strip.text=element_text(face="bold", size=12, color="black"),
        strip.background=element_rect(fill="grey92", color=NA),
        axis.text=element_text(face="bold", color="black", size=11),
        axis.text.x=element_text(lineheight=0.9),
        plot.title=element_text(face="bold", size=14, hjust=0.5),
        plot.subtitle=element_text(size=10, hjust=0.5, color="grey40"),
        legend.position="right", panel.spacing=unit(1,"lines"))

# ---- PANEL B: per-class recall ---------------------------------------
perclass_acc <- bind_rows(lapply(names(all_results), function(mid) {
  r <- all_results[[mid]]; parts <- strsplit(mid, "_")[[1]]
  r$conf_df %>% group_by(Actual) %>% mutate(Recall = n/sum(n)) %>%
    filter(Actual == Predicted) %>% select(Actual, Recall) %>%
    mutate(Year = parts[1], Compartment = compartment_labels[compartments[comp_suffix == parts[2]]],
           Kingdom = kingdom_labels[parts[3]], model_id = mid)
}))
p_perclass <- perclass_acc %>%
  mutate(Kingdom = factor(Kingdom, levels=c("Bacterial Microbiome","Fungal Microbiome","Combined Microbiome")),
         Compartment = factor(Compartment, levels=c("All Compartments","Aboveground","Belowground")),
         Year_label = paste0(Year, "\n", case_when(Year=="2018"~"Soy", Year=="2019"~"Wheat", Year=="2020"~"Maize"))) %>%
  ggplot(aes(x=Year_label, y=Recall, fill=Actual)) +
  geom_col(position="dodge", width=0.7) +
  geom_hline(yintercept=0.33, linetype="dashed", color="grey50", linewidth=0.5) +
  scale_fill_manual(values=management_colors, name="Management") +
  scale_y_continuous(limits=c(0,1), labels=scales::percent) +
  facet_grid(Compartment ~ Kingdom) +
  labs(title="Per-Class Recall by Model", subtitle="Dashed line = random chance (0.33)", x=NULL, y="Recall") +
  rf_theme + theme(axis.text.x = element_text(angle=45, hjust=1, size=7))

# ---- PANELS C + D: cross-year consistent features (BY SEQUENCE) ------
make_consistent <- function(comp_filter) {
  imp_all %>%
    filter(Kingdom != "Combined Microbiome", Compartment == comp_filter, Importance > 0, !is.na(Sequence)) %>%
    group_by(Sequence, Management, Kingdom) %>%
    summarise(mean_importance = mean(Importance, na.rm = TRUE), n_years = n_distinct(Year),
              taxon = sub(" \\([^)]*\\)$", "", label[which.max(Importance)]), .groups = "drop") %>%
    filter(n_years >= 2) %>%
    group_by(Management, Kingdom) %>% slice_max(mean_importance, n = 10) %>% ungroup() %>%
    mutate(taxon = factor(taxon, levels = unique(taxon[order(mean_importance)])),
           Kingdom = factor(Kingdom, levels = c("Bacterial Microbiome","Fungal Microbiome")),
           Years = factor(n_years, levels = c(2,3), labels = c("2 years","3 years"))) %>%
    ggplot(aes(x = mean_importance, y = taxon, color = Management)) +
    geom_segment(aes(x=0, xend=mean_importance, y=taxon, yend=taxon), linewidth=0.7, alpha=0.6) +
    geom_point(aes(size = Years), alpha = 0.9) +
    scale_color_manual(values = management_colors) +
    scale_size_manual(values = c("2 years"=2.5,"3 years"=5), drop = FALSE, name = "Consistency") +
    facet_grid(Kingdom ~ Management, scales = "free") +
    labs(title = paste0("Cross-Year Consistent Features (by sequence), ", comp_filter),
         subtitle = "Top 10 organisms with positive importance in >=2 years",
         x = "Mean Permutation Importance", y = NULL) + rf_theme
}
p_consistent_above <- make_consistent("Aboveground")
p_consistent_below <- make_consistent("Belowground")

# ---- PERSISTENCE (BY SEQUENCE) -> SUPPLEMENTAL D ---------------------
persistence_df <- imp_all %>%
  filter(Kingdom != "Combined Microbiome", Importance > 0, !is.na(Sequence)) %>%
  group_by(Sequence, Management, Kingdom, Compartment) %>%
  summarise(n_years = n_distinct(Year), .groups = "drop") %>%
  mutate(Persistence = factor(c("1 year","2 years","3 years")[n_years], levels=c("1 year","2 years","3 years")),
         Kingdom = factor(Kingdom, levels=c("Bacterial Microbiome","Fungal Microbiome")),
         Compartment = factor(Compartment, levels=c("All Compartments","Aboveground","Belowground")))
p_persistence <- persistence_df %>%
  filter(as.character(Compartment) != "All Compartments") %>%
  group_by(Kingdom, Compartment, Management, Persistence) %>% summarise(n = n(), .groups = "drop") %>%
  ggplot(aes(x = Management, y = n, fill = Persistence)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values=c("1 year"="#f5f1e8","2 years"="#7A6E3C","3 years"="tomato4"), name="Years as\nclassifier") +
  facet_grid(Compartment ~ Kingdom) +
  labs(title="Classifier Taxa Persistence Across Crop Rotation (by sequence)",
       subtitle="Number of organisms with positive importance in 1, 2, or 3 years", x=NULL, y="Number of taxa") +
  rf_theme + theme(axis.text.x = element_text(angle=45, hjust=1))

# ---- KINGDOM CONTRIBUTION -> MAIN E ----------------------------------
p_kingdom_contrib <- imp_all %>%
  filter(Kingdom == "Combined Microbiome", Importance > 0) %>%
  mutate(kingdom_type = ifelse(grepl("^bac_", ASV), "Bacteria", "Fungi")) %>%
  group_by(Year, Compartment, Management, kingdom_type) %>% summarise(total_importance = sum(Importance), .groups="drop") %>%
  group_by(Year, Compartment, Management) %>% mutate(pct = 100*total_importance/sum(total_importance)) %>% ungroup() %>%
  mutate(Compartment = factor(Compartment, levels=c("All Compartments","Aboveground","Belowground"))) %>%
  ggplot(aes(x=Management, y=pct, fill=kingdom_type)) +
  geom_col(position="stack", width=0.7) +
  geom_hline(yintercept=50, linetype="dashed", color="grey40", linewidth=0.5) +
  scale_fill_manual(values=c("Bacteria"="#83331D","Fungi"="#22569F"), name="Kingdom") +
  facet_grid(Compartment ~ Year) +
  labs(title="Kingdom Contribution to Combined Model Importance",
       subtitle="% of total permutation importance from each kingdom", x=NULL, y="% Total Importance") +
  rf_theme + theme(axis.text.x = element_text(angle=45, hjust=1))

# ---- SUPPLEMENTAL A/B/C ----------------------------------------------
p_oob <- ggplot(oob_df, aes(x=CV_accuracy, y=OOB_accuracy, color=Kingdom)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey60") +
  geom_point(aes(shape=Compartment), size=3, alpha=0.85) +
  scale_color_manual(values=kingdom_colors) +
  scale_shape_manual(values=c("All Compartments"=16,"Aboveground"=17,"Belowground"=15)) +
  labs(title="OOB vs CV Accuracy", subtitle="Points near dashed line = no overfitting",
       x="CV Accuracy (%)", y="OOB Accuracy (%)") + rf_theme

p_ab_accuracy <- ggplot(acc_ab, aes(x=Aboveground, y=Belowground, color=Kingdom)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey60") +
  geom_point(aes(shape=Year), size=4, alpha=0.9) +
  scale_color_manual(values=kingdom_colors) +
  scale_shape_manual(values=c("2018"=16,"2019"=17,"2020"=15)) +
  labs(title="Aboveground vs Belowground Classification Accuracy",
       subtitle="Points above dashed line = belowground more classifiable",
       x="Aboveground CV Accuracy (%)", y="Belowground CV Accuracy (%)") + rf_theme

p_gain <- ggplot(acc_gain, aes(x=Year, y=Gain_combined, fill=Compartment)) +
  geom_col(position="dodge", width=0.7) +
  geom_hline(yintercept=0, linetype="dashed", color="grey40") +
  scale_fill_manual(values=compartment_colors) +
  labs(title="Combined Model Accuracy Gain over Best Single-Kingdom",
       subtitle="Positive = combined outperforms best single-kingdom model", x=NULL, y="Accuracy Gain (%)") + rf_theme

# ---- TABLES ----------------------------------------------------------
table_rf1 <- summary_df %>%
  mutate(OOB_Accuracy = sapply(Model, function(mid) round((1 - all_results[[mid]]$model_ranger$prediction.error)*100,1)),
         Overfit_Delta = OOB_Accuracy - Accuracy,
         Kingdom = factor(Kingdom, levels=c("bac","fun","combined"), labels=c("Bacteria","Fungi","Combined")),
         Crop = case_when(Year=="2018"~"Soybean", Year=="2019"~"Wheat", Year=="2020"~"Maize")) %>%
  select(Year, Crop, Compartment, Kingdom, CV_Accuracy = Accuracy, OOB_Accuracy, Overfit_Delta,
         N_features_RFE = N_features_rfe, Best_ntrees, Best_mtry, Best_nodesize) %>%
  arrange(Year, Compartment, Kingdom)
write.csv(table_rf1, "Supplementary_Table_RF1_ModelPerformance.csv", row.names = FALSE)
cat("Saved Table RF-1:", nrow(table_rf1), "rows\n")

table_rf2 <- imp_all %>%
  filter(Kingdom != "Combined Microbiome", Importance > 0, !is.na(Sequence)) %>%
  group_by(Sequence, Management, Kingdom, Compartment) %>%
  summarise(Taxonomy = sub(" \\([^)]*\\)$", "", label[which.max(Importance)]),
            Representative_ASV = ASV[which.max(Importance)],
            Mean_Importance = round(mean(Importance, na.rm = TRUE), 6),
            Max_Importance  = round(max(Importance, na.rm = TRUE), 6),
            N_years = n_distinct(Year),
            Years_present = paste(sort(unique(Year)), collapse = ", "), .groups = "drop") %>%
  mutate(Persistence = case_when(N_years==3 ~ "3 years - full rotation",
                                 N_years==2 ~ "2 years - partial", N_years==1 ~ "1 year - specific"),
         Kingdom = factor(Kingdom, levels=c("Bacterial Microbiome","Fungal Microbiome")),
         Compartment = factor(Compartment, levels=c("All Compartments","Aboveground","Belowground"))) %>%
  arrange(Kingdom, Compartment, Management, desc(N_years), desc(Mean_Importance)) %>%
  select(Kingdom, Compartment, Management, Taxonomy, Representative_ASV,
         Persistence, N_years, Years_present, Mean_Importance, Max_Importance)
write.csv(table_rf2, "Supplementary_Table_RF2_ClassifierTaxa.csv", row.names = FALSE)
cat("Saved Table RF-2 (by sequence):", nrow(table_rf2), "rows\n")

# ---- ASSEMBLE (E = kingdom contrib; sup D = persistence) -------------
row1_fig4 <- plot_grid(p_accuracy, p_perclass, ncol=2, rel_widths=c(1,1.2),
                       labels=c("A","B"), label_size=14, label_fontface="bold")
row2_fig4 <- plot_grid(p_consistent_above, p_consistent_below, ncol=2, rel_widths=c(1,1),
                       labels=c("C","D"), label_size=14, label_fontface="bold")
row3_fig4 <- plot_grid(p_kingdom_contrib, ncol=1, labels=c("E"), label_size=14, label_fontface="bold")
fig4 <- plot_grid(row1_fig4, row2_fig4, row3_fig4, ncol=1, rel_heights=c(1,1.4,0.9))

row1_sup <- plot_grid(p_oob, p_ab_accuracy, ncol=2, rel_widths=c(1,1),
                      labels=c("A","B"), label_size=14, label_fontface="bold")
row2_sup <- plot_grid(p_gain, p_persistence, ncol=2, rel_widths=c(0.8,1.2),
                      labels=c("C","D"), label_size=14, label_fontface="bold")
fig4_sup <- plot_grid(row1_sup, row2_sup, ncol=1, rel_heights=c(1,1.2))

message("===== FIG 4 (sequence-keyed C/D, kingdom-contrib E) ====="); print(fig4)
message("===== SUPPLEMENTAL RF (persistence at D) ====="); print(fig4_sup)
message("===== COMPLETE — cross-year analysis is SEQUENCE-based =====")

# =============================================================================
# =============================================================================

#  SECTION 14
#  Differential abundance tables and volcano plots (unused in manuscript but may as well keep code)
#
#  Supplemental Tables 7 and 8.
#  Writes the two DESeq2_*_SignificantASVs.csv files, read by section 15.

# ====== ENSURE per-year UNRAREFIED objects exist ======

if (!exists("ps_unrare_2018_bac") && exists("ps_18")) ps_unrare_2018_bac <- ps_18
if (!exists("ps_unrare_2019_bac") && exists("ps_19")) ps_unrare_2019_bac <- ps_19
if (!exists("ps_unrare_2020_bac") && exists("ps_20")) ps_unrare_2020_bac <- ps_20

if (!all(sapply(c("ps_unrare_2018","ps_unrare_2019","ps_unrare_2020"), exists))) {
  if (exists("process_year_fun") && all(sapply(c("main_2018","main_2020","ps_2019"), exists))) {
    cat("Rebuilding fungal ps_unrare_* ...\n")
    body(process_year_fun)[[length(body(process_year_fun))]] <- quote(list(rare = ps_rare, unrare = ps))
    .t <- process_year_fun(main_2018, "2018"); physeq_rare_2018 <- .t$rare; ps_unrare_2018 <- .t$unrare
    .t <- process_year_fun(main_2020, "2020"); physeq_rare_2020 <- .t$rare; ps_unrare_2020 <- .t$unrare
    .t <- process_year_fun(ps_2019,   "2019"); physeq_rare_2019 <- .t$rare; ps_unrare_2019 <- .t$unrare
  } else {
    stop("Fungal ps_unrare_* missing and can't rebuild (need process_year_fun + main_2018/main_2020/ps_2019).",
         "\n -> source the fungal stage script first.")
  }
}

.need <- c("ps_unrare_2018","ps_unrare_2019","ps_unrare_2020",
           "ps_unrare_2018_bac","ps_unrare_2019_bac","ps_unrare_2020_bac")
.miss <- .need[!sapply(.need, exists)]
if (length(.miss) > 0) stop("Still missing: ", paste(.miss, collapse=", "))
cat("All unrarefied objects present.\n")

# ====== THRESHOLDS ======
PADJ_CUTOFF   <- 0.05
LOG2FC_CUTOFF <- 1

# normalize stage spelling
prep_year <- function(ps) {
  sd <- as.data.frame(sample_data(ps))
  sd$Growth_Stage_Description <- gsub("Infloresence","Inflorescence",
                                      as.character(sd$Growth_Stage_Description), ignore.case = TRUE)
  sample_data(ps) <- sd
  ps
}

# ====== DESeq2 on a phyloseq object (one year, two levels) ======
run_deseq2 <- function(ps, factor_col, level1, level2) {
  sd_df <- as.data.frame(sample_data(ps))
  keep  <- sd_df[[factor_col]] %in% c(level1, level2)
  ps_sub <- prune_samples(rownames(sd_df)[keep], ps)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  
  sd_sub <- as.data.frame(sample_data(ps_sub))
  sd_sub[[factor_col]] <- factor(sd_sub[[factor_col]], levels = c(level1, level2))
  sample_data(ps_sub) <- sd_sub
  
  dds <- phyloseq_to_deseq2(ps_sub, as.formula(paste("~", factor_col)))
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- tryCatch(DESeq(dds, test="Wald", fitType="local", sfType="poscounts"),
                  error = function(e) tryCatch(DESeq(dds, test="Wald", fitType="parametric", sfType="poscounts"),
                                               error = function(e2) DESeq(dds, test="Wald", fitType="mean", sfType="poscounts")))
  
  res <- results(dds, contrast = c(factor_col, level2, level1), alpha = PADJ_CUTOFF)
  as.data.frame(res) %>%
    rownames_to_column("ASV") %>%
    mutate(
      comparison = paste0(level1, "_vs_", level2),
      factor_type = factor_col, level1 = level1, level2 = level2,
      significant = !is.na(padj) & padj < PADJ_CUTOFF & abs(log2FoldChange) > LOG2FC_CUTOFF,
      enriched_in = case_when(!significant ~ "NS",
                              log2FoldChange > 0 ~ level2,
                              log2FoldChange < 0 ~ level1, TRUE ~ "NS"))
}

# ====== volcano plot ======
make_volcano <- function(res_df, title = "", factor_type = "Management") {
  color_palette <- if (factor_type == "Management") management_colors else growth_stage_colors
  res_df <- res_df %>% mutate(neg_log10_padj = -log10(padj)) %>% filter(!is.na(padj))
  level1 <- res_df$level1[1]; level2 <- res_df$level2[1]
  n_level2 <- sum(res_df$enriched_in == level2, na.rm = TRUE)
  n_level1 <- sum(res_df$enriched_in == level1, na.rm = TRUE)
  subtitle <- paste0("\u2191", level2, ": ", n_level2, " | \u2191", level1, ": ", n_level1)
  res_df$enriched_in <- factor(res_df$enriched_in, levels = c(level1, level2, "NS"))
  
  ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_padj, color = enriched_in)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", color = "grey50") +
    scale_color_manual(values = color_palette, breaks = c(level1, level2), name = "Enriched in") +
    labs(title = title, subtitle = subtitle, x = "log2 Fold Change", y = "-log10(adjusted p-value)") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 8), legend.text = element_text(size = 7),
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 8, hjust = 0.5),
          panel.grid.minor = element_blank())
}

# ====== comparisons ======
years <- c("2018", "2019", "2020")
mgmt_comparisons <- list(c("Conventional","No-Till"), c("Conventional","Organic"), c("No-Till","Organic"))
gs_comparisons   <- list(c("Vegetative","Inflorescence"), c("Vegetative","Reproductive"), c("Inflorescence","Reproductive"))

# ====== run all comparisons — each year from its OWN unrarefied object ======
run_all_deseq <- function(kingdom_name) {
  all_results <- list(); all_volcanos <- list()
  obj_for <- function(yr) {
    nm <- if (kingdom_name == "Bacteria") paste0("ps_unrare_", yr, "_bac") else paste0("ps_unrare_", yr)
    prep_year(get(nm))
  }
  for (yr in years) {
    cat("\n=== Processing", kingdom_name, yr, "===\n")
    ps_yr <- prune_taxa(taxa_sums(obj_for(yr)) > 0, obj_for(yr))
    cat("Samples:", nsamples(ps_yr), "| ASVs:", ntaxa(ps_yr), "\n")
    
    for (comp in mgmt_comparisons) {
      cn <- paste0(yr, "_Mgmt_", comp[1], "_vs_", comp[2]); cat("  Running:", cn, "\n")
      tryCatch({
        res <- run_deseq2(ps_yr, "Management", comp[1], comp[2])
        res$Year <- yr; res$Kingdom <- kingdom_name; all_results[[cn]] <- res
        all_volcanos[[cn]] <- make_volcano(res, paste0(yr, ": ", comp[1], " vs ", comp[2]), "Management")
        cat("    Significant ASVs:", sum(res$significant), "\n")
      }, error = function(e) cat("    ERROR:", e$message, "\n"))
    }
    for (comp in gs_comparisons) {
      cn <- paste0(yr, "_GS_", comp[1], "_vs_", comp[2]); cat("  Running:", cn, "\n")
      tryCatch({
        res <- run_deseq2(ps_yr, "Growth_Stage_Description", comp[1], comp[2])
        res$Year <- yr; res$Kingdom <- kingdom_name; all_results[[cn]] <- res
        all_volcanos[[cn]] <- make_volcano(res, paste0(yr, ": ", comp[1], " vs ", comp[2]), "GrowthStage")
        cat("    Significant ASVs:", sum(res$significant), "\n")
      }, error = function(e) cat("    ERROR:", e$message, "\n"))
    }
  }
  list(results = bind_rows(all_results), volcanos = all_volcanos)
}

cat("\n", strrep("=", 50), "\nRUNNING FUNGI DESeq2\n", strrep("=", 50), "\n")
fungi_deseq    <- run_all_deseq("Fungi")
cat("\n", strrep("=", 50), "\nRUNNING BACTERIA DESeq2\n", strrep("=", 50), "\n")
bacteria_deseq <- run_all_deseq("Bacteria")

# ====== summary tables — taxonomy pulled per-ASV from each year's own object ======
build_tax_lookup <- function(kingdom_name) {
  objs <- if (kingdom_name == "Bacteria")
    list(ps_unrare_2018_bac, ps_unrare_2019_bac, ps_unrare_2020_bac)
  else list(ps_unrare_2018, ps_unrare_2019, ps_unrare_2020)
  tx <- bind_rows(lapply(objs, function(ps)
    as.data.frame(tax_table(ps)) %>% rownames_to_column("ASV")))
  distinct(tx, ASV, .keep_all = TRUE)
}

create_summary_table <- function(results_df, tax_lookup) {
  if (is.null(results_df) || nrow(results_df) == 0) { cat("No results!\n"); return(NULL) }
  sig <- results_df %>% filter(significant == TRUE) %>%
    select(ASV, Year, comparison, enriched_in, factor_type)
  if (nrow(sig) == 0) { cat("No significant results!\n"); return(NULL) }
  sig_wide <- sig %>%
    mutate(col_name = paste0(Year, "_", gsub("_vs_", "/", comparison)),
           value = paste0("\u2191", enriched_in)) %>%
    select(ASV, col_name, value) %>% distinct() %>%
    pivot_wider(names_from = col_name, values_from = value, values_fill = "-")
  tax_lookup %>% inner_join(sig_wide, by = "ASV") %>%
    arrange(Phylum, Class, Order, Family, Genus)
}

cat("\n=== Creating Summary Tables ===\n")
fungi_table    <- create_summary_table(fungi_deseq$results,    build_tax_lookup("Fungi"))
bacteria_table <- create_summary_table(bacteria_deseq$results, build_tax_lookup("Bacteria"))
cat("Fungi:",    if (is.null(fungi_table)) 0 else nrow(fungi_table),    "ASVs with >=1 significant comparison\n")
cat("Bacteria:", if (is.null(bacteria_table)) 0 else nrow(bacteria_table), "ASVs with >=1 significant comparison\n")

# ====== assemble volcano grids ======
assemble_volcano_grid <- function(volcano_list, comparison_type = "Management") {
  if (comparison_type == "Management") {
    comp_order  <- c("Mgmt_Conventional_vs_No-Till","Mgmt_Conventional_vs_Organic","Mgmt_No-Till_vs_Organic")
    comp_labels <- c("Conv vs NT","Conv vs Org","NT vs Org")
  } else {
    comp_order  <- c("GS_Vegetative_vs_Inflorescence","GS_Vegetative_vs_Reproductive","GS_Inflorescence_vs_Reproductive")
    comp_labels <- c("Veg vs Inf","Veg vs Rep","Inf vs Rep")
  }
  plot_list <- list()
  for (yr in years) for (i in seq_along(comp_order)) {
    key <- paste0(yr, "_", comp_order[i])
    plot_list[[paste0(yr, "_", i)]] <- if (key %in% names(volcano_list))
      volcano_list[[key]] + theme(legend.position = "none")
    else ggplot() + theme_void() + labs(title = paste(yr, comp_labels[i], "- No data"))
  }
  wrap_plots(plot_list, ncol = 3, nrow = 3, byrow = TRUE)
}

cat("\n=== Assembling Volcano Figures ===\n")
fungi_volcano_mgmt <- assemble_volcano_grid(fungi_deseq$volcanos, "Management") +
  plot_annotation(title = "Fungi - Management Comparisons",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)))
fungi_volcano_gs <- assemble_volcano_grid(fungi_deseq$volcanos, "GrowthStage") +
  plot_annotation(title = "Fungi - Growth Stage Comparisons",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)))
bacteria_volcano_mgmt <- assemble_volcano_grid(bacteria_deseq$volcanos, "Management") +
  plot_annotation(title = "Bacteria - Management Comparisons",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)))
bacteria_volcano_gs <- assemble_volcano_grid(bacteria_deseq$volcanos, "GrowthStage") +
  plot_annotation(title = "Bacteria - Growth Stage Comparisons",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)))

cat("\nDisplaying volcano plots...\n")
print(fungi_volcano_mgmt)
print(fungi_volcano_gs)
print(bacteria_volcano_mgmt)
print(bacteria_volcano_gs)

write.csv(fungi_table,    "DESeq2_Fungi_SignificantASVs.csv",    row.names = FALSE)
write.csv(bacteria_table, "DESeq2_Bacteria_SignificantASVs.csv", row.names = FALSE)
cat("\n=== DESeq2 Volcano Analysis Complete ===\n")

# =============================================================================
# =============================================================================
#  SECTION 15
#  Four-way analysis overlap
#
#  Supplemental Table 13.
#  Reads the CSVs written by sections 5, 12, 13 and 14, so it runs last.

RF_TOP <- 1000   # arbitrary cutoff; MUST be stated in the caption

# strip prefixes/suffixes so ids match across tables
#   bac_Zotu20_root -> Zotu20 ; Zotu54_main -> Zotu54
zid <- function(x) sub("^.*?([Zz]otu[0-9]+).*$", "\\1", as.character(x))

isa <- read.csv("isa_ALL_long.csv", stringsAsFactors = FALSE)
hub <- read.csv("hub_taxa_YrCompMgmt.csv", stringsAsFactors = FALSE)
rf  <- read.csv("Supplementary_Table_RF2_ClassifierTaxa.csv", stringsAsFactors = FALSE)
dab <- read.csv("DESeq2_Bacteria_SignificantASVs.csv", check.names = FALSE, stringsAsFactors = FALSE)
daf <- read.csv("DESeq2_Fungi_SignificantASVs.csv",   check.names = FALSE, stringsAsFactors = FALSE)

years <- c("2018", "2019", "2020")

# ---- ISA: pool tissue, keep kingdom x year -------------------
isa_set <- function(km, yr) {
  unique(zid(isa$ASV[isa$Kingdom == km & as.character(isa$Year) == yr]))
}

# ---- DA: any management contrast that year -------------------
# bacteria file writes "Organic"; fungi file writes "^Organic" with an
# arrow. Treat anything that is not "-" / blank / NA as a hit.
da_set <- function(d, yr) {
  mg   <- grep("Conventional|No-Till/Organic", names(d), value = TRUE)
  cols <- grep(paste0("^", yr), mg, value = TRUE)
  if (length(cols) == 0) return(character(0))
  vals <- as.data.frame(lapply(d[cols], function(x) trimws(as.character(x))),
                        stringsAsFactors = FALSE)
  hit  <- Reduce(`|`, lapply(vals, function(x) !is.na(x) & x != "-" & x != ""))
  unique(zid(d$ASV[hit]))
}

# ---- RF: expand Years_present, then top N by importance ------
rf$asv <- zid(rf$Representative_ASV)
rf$km  <- ifelse(grepl("Bacter", rf$Kingdom, ignore.case = TRUE), "Bacteria", "Fungi")

rf_set <- function(km, yr) {
  sub <- rf[rf$km == km & grepl(yr, rf$Years_present, fixed = TRUE), ]
  if (nrow(sub) == 0) return(character(0))
  sub <- sub[order(-sub$Mean_Importance), ]
  unique(sub$asv)[seq_len(min(RF_TOP, length(unique(sub$asv))))]
}

# ---- Hubs: pool compartment and management ------------------
hub$km <- ifelse(grepl("Bacter", hub$Kingdom, ignore.case = TRUE), "Bacteria", "Fungi")
hub_set <- function(km, yr) {
  unique(zid(hub$ASV[hub$km == km & as.character(hub$Year) == yr]))
}

# ---- assemble ------------------------------------------------
get_sets <- function(km, yr) {
  da <- if (km == "Bacteria") da_set(dab, yr) else da_set(daf, yr)
  list(ISA = isa_set(km, yr), DA = da, RF = rf_set(km, yr), Hub = hub_set(km, yr))
}

cat("\n================ SET SIZES ================\n")
for (km in c("Fungi", "Bacteria")) for (yr in years) {
  s <- get_sets(km, yr)
  cat(sprintf("%-9s %s | ISA %5d | DA %5d | RF %4d | Hub %4d | union %5d\n",
              km, yr, length(s$ISA), length(s$DA), length(s$RF), length(s$Hub),
              length(Reduce(union, s))))
}

# ---- Euler panels -------------------------------------------
set_cols <- c(ISA = "#6f8740", DA = "#bc8400", RF = "#bd461d", Hub = "#4A2C7A")

mk_panel <- function(km, yr) {
  s <- get_sets(km, yr)
  s <- s[sapply(s, length) > 0]            # drop empty sets, they break euler()
  if (length(s) < 2) return(nullGrob())
  fit <- euler(s)
  cat(sprintf("  %-9s %s  diagError = %.3f  stress = %.3f\n",
              km, yr, fit$diagError, fit$stress))
  crop <- c("2018" = "Soybean", "2019" = "Wheat", "2020" = "Maize")[yr]
  as.grob(plot(fit,
               fills      = list(fill = unname(set_cols[names(s)]), alpha = 0.55),
               labels     = list(fontsize = 12, fontface = "bold", col = "black"),
               quantities = list(fontsize = 10, fontface = "bold", col = "black"),
               main       = list(label = sprintf("%s  %s (%s)", km, crop, yr),
                                 fontsize = 14, fontface = "bold", col = "black")))
}

cat("\n============ EULER FIT DIAGNOSTICS ============\n")
cat("diagError > ~0.05 means the areas are NOT trustworthy;\n")
cat("use an UpSet plot instead if that happens.\n\n")

panels <- list()
for (km in c("Fungi", "Bacteria")) for (yr in years) {
  panels[[length(panels) + 1]] <- mk_panel(km, yr)
}

overlap_figure <- arrangeGrob(grobs = panels, nrow = 2, ncol = 3)

grid.newpage(); grid.draw(overlap_figure)

pdf("SupFig_AnalysisOverlap.pdf", width = 15, height = 10)
grid.newpage(); grid.draw(overlap_figure)
dev.off()



# ─────────────────────────────────────────────────────────────
# SUPPLEMENTAL TABLE — ASV overlap across the four analyses
#   ISA | DESeq2 (DA) | Random Forest | Network hubs
# Pooled across compartments (DESeq2 has no compartment field):
#   2 kingdoms x 3 years = 6 strata
#
# Writes two CSVs:
#   overlap_intersections.csv — every non-empty combination
#   overlap_summary.csv       — per method: set size, unique, shared
# ─────────────────────────────────────────────────────────────

RF_TOP <- Inf   # no cap; the table below confirms it never binds

zid <- function(x) sub("^.*?([Zz]otu[0-9]+).*$", "\\1", as.character(x))

isa <- read.csv("isa_ALL_long.csv", stringsAsFactors = FALSE)
hub <- read.csv("hub_taxa_YrCompMgmt.csv", stringsAsFactors = FALSE)
rf  <- read.csv("Supplementary_Table_RF2_ClassifierTaxa.csv", stringsAsFactors = FALSE)
dab <- read.csv("DESeq2_Bacteria_SignificantASVs.csv", check.names = FALSE, stringsAsFactors = FALSE)
daf <- read.csv("DESeq2_Fungi_SignificantASVs.csv",   check.names = FALSE, stringsAsFactors = FALSE)

years <- c("2018", "2019", "2020")
crops <- c("2018" = "Soybean", "2019" = "Wheat", "2020" = "Maize")

isa_set <- function(km, yr)
  unique(zid(isa$ASV[isa$Kingdom == km & as.character(isa$Year) == yr]))

da_set <- function(d, yr) {
  mg   <- grep("Conventional|No-Till/Organic", names(d), value = TRUE)
  cols <- grep(paste0("^", yr), mg, value = TRUE)
  if (length(cols) == 0) return(character(0))
  vals <- as.data.frame(lapply(d[cols], function(x) trimws(as.character(x))),
                        stringsAsFactors = FALSE)
  hit  <- Reduce(`|`, lapply(vals, function(x) !is.na(x) & x != "-" & x != ""))
  unique(zid(d$ASV[hit]))
}

rf$asv <- zid(rf$Representative_ASV)
rf$km  <- ifelse(grepl("Bacter", rf$Kingdom, ignore.case = TRUE), "Bacteria", "Fungi")
rf_set <- function(km, yr) {
  sub <- rf[rf$km == km & grepl(yr, rf$Years_present, fixed = TRUE), ]
  if (nrow(sub) == 0) return(character(0))
  sub <- sub[order(-sub$Mean_Importance), ]
  u <- unique(sub$asv)
  u[seq_len(min(RF_TOP, length(u)))]
}

hub$km <- ifelse(grepl("Bacter", hub$Kingdom, ignore.case = TRUE), "Bacteria", "Fungi")
hub_set <- function(km, yr)
  unique(zid(hub$ASV[hub$km == km & as.character(hub$Year) == yr]))

get_sets <- function(km, yr) {
  da <- if (km == "Bacteria") da_set(dab, yr) else da_set(daf, yr)
  list(ISA = isa_set(km, yr), DA = da, RF = rf_set(km, yr), Hub = hub_set(km, yr))
}

methods <- c("ISA", "DA", "RF", "Hub")

# ---- TABLE 1: intersection counts ---------------------------
# For every ASV, record which of the four sets it belongs to, then
# tabulate the resulting membership patterns. Every ASV lands in
# exactly one pattern, so the counts sum to the union.
intersections <- do.call(rbind, lapply(c("Fungi", "Bacteria"), function(km) {
  do.call(rbind, lapply(years, function(yr) {
    s   <- get_sets(km, yr)
    all <- Reduce(union, s)
    if (length(all) == 0) return(NULL)
    memb <- sapply(methods, function(m) all %in% s[[m]])
    pat  <- apply(memb, 1, function(r) paste(methods[r], collapse = " + "))
    tb   <- as.data.frame(table(pat), stringsAsFactors = FALSE)
    names(tb) <- c("Combination", "n_ASV")
    tb$n_methods <- lengths(strsplit(tb$Combination, " \\+ "))
    tb <- tb[order(-tb$n_methods, -tb$n_ASV), ]
    data.frame(Kingdom = km, Year = yr, Crop = unname(crops[yr]),
               Combination = tb$Combination,
               n_methods   = tb$n_methods,
               n_ASV       = tb$n_ASV,
               pct_of_union = round(100 * tb$n_ASV / length(all), 1),
               union_total = length(all),
               stringsAsFactors = FALSE)
  }))
}))

# ---- TABLE 2: per-method summary ----------------------------
# How big is each set, how much is unique to it, how much is shared.
summary_tbl <- do.call(rbind, lapply(c("Fungi", "Bacteria"), function(km) {
  do.call(rbind, lapply(years, function(yr) {
    s <- get_sets(km, yr)
    do.call(rbind, lapply(methods, function(m) {
      others <- Reduce(union, s[setdiff(methods, m)])
      uniq   <- setdiff(s[[m]], others)
      data.frame(Kingdom = km, Year = yr, Crop = unname(crops[yr]),
                 Method = m,
                 set_size    = length(s[[m]]),
                 unique_to_method = length(uniq),
                 shared      = length(s[[m]]) - length(uniq),
                 pct_shared  = if (length(s[[m]]) == 0) NA else
                   round(100 * (length(s[[m]]) - length(uniq)) / length(s[[m]]), 1),
                 stringsAsFactors = FALSE)
    }))
  }))
}))

write.csv(intersections, "overlap_intersections.csv", row.names = FALSE)
write.csv(summary_tbl,   "overlap_summary.csv",       row.names = FALSE)

cat("\n===== PER-METHOD SUMMARY =====\n")
print(summary_tbl, row.names = FALSE)

cat("\n===== INTERSECTIONS (3+ methods only, for readability) =====\n")
print(intersections[intersections$n_methods >= 3, ], row.names = FALSE)

cat("\nWrote overlap_intersections.csv (", nrow(intersections), " rows)\n", sep = "")
cat("Wrote overlap_summary.csv (", nrow(summary_tbl), " rows)\n", sep = "")

# ---- sanity: does RF_TOP ever bind? -------------------------
cat("\n===== RF cap check =====\n")
for (km in c("Fungi","Bacteria")) for (yr in years) {
  sub <- rf[rf$km == km & grepl(yr, rf$Years_present, fixed = TRUE), ]
  cat(sprintf("%-9s %s | available %5d | used %5d\n",
              km, yr, length(unique(sub$asv)), length(rf_set(km, yr))))
}




# ─────────────────────────────────────────────────────────────
# SUPPLEMENTAL TABLES — ASV overlap across the four analyses
#   ISA | DESeq2 (DA) | Random Forest | Network hubs
# Pooled across compartments (DESeq2 has no compartment field):
#   2 kingdoms x 3 years = 6 strata
#
# Writes three CSVs:
#   Table_S_overlap_summary.csv       — per method: size, unique, shared
#   Table_S_overlap_intersections.csv — every membership combination
#   Table_S_overlap_shared_ASVs.csv   — the shared ASVs themselves, named
# ─────────────────────────────────────────────────────────────

RF_TOP <- Inf   # no cap; the check at the bottom confirms it never binds

zid <- function(x) sub("^.*?([Zz]otu[0-9]+).*$", "\\1", as.character(x))

isa <- read.csv("isa_ALL_long.csv", stringsAsFactors = FALSE)
hub <- read.csv("hub_taxa_YrCompMgmt.csv", stringsAsFactors = FALSE)
rf  <- read.csv("Supplementary_Table_RF2_ClassifierTaxa.csv", stringsAsFactors = FALSE)
dab <- read.csv("DESeq2_Bacteria_SignificantASVs.csv", check.names = FALSE, stringsAsFactors = FALSE)
daf <- read.csv("DESeq2_Fungi_SignificantASVs.csv",   check.names = FALSE, stringsAsFactors = FALSE)

years <- c("2018", "2019", "2020")
crops <- c("2018" = "Soybean", "2019" = "Wheat", "2020" = "Maize")

# ---- taxonomy lookup, so shared ASVs get names not just ids ---
tax_lookup <- local({
  from_isa <- data.frame(
    asv  = zid(isa$ASV),
    name = trimws(sub("\\s*\\(.*\\)$", "", isa$ASV_label)),   # strip " (ZotuN)"
    stringsAsFactors = FALSE)
  
  best <- function(d, idcol, ranks) {
    nm <- apply(d[ranks], 1, function(r) {
      r <- r[!is.na(r) & trimws(r) != "" & trimws(r) != "-"]
      if (length(r) == 0) NA_character_ else r[1]
    })
    data.frame(asv = zid(d[[idcol]]), name = nm, stringsAsFactors = FALSE)
  }
  from_dab <- best(dab, "ASV", c("Genus","Family","Order","Class","Phylum"))
  from_daf <- best(daf, "ASV", c("Genus","Family","Order","Class","Phylum"))
  from_hub <- best(hub, "ASV", c("Genus","Family","Phylum"))
  
  all <- rbind(from_isa, from_dab, from_daf, from_hub)
  all <- all[!is.na(all$name) & all$name != "", ]
  all <- all[!duplicated(all$asv), ]
  setNames(all$name, all$asv)
})
tax_of <- function(a) ifelse(a %in% names(tax_lookup), unname(tax_lookup[a]), "Unclassified")

# ---- set builders --------------------------------------------
isa_set <- function(km, yr)
  unique(zid(isa$ASV[isa$Kingdom == km & as.character(isa$Year) == yr]))

da_set <- function(d, yr) {
  mg   <- grep("Conventional|No-Till/Organic", names(d), value = TRUE)
  cols <- grep(paste0("^", yr), mg, value = TRUE)
  if (length(cols) == 0) return(character(0))
  vals <- as.data.frame(lapply(d[cols], function(x) trimws(as.character(x))),
                        stringsAsFactors = FALSE)
  hit  <- Reduce(`|`, lapply(vals, function(x) !is.na(x) & x != "-" & x != ""))
  unique(zid(d$ASV[hit]))
}

rf$asv <- zid(rf$Representative_ASV)
rf$km  <- ifelse(grepl("Bacter", rf$Kingdom, ignore.case = TRUE), "Bacteria", "Fungi")
rf_set <- function(km, yr) {
  sub <- rf[rf$km == km & grepl(yr, rf$Years_present, fixed = TRUE), ]
  if (nrow(sub) == 0) return(character(0))
  sub <- sub[order(-sub$Mean_Importance), ]
  u <- unique(sub$asv)
  u[seq_len(min(RF_TOP, length(u)))]
}

hub$km <- ifelse(grepl("Bacter", hub$Kingdom, ignore.case = TRUE), "Bacteria", "Fungi")
hub_set <- function(km, yr)
  unique(zid(hub$ASV[hub$km == km & as.character(hub$Year) == yr]))

get_sets <- function(km, yr) {
  da <- if (km == "Bacteria") da_set(dab, yr) else da_set(daf, yr)
  list(`Indicator species` = isa_set(km, yr),
       `Differential abundance` = da,
       `Random forest` = rf_set(km, yr),
       `Network hub` = hub_set(km, yr))
}
methods <- c("Indicator species", "Differential abundance", "Random forest", "Network hub")

# TABLE 1 — per-method summary
summary_tbl <- do.call(rbind, lapply(c("Fungi", "Bacteria"), function(km) {
  do.call(rbind, lapply(years, function(yr) {
    s <- get_sets(km, yr)
    do.call(rbind, lapply(methods, function(m) {
      others <- Reduce(union, s[setdiff(methods, m)])
      uniq   <- setdiff(s[[m]], others)
      n      <- length(s[[m]])
      data.frame(
        Kingdom                  = km,
        `Crop year`              = paste0(unname(crops[yr]), " (", yr, ")"),
        Analysis                 = m,
        `ASVs recovered`         = n,
        `Shared with other analyses` = n - length(uniq),
        `Unique to this analysis`    = length(uniq),
        `Percent shared`         = if (n == 0) NA else round(100 * (n - length(uniq)) / n, 1),
        check.names = FALSE, stringsAsFactors = FALSE)
    }))
  }))
}))

# TABLE 2 — intersection counts
intersections <- do.call(rbind, lapply(c("Fungi", "Bacteria"), function(km) {
  do.call(rbind, lapply(years, function(yr) {
    s   <- get_sets(km, yr)
    all <- Reduce(union, s)
    if (length(all) == 0) return(NULL)
    memb <- sapply(methods, function(m) all %in% s[[m]])
    pat  <- apply(memb, 1, function(r) paste(methods[r], collapse = " + "))
    tb   <- as.data.frame(table(pat), stringsAsFactors = FALSE)
    names(tb) <- c("Combination", "n")
    tb$k <- lengths(strsplit(tb$Combination, " \\+ "))
    tb   <- tb[order(-tb$k, -tb$n), ]
    data.frame(
      Kingdom                = km,
      `Crop year`            = paste0(unname(crops[yr]), " (", yr, ")"),
      `Analyses in agreement` = tb$Combination,
      `Number of analyses`   = tb$k,
      `ASVs`                 = tb$n,
      `Percent of total`     = round(100 * tb$n / length(all), 1),
      `Total ASVs recovered` = length(all),
      check.names = FALSE, stringsAsFactors = FALSE)
  }))
}))

# TABLE 3 — the shared ASVs themselves
# Any ASV recovered by two or more analyses, with which ones and its finest available taxonomic assignment.
shared_asvs <- do.call(rbind, lapply(c("Fungi", "Bacteria"), function(km) {
  do.call(rbind, lapply(years, function(yr) {
    s   <- get_sets(km, yr)
    all <- Reduce(union, s)
    if (length(all) == 0) return(NULL)
    memb <- sapply(methods, function(m) all %in% s[[m]])
    k    <- rowSums(memb)
    keep <- k >= 2
    if (!any(keep)) return(NULL)
    out <- data.frame(
      Kingdom                = km,
      `Crop year`            = paste0(unname(crops[yr]), " (", yr, ")"),
      ASV                    = all[keep],
      Taxon                  = tax_of(all[keep]),
      `Number of analyses`   = k[keep],
      `Analyses in agreement` = apply(memb[keep, , drop = FALSE], 1,
                                      function(r) paste(methods[r], collapse = " + ")),
      check.names = FALSE, stringsAsFactors = FALSE)
    out[order(-out$`Number of analyses`, out$Taxon), ]
  }))
}))


write.csv(summary_tbl,   "Table_S_overlap_summary.csv",       row.names = FALSE)
write.csv(intersections, "Table_S_overlap_intersections.csv", row.names = FALSE)
write.csv(shared_asvs,   "Table_S_overlap_shared_ASVs.csv",   row.names = FALSE)

cat("\n===== PER-ANALYSIS SUMMARY =====\n")
print(summary_tbl, row.names = FALSE)

cat("\n===== ASVs FOUND BY ALL FOUR ANALYSES =====\n")
four <- shared_asvs[shared_asvs$`Number of analyses` == 4, ]
if (nrow(four)) print(four[, c("Kingdom","Crop year","ASV","Taxon")], row.names = FALSE) else
  cat("(none)\n")

cat("\nWrote Table_S_overlap_summary.csv (", nrow(summary_tbl), " rows)\n", sep = "")
cat("Wrote Table_S_overlap_intersections.csv (", nrow(intersections), " rows)\n", sep = "")
cat("Wrote Table_S_overlap_shared_ASVs.csv (", nrow(shared_asvs), " rows)\n", sep = "")

# ---- RF cap check --------------------------------------------
cat("\n===== RF cap check (available vs used) =====\n")
for (km in c("Fungi","Bacteria")) for (yr in years) {
  sub <- rf[rf$km == km & grepl(yr, rf$Years_present, fixed = TRUE), ]
  cat(sprintf("%-9s %s | available %5d | used %5d\n",
              km, yr, length(unique(sub$asv)), length(rf_set(km, yr))))
  

# SUPPLEMENTAL TABLE — one row per stratum, wide format
  
  overlap_wide <- do.call(rbind, lapply(c("Fungi", "Bacteria"), function(km) {
    do.call(rbind, lapply(years, function(yr) {
      s   <- get_sets(km, yr)
      all <- Reduce(union, s)
      # membership matrix: rows = ASVs in the union, cols = the 4 analyses
      memb <- sapply(methods, function(m) all %in% s[[m]])
      k    <- rowSums(memb)          # how many analyses recovered each ASV
      
      data.frame(
        Kingdom                 = km,
        `Crop year`             = paste0(unname(crops[yr]), " (", yr, ")"),
        `Indicator species`     = length(s[["Indicator species"]]),
        `Differential abundance`= length(s[["Differential abundance"]]),
        `Random forest`         = length(s[["Random forest"]]),
        `Network hub`           = length(s[["Network hub"]]),
        `Total ASVs recovered`  = length(all),
        `Shared by 2 or more`   = sum(k >= 2),
        `Shared by 3 or more`   = sum(k >= 3),
        `Shared by all four`    = sum(k == 4),
        `Percent shared`        = round(100 * sum(k >= 2) / length(all), 1),
        check.names = FALSE, stringsAsFactors = FALSE)
    }))
  }))
  
  write.csv(overlap_wide, "Table_S_overlap_wide.csv", row.names = FALSE)
  
  cat("\n===== OVERLAP SUMMARY (wide) =====\n")
  print(overlap_wide, row.names = FALSE)
  cat("\nWrote Table_S_overlap_wide.csv (", nrow(overlap_wide), " rows)\n", sep = "")
}


