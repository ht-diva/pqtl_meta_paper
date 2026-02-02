#Question Claudia
#why multigenes removal from all panel (including LOD) and including pQTL-> add them as missing instead of removing them (make sure that legend is clear)
#why remove NA in below LOD -> wrong percentage -> solved
#why constraining scale for CV plots? -> currently log scale of percent... ->solved % for sf1  log raw value for the current ones (sf4 sf7)
#below LOD pattern
#sum of prot in panel d e f is not equal to the ones from other panels-> remove them OK
##problem of n in models?
##binomial model (outcome always 1 0) ->correctly specified, just have to change N in the reporting table
#re-add models with QC metrics as outcome (+ heterogenety) -> poisson/ZIP




############################################################
## 0. Libraries, helpers, paths
############################################################
# options: make_figure1, make_figure2, run_models
make_figure1 = TRUE
make_figure2 = TRUE
run_models   = TRUE
add_complexity_models = TRUE

# output file for models
models_file = "Supplementary_Table_Regressions_cis_discovery_and_complexity.csv"

# input: need these files in wd:
# supplementary_table_2.xlsx
# supplementary_table_3.xlsx
# supplementary_table_4.xlsx
# setwd("/Users/c.giambartolomei/Downloads/macbook_change_jan2026/interval_chris_manuscript/annotations_of_proteins_oct222025/manuscript_tables/")
setwd("/group/diangelantonio/users/Solene/2026_meta_files/")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(tidyr)
  library(cowplot)
  library(RColorBrewer)
})



## Helper: mode for categorical annotations
unique_or_na <- function(x) {
  ux <- unique(na.omit(x))
  if (length(ux) == 0) return(NA_character_)
  if (length(ux) > 1) {
    stop("Multiple conflicting values found: ", paste(ux, collapse = ", "))
  }
  ux
}

## SeqIDs may target multiple proteins; a representative UniProt ID is selected
## for bookkeeping and plotting purposes only.
representative_uniprot <- function(x) {
  ux <- na.omit(x)
  if (length(ux) == 0) return(NA_character_)
  names(sort(table(ux), decreasing = TRUE))[1]
}

protein_panel_colors <- c(
  "Present in 5k" = "#404080",
  "New in 7k"     = "#69b3a2"
)

############################################################
## 1. Load input tables
############################################################

annot_raw <- read_excel("supplementary_table_2.xlsx", sheet = 2, skip = 1)
##########
##adding real number of samples below LOD
annot_raw$INTERVAL_Number_under_LOD_batch1<-round(annot_raw$INTERVAL_Percent_under_LOD_batch1*4933/100)
annot_raw$INTERVAL_Number_under_LOD_batch2<-round(annot_raw$INTERVAL_Percent_under_LOD_batch2*4724/100)
annot_raw$CHRIS_Number_under_LOD<-round(annot_raw$CHRIS_Percent_under_LOD*4204/100)
###################
pqtl_raw  <- read_excel("supplementary_table_3.xlsx", sheet = "ST3", skip = 1)
dim(annot_raw) #7391   48
table(duplicated(annot_raw$SeqID))
# FALSE  TRUE 
#7289   102 
##duplicated seqID -> uncollapsed table
#7289 matches number of human protein in the assay

## Clean column names
names(annot_raw) <- make.names(names(annot_raw))
names(pqtl_raw)  <- make.names(names(pqtl_raw))

## Harmonise SeqID
# annot_raw <- annot_raw %>% rename(SeqID = SeqId)

## Aptamer new vs old (7k vs 5k)
annot <- annot_raw %>%
  mutate(
    aptamer_new = !aptamer_in_somascan5k
  )

## ---------------------------------------------------------
## Diagnostics: mapping multiplicity
## ---------------------------------------------------------
# Protein-level annotations (must be unique)
# Gene-level annotations (remove 4 cases where A UniProt ID maps to >1 HGNC symbol)
# Aptamer bookkeeping (representative mapping OK)

cat("SeqIDs mapping to >1 UniProt:",
    annot %>% distinct(SeqID, UniProt_ID) %>% count(SeqID) %>% filter(n > 1) %>% nrow(),
    "\n")

cat("UniProt IDs mapping to >1 RNA tissue specificity:",
    annot %>% distinct(UniProt_ID, RNA.tissue.specificity) %>%
      filter(!is.na(RNA.tissue.specificity), RNA.tissue.specificity != "Not detected") %>%
      count(UniProt_ID) %>% filter(n > 1) %>% nrow(),
    "\n")

## Identify UniProt IDs mapping to multiple HGNC symbols
protein_multigene <- annot %>%
  distinct(UniProt_ID, HGNC_Symbol) %>%
  group_by(UniProt_ID) %>%
  summarise(
    n_genes = n_distinct(HGNC_Symbol),
    has_multigene = n_genes > 1,
    .groups = "drop"
  )

cat("Proteins with >1 HGNC symbol:", sum(protein_multigene$has_multigene), "\n")

## Identify UniProt IDs with RNA annotation conflicts
protein_rna_conflict <- annot %>%
  distinct(UniProt_ID, RNA.tissue.distribution, RNA.tissue.specificity) %>%
  group_by(UniProt_ID) %>%
  summarise(
    n_dist = n_distinct(RNA.tissue.distribution, na.rm = TRUE),
    n_spec = n_distinct(RNA.tissue.specificity,  na.rm = TRUE),
    has_rna_conflict = (n_dist > 1 | n_spec > 1),
    .groups = "drop"
  )
protein_flags <- protein_multigene %>%
  left_join(protein_rna_conflict, by = "UniProt_ID") %>%
  mutate(
    has_rna_conflict = replace_na(has_rna_conflict, FALSE)
  )

table(
  multigene    = protein_flags$has_multigene,
  rna_conflict = protein_flags$has_rna_conflict
)

proteins_to_exclude <- protein_flags %>%
  filter(has_multigene) %>%   # ← this is the decisive rule
  pull(UniProt_ID)

cat(
  "Total proteins excluded as multi-gene complexes:",
  length(proteins_to_exclude),
  "\n"
)
#4
annot_prot <- annot %>%
  filter(!UniProt_ID %in% proteins_to_exclude)

# 4 multigene proteins, 2 of them also RNA-conflicting
View(annot_raw[annot_raw$UniProt_ID%in%protein_multigene$UniProt_ID[protein_multigene$has_multigene==TRUE],]) #View excluded raws
unique(annot_raw[annot_raw$UniProt_ID%in%protein_multigene$UniProt_ID[protein_multigene$has_multigene==TRUE],]$SeqID) #View excluded raws
#4 seqID
dim(annot) #7375 49 #16 raws excluded
#
## ---------------------------------------------------------
## CLEAN RNA ANNOTATIONS (before collapsing)
## ---------------------------------------------------------
table(annot_prot$Secretion.pathway)
table(annot_prot$RNA.tissue.specificity)
table(annot_prot$RNA.tissue.distribution)
# Fix this: RNA tissue expression annotations from the Human Protein Atlas occasionally appeared as missing (‘Not detected’) in Supplementary Table 2
annot_prot <- annot_prot %>%
  mutate(
    RNA.tissue.specificity = ifelse(
      RNA.tissue.specificity %in% c("Not detected", "NA", "", NA),
      NA,
      RNA.tissue.specificity
    ),
    RNA.tissue.distribution = ifelse(
      RNA.tissue.distribution %in% c("Not detected", "NA", "", NA),
      NA,
      RNA.tissue.distribution
    )
  )

############################################################
## 2. Protein-level summarisation
##    One row per UniProt_ID
############################################################
df_uniprot1 <- annot_prot %>%
  group_by(UniProt_ID) %>%
  summarise(
    HGNC_Symbol             = unique_or_na(HGNC_Symbol),
    Secretion.pathway       = unique_or_na(Secretion.pathway),
    RNA.tissue.distribution = unique_or_na(`RNA.tissue.distribution`),
    RNA.tissue.specificity  = unique_or_na(`RNA.tissue.specificity`),
    Blood_conc_pgL          = median(as.numeric(Blood_conc_pgL), na.rm = TRUE),
    uniprot_new_in_somascan7k_vs5k =
      any(uniprot_new_in_somascan7k_vs5k, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    group_lbl = factor(
      ifelse(uniprot_new_in_somascan7k_vs5k, "New in 7k", "Present in 5k"),
      levels = c("Present in 5k", "New in 7k")
    )
  )

############################################################
## 3. pQTL expansion (for Figure 2 cis/no-cis stratification)
############################################################

pqtl_expanded <- pqtl_raw %>%
  separate_rows(UniProt_ID, sep = "\\|") %>%
  mutate(UniProt_ID = trimws(UniProt_ID)) %>%
  filter(UniProt_ID != "")

cis_protein <- pqtl_expanded %>%
  filter(cis_or_trans == "cis") %>%
  distinct(UniProt_ID) %>%
  mutate(has_cis = TRUE)

pqtl_status_prot <- df_uniprot1 %>%
  left_join(cis_protein, by = "UniProt_ID") %>%
  mutate(
    has_cis = ifelse(is.na(has_cis), FALSE, has_cis),
    pqtl_group = factor(ifelse(has_cis, "Has ≥1 cis-pQTL", "No cis-pQTL"),
                        levels = c("No cis-pQTL", "Has ≥1 cis-pQTL"))
  )

############################################################
## 4. Aptamer-level dataset (for QC panels)
##    One row per SeqID, QC per aptamer × batch
############################################################
df_aptamer <- annot %>%
  group_by(SeqID) %>%
  summarise(
    UniProt_ID  = representative_uniprot(UniProt_ID),
    HGNC_Symbol = representative_uniprot(HGNC_Symbol),
    aptamer_new = any(aptamer_new),

    # one value per SeqID × batch (take the first non-NA) 
    INTERVAL_LOD_batch1               = dplyr::first(na.omit(INTERVAL_LOD_batch1)),
    INTERVAL_LOD_batch2               = dplyr::first(na.omit(INTERVAL_LOD_batch2)),
    INTERVAL_Percent_under_LOD_batch1 = dplyr::first(na.omit(INTERVAL_Percent_under_LOD_batch1)),
    INTERVAL_Percent_under_LOD_batch2 = dplyr::first(na.omit(INTERVAL_Percent_under_LOD_batch2)),
    INTERVAL_Number_under_LOD_batch1 = dplyr::first(na.omit(INTERVAL_Number_under_LOD_batch1)),
    INTERVAL_Number_under_LOD_batch2 = dplyr::first(na.omit(INTERVAL_Number_under_LOD_batch2)),
    INTERVAL_CV_batch1                = dplyr::first(na.omit(INTERVAL_CV_batch1)),
    INTERVAL_CV_batch2                = dplyr::first(na.omit(INTERVAL_CV_batch2)),
    CHRIS_LOD                         = dplyr::first(na.omit(CHRIS_LOD)),
    CHRIS_Percent_under_LOD           = dplyr::first(na.omit(CHRIS_Percent_under_LOD)),
    CHRIS_CV                          = dplyr::first(na.omit(CHRIS_CV)),
    CHRIS_Number_under_LOD = dplyr::first(na.omit(CHRIS_Number_under_LOD)),
    .groups = "drop"
  ) %>%
  mutate(
    group_lbl = factor(
      ifelse(aptamer_new, "New in 7k", "Present in 5k"),
      levels = c("Present in 5k", "New in 7k")
    )
  )
dim(df_aptamer) #7285

## Join cis/no-cis for aptamers (used for Figure 2)
cis_apt <- pqtl_raw %>%
  filter(cis_or_trans == "cis") %>%
  distinct(SeqID) %>%
  mutate(has_cis = TRUE)

df_aptamer2 <- df_aptamer %>%
  left_join(cis_apt, by = "SeqID") %>%
  mutate(
    has_cis = ifelse(is.na(has_cis), FALSE, has_cis),
    pqtl_group = factor(ifelse(has_cis, "Has ≥1 cis-pQTL", "No cis-pQTL"),
                        levels = c("No cis-pQTL", "Has ≥1 cis-pQTL"))
  )

############################################################
## PART 2 — Unified QC system (APPLIES TO ALL PANELS D/E/F)
############################################################

## ---------------------------------------------------------
## 2.1 QC value transformations
## ---------------------------------------------------------

transform_LOD <- function(x) {
  x <- as.numeric(x)
  x <- pmax(x, 1e-4) #to allow log transfor for 0 
  log10(x)
}

transform_CV <- function(x) {
  as.numeric(x) * 100     # convert fraction → %
}

transform_pctLOD <- function(x) {
  as.numeric(x)
}

## ---------------------------------------------------------
## 2.2 UNIFIED RESHAPE ENGINE FOR ALL QC PANELS
##     (This fixes all CHRIS/INTERVAL mixing issues)
## ---------------------------------------------------------

reshape_QC <- function(df, batches, transform_fun = identity) {

  long_df <- bind_rows(
    lapply(names(batches), function(b) {
      tibble(
        Batch      = b,
        SeqID      = df$SeqID,
        UniProt_ID = df$UniProt_ID,
        group_lbl  = df$group_lbl,
        value_raw  = df[[ batches[[b]] ]]
      )
    })
  ) %>%
    # filter(!is.na(value_raw)) %>%      # IMPORTANT: remove batch-specific NA-->?
    #Commented to prevent removing of batch specific NA: LOD should not be NA, Cv as well, and below LOD NA corresponds to 0%
    mutate(
      value = transform_fun(value_raw)
    )

  return(long_df)
}

## ---------------------------------------------------------
## 2.3 Unified QC violin plotter (now simplified)
## ---------------------------------------------------------

batch_colors <- c(
  "INTERVAL Batch 1" = "#1f78b4",
  "INTERVAL Batch 2" = "#a6cee3",
  "CHRIS"            = "#ff7f00"
)


make_qc_violin <- function(long_df, ylab, y_limits = NULL,y_scale="linear") {

  # Per-batch Ns (apt first)
  counts <- long_df %>%
    group_by(group_lbl, Batch) %>%
    summarise(
      aptamers = n_distinct(SeqID),
      proteins = n_distinct(UniProt_ID),
      .groups  = "drop"
    ) %>%
     mutate(label = paste0(aptamers, " aptamers\n", proteins, " proteins")) 

  counts_center <- counts %>% filter(Batch == "INTERVAL Batch 2") # print only middle since they are all same (non-missing CV and LOD across all cohorts) ??

  p <- ggplot(long_df, aes(x = group_lbl, y = value, fill = Batch)) +
    geom_violin(alpha = 0.5, scale = "area",
                position = position_dodge(width = 0.9)) +
    geom_boxplot(width = 0.12, outlier.size = 0.4,
                 position = position_dodge(width = 0.9)) +
    geom_text(
      data = counts_center,
      aes(x = group_lbl, y = Inf, label = label, group = Batch),
      position    = position_dodge(width = 0.9),
      vjust       = 1.1,
      size        = 3,
      fontface    = "italic",
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = batch_colors) +
    labs(x = NULL, y = ylab) +
    theme_bw(14) +
    theme(legend.position = "none")
  
  if (y_scale == "log10") {
    p <- p +
      scale_y_log10(
        expand = expansion(mult = c(0.01, 0.10))
      )
  } else {
    p <- p +
      scale_y_continuous(
        limits = y_limits,
        expand = expansion(mult = c(0.01, 0.10))
      )
  }
  return(p)
}

## 1.1 Stacked bar helper (A–C), 100% stacked to 100
make_panel1 <- function(df, feature, y_max, valid_levels) {

  df_filtered <- df %>%
    group_by(UniProt_ID, HGNC_Symbol, group_lbl) %>%
    summarise(
      feature_val = .data[[feature]],
      .groups = "drop"
    ) %>%
    filter(!is.na(feature_val), feature_val %in% valid_levels) %>%
    mutate(feature_val = factor(feature_val, levels = valid_levels))

  panel_data <- df_filtered %>%
    count(group_lbl, feature_val) %>%
    group_by(group_lbl) %>%
    mutate(
      total   = sum(n),
      percent = n / total * 100,
      label   = paste0(n, " (", round(percent, 1), "%)")
    )

  summary_counts <- df_filtered %>%
    group_by(group_lbl) %>%
    summarise(
      proteins_with_data = n_distinct(UniProt_ID),
      proteins_total     = n_distinct(df$UniProt_ID[df$group_lbl == unique(group_lbl)]),
      proteins_missing   = proteins_total - proteins_with_data,
      label = paste0(
        proteins_with_data, " proteins\n(",
        proteins_missing, " missing)"
      ),
      .groups = "drop"
    )
#      label    = paste0(proteins, " proteins\n", genes, " genes"),


  ggplot(panel_data, aes(x = group_lbl, y = percent, fill = feature_val)) +
    geom_col(position = "stack", alpha = 0.9) +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 3.5, fontface = "bold", color = "white"
    ) +
    geom_text(
      data = summary_counts,
      aes(x = group_lbl, y = y_max, label = label),
      vjust = 0,
      inherit.aes = FALSE,
      size = 3,
      fontface = "italic"
    ) +
    scale_y_continuous(
  	limits = c(0, y_max * 1.02),
  	breaks = c(0, 25, 50, 75, 100),
  	expand = expansion(mult = c(0, 0.15))) +
    scale_fill_brewer(palette = "Set2", name = NULL) +
    labs(x = NULL, y = "Percentage of proteins") +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom")
}

if (make_figure1) {

  ## 1.2 Panels A–C (biology)
  p1_a <- make_panel1(
    df_uniprot1,
    feature      = "Secretion.pathway",
    y_max        = 110,
    valid_levels = c("intracellular", "membrane", "secreted")
  )

  p1_b <- make_panel1(
    df_uniprot1,
    feature      = "RNA.tissue.distribution",
    y_max        = 105,
    valid_levels = c("Detected in all tissues",
                     "Detected in many",
                     "Detected in some",
                     "Detected in single")
  )

  p1_c <- make_panel1(
    df_uniprot1,
    feature      = "RNA.tissue.specificity",
    y_max        = 105,
    valid_levels = c("Group enriched",
                     "Low tissue specificity",
                     "Tissue enhanced",
                     "Tissue enriched")
  )

  ## 1.3 Blood concentration panel G
  df_blood <- df_uniprot1 %>%
    filter(!is.na(Blood_conc_pgL), Blood_conc_pgL > 0)

  summary_counts_g <- df_uniprot1 %>%
    group_by(group_lbl) %>%
    summarise(
      proteins_total   = n_distinct(UniProt_ID),
      genes            = n_distinct(HGNC_Symbol),
      proteins_nonmiss = n_distinct(UniProt_ID[!is.na(Blood_conc_pgL) & Blood_conc_pgL > 0]),
      proteins_missing = proteins_total - proteins_nonmiss,
      label = paste0(
        proteins_nonmiss, " proteins\n(",
        proteins_missing, " missing)"
      ),
    .groups = "drop"
  )
#        "N=", proteins_nonmiss, " proteins\n(",
#        genes, " genes; ", proteins_missing, " missing)"


  p1_g <- ggplot(df_blood,
                 aes(x = group_lbl, y = Blood_conc_pgL, fill = group_lbl)) +
    geom_violin(alpha = 0.4, scale = "area") +
    geom_boxplot(width = 0.15, outlier.size = 0.4) +
    geom_text(
      data = summary_counts_g,
      aes(x = group_lbl, y = Inf, label = label),
      vjust = 1.2,
      size = 3,
      fontface = "italic",
      inherit.aes = FALSE
    ) +
    scale_y_log10(
      expand = expansion(mult = c(0.02, 0.20))
    ) +
    scale_fill_manual(values = protein_panel_colors, name = NULL) +
    labs(x = NULL, y = "Blood conc. (log10 pg/L)") +
    theme_bw(14) +
    theme(
      legend.position = "none",
      axis.title.y = element_text(margin = margin(r = 30))
    )

  ##########################################################
  ## 2.  QC panels (D–F) – per-batch QC
  ##########################################################
  ## D – LOD

long_LOD <- reshape_QC(
  df = df_aptamer,
  batches = list(
    "INTERVAL Batch 1" = "INTERVAL_LOD_batch1",
    "INTERVAL Batch 2" = "INTERVAL_LOD_batch2",
    "CHRIS"            = "CHRIS_LOD"
  ),
  transform_fun = transform_LOD
)

p1_d <- make_qc_violin(long_LOD, ylab = "LOD (log10 RFU)")



  ## E – CV (%) 0–50
long_CV <- reshape_QC(
  df = df_aptamer,
  batches = list(
    "INTERVAL Batch 1" = "INTERVAL_CV_batch1",
    "INTERVAL Batch 2" = "INTERVAL_CV_batch2",
    "CHRIS"            = "CHRIS_CV"
  ),
  transform_fun = transform_CV
)

p1_e <- make_qc_violin(long_CV, ylab = "CV (% - log scale)", y_scale="log10")



## F – % aptamers with >20% samples below LOD (per batch)

long_pctLOD <- bind_rows(
  tibble(
    Batch = "INTERVAL Batch 1",
    SeqID = df_aptamer$SeqID,
    group_lbl = df_aptamer$group_lbl,
    value = df_aptamer$INTERVAL_Percent_under_LOD_batch1
  ),
  tibble(
    Batch = "INTERVAL Batch 2",
    SeqID = df_aptamer$SeqID,
    group_lbl = df_aptamer$group_lbl,
    value = df_aptamer$INTERVAL_Percent_under_LOD_batch2
  ),
  tibble(
    Batch = "CHRIS",
    SeqID = df_aptamer$SeqID,
    group_lbl = df_aptamer$group_lbl,
    value = df_aptamer$CHRIS_Percent_under_LOD
  )
) %>%
  # filter(!is.na(value)) %>%            # critical: must drop NA batch entries
   mutate(value = ifelse(is.na(value),0,value)) %>%    
  mutate(high = value > 20)
#problem: not all values, the NA are

summary_pct <- long_pctLOD %>%
  group_by(Batch, group_lbl) %>%
  summarise(
    total = n_distinct(SeqID),
    n_hi  = n_distinct(SeqID[high]),
    pct   = 100 * n_hi / total,
    label = paste0(n_hi, "/", total, "\n(", sprintf("%.1f", pct), "%)"),
    .groups = "drop"
  )

ymax_f <- max(summary_pct$pct) * 1.25

p1_f <- ggplot(summary_pct,
               aes(x = group_lbl, y = pct, fill = Batch)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.55) +
  geom_text(
    aes(label = label),
    vjust = -0.2,
    position = position_dodge(width = 0.7),
    size = 3
  ) +
  scale_y_continuous(limits = c(0, ymax_f),
                     expand = expansion(mult = c(0,0.12))) +
  scale_fill_manual(values = batch_colors) +
  labs(x = NULL, y = "aptamers with >20% below LOD") +
  theme_bw(14) + 
  theme(legend.position = "bottom", legend.title    = element_blank())

#######################
##poisson number of samples below LOD
df_aptamer<-df_aptamer %>%
mutate(INTERVAL_Number_under_LOD_batch1 = ifelse(is.na(INTERVAL_Number_under_LOD_batch1),0,INTERVAL_Number_under_LOD_batch1),
       INTERVAL_Number_under_LOD_batch2 = ifelse(is.na(INTERVAL_Number_under_LOD_batch2),0,INTERVAL_Number_under_LOD_batch2),
       CHRIS_Number_under_LOD = ifelse(is.na(CHRIS_Number_under_LOD),0,CHRIS_Number_under_LOD)                                          )


mod <- glm(INTERVAL_Number_under_LOD_batch1 ~ group_lbl, family="poisson", data=df_aptamer)
summary(mod)
mod <- glm(INTERVAL_Number_under_LOD_batch2 ~ group_lbl, family="poisson", data=df_aptamer)
summary(mod)
mod <- glm(CHRIS_Number_under_LOD ~ group_lbl, family="poisson", data=df_aptamer)
summary(mod)

  ##########################################################
  ## 3. Legend for QC panels
  ##########################################################

  legend_df <- data.frame(
    Batch = factor(
      c("INTERVAL Batch 1", "INTERVAL Batch 2", "CHRIS"),
      levels = c("INTERVAL Batch 1", "INTERVAL Batch 2", "CHRIS")
    )
  )

  legend_qc <- cowplot::get_legend(
    ggplot(legend_df, aes(x = 1, y = 1, fill = Batch)) +
      geom_point(size = 4, shape = 21) +
      scale_fill_manual(values = batch_colors, name = NULL) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "bottom",
        legend.title    = element_blank(),
        axis.text       = element_blank(),
        axis.title      = element_blank(),
        axis.ticks      = element_blank(),
        panel.grid      = element_blank()
      )
  )

  ##########################################################
  ## 4. Assemble Figure 1
  ##########################################################

  row1 <- plot_grid(p1_a, p1_d, ncol = 2, rel_widths = c(1.4, 1),
                    labels = c("a.", "d."), label_size = 14,
                    label_fontface = "bold")

  row2 <- plot_grid(p1_b, p1_e, ncol = 2, rel_widths = c(1.4, 1),
                    labels = c("b.", "e."), label_size = 14,
                    label_fontface = "bold")

  row3 <- plot_grid(p1_c, p1_f, ncol = 2, rel_widths = c(1.4, 1),
                    labels = c("c.", "f."), label_size = 14,
                    label_fontface = "bold")

  legend_row <- plot_grid(legend_qc, ncol = 1)

  figure1 <- plot_grid(
    row1,
    row2,
    row3,
    plot_grid(p1_g, labels = "g.", label_size = 14, label_fontface = "bold"),
    legend_row,
    ncol = 1,
    rel_heights = c(1.05, 1.05, 1.05, 1.0, 0.25)
  )

  ggsave("Figure_1_FINAL.png", figure1, width = 13, height = 15, dpi = 300)
  cat("Figure 1 saved as Figure_1_FINAL.png\n")

} # make_figure1 = TRUE




############################################################
## FIGURE 2 — CIS-pQTL STRATIFIED PANELS: here we need to average across QC values (not by cohort)
############################################################
# LOD_log_mean (mean log-LOD across batches)
# CV_mean_pct (mean CV ×100 across batches)
# pct_underLOD_mean (mean %<LOD across batches)

############################################################
## A. Unified aptamer universe (Table 2 denominator)
############################################################
# pass_qc is per-row in Table 2 → collapse to SeqID
# Collapse Supplementary Table 2 to one row per SeqID
## SeqIDs may target multiple proteins; a representative UniProt ID is selected
## for bookkeeping and plotting purposes only.

apt_passqc <- annot %>%
  group_by(SeqID) %>%
  summarise(
    pass_qc     = any(pass_qc),
    UniProt_ID  = representative_uniprot(UniProt_ID),
    HGNC_Symbol = representative_uniprot(HGNC_Symbol),
    aptamer_new = any(aptamer_new),
    .groups = "drop"
  ) %>%
  filter(pass_qc) %>%
  mutate(
    group_lbl = factor(
      ifelse(aptamer_new, "New in 7k", "Present in 5k"),
      levels = c("Present in 5k", "New in 7k")
    )
  )

############################################################
## B. pQTL membership flags from Table 3
############################################################

apt_has_any <- pqtl_raw %>%
  distinct(SeqID) %>%
  mutate(has_any_pqtl = TRUE)

apt_has_cis <- pqtl_raw %>%
  filter(cis_or_trans == "cis") %>%
  distinct(SeqID) %>%
  mutate(has_cis_pqtl = TRUE)

############################################################
## C. QC summarisation (mean across batches)
############################################################

df_qc_summary <- df_aptamer %>%
  semi_join(apt_passqc, by = "SeqID") %>%
  mutate(
    LOD_log_mean = rowMeans(cbind(
      transform_LOD(INTERVAL_LOD_batch1),
      transform_LOD(INTERVAL_LOD_batch2),
      transform_LOD(CHRIS_LOD)
    ), na.rm = TRUE),

    CV_mean_pct = rowMeans(cbind(
      transform_CV(INTERVAL_CV_batch1),
      transform_CV(INTERVAL_CV_batch2),
      transform_CV(CHRIS_CV)
    ), na.rm = TRUE),

    pct_underLOD_mean = rowMeans(cbind(
      INTERVAL_Percent_under_LOD_batch1,
      INTERVAL_Percent_under_LOD_batch2,
      CHRIS_Percent_under_LOD
    ), na.rm = TRUE)
  ) %>%
  select(
    SeqID, UniProt_ID, HGNC_Symbol, aptamer_new, group_lbl,
    LOD_log_mean, CV_mean_pct, pct_underLOD_mean
  )

############################################################
## D. Final master table (THIS is the alignment point)
############################################################

apt_master <- df_qc_summary %>%
  left_join(apt_has_any, by = "SeqID") %>%
  left_join(apt_has_cis, by = "SeqID") %>%
  mutate(
    has_any_pqtl = replace_na(has_any_pqtl, FALSE),
    has_cis_pqtl = replace_na(has_cis_pqtl, FALSE),
    pqtl_group = factor(
      ifelse(has_cis_pqtl, "Has ≥1 cis-pQTL", "No cis-pQTL"),
      levels = c("No cis-pQTL", "Has ≥1 cis-pQTL")
    )
  )

# sanity checks
table(apt_master$has_any_pqtl)
table(apt_master$has_cis_pqtl)

# UniProt-level cis status based on QC-passing aptamers only
cis_protein_qc <- apt_master %>%
  filter(has_cis_pqtl) %>%
  distinct(UniProt_ID) %>%
  mutate(pqtl_group = "Has ≥1 cis-pQTL")

df_uniprot2 <- df_uniprot1 %>%
  left_join(cis_protein_qc, by = "UniProt_ID") %>%
  mutate(
    pqtl_group = ifelse(
      is.na(pqtl_group),
      "No cis-pQTL",
      pqtl_group
    ),
    pqtl_group = factor(
      pqtl_group,
      levels = c("No cis-pQTL", "Has ≥1 cis-pQTL")
    )
  )

# df_uniprot2 has: UniProt_ID, HGNC_Symbol, features, Blood_conc_pgL,
# group_lbl (Present in 5k / New in 7k), pqtl_group (No cis-pQTL / Has ≥1 cis-pQTL)

############################################################
## 5.1 Stacked bar helper for Figure 2 (A–C), 100% stacked
##     Facetted by pqtl_group, same style as Figure 1
############################################################

make_panel2 <- function(df, feature, y_max, valid_levels) {

  df_filtered <- df %>%
    group_by(UniProt_ID, HGNC_Symbol, group_lbl, pqtl_group) %>%
    summarise(
      feature_val = .data[[feature]],
      .groups = "drop"
    ) %>%
    filter(!is.na(feature_val), feature_val %in% valid_levels) %>%
    mutate(feature_val = factor(feature_val, levels = valid_levels))

  if (nrow(df_filtered) == 0) {
    empty <- df %>% distinct(pqtl_group, group_lbl)
    return(
      ggplot(empty, aes(x = group_lbl, y = 0)) +
        geom_blank() +
        facet_wrap(~ pqtl_group, nrow = 1) +
        labs(x = NULL, y = "Percentage of proteins") +
        theme_bw(14)
    )
  }

  panel_data <- df_filtered %>%
    count(pqtl_group, group_lbl, feature_val) %>%
    group_by(pqtl_group, group_lbl) %>%
    mutate(
      total   = sum(n),
      percent = 100 * n / total,
      label   = paste0(n, " (", round(percent, 1), "%)")
    )

  summary_counts <- df_filtered %>%
    group_by(pqtl_group, group_lbl) %>%
    summarise(
      proteins_with_data = n_distinct(UniProt_ID),
      proteins_total = n_distinct(df$UniProt_ID[df$group_lbl == unique(group_lbl) &
                                 df$pqtl_group == unique(pqtl_group)]),
      proteins_missing = proteins_total - proteins_with_data,
      genes    = n_distinct(HGNC_Symbol),
      label    = paste0(proteins_with_data, " proteins\n(", proteins_missing, " missing)"),
      .groups  = "drop"
    )

  ggplot(panel_data,
         aes(x = group_lbl, y = percent, fill = feature_val)) +
    geom_col(position = "stack", alpha = 0.9) +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 3.3, fontface = "bold", colour = "white"
    ) +
    geom_text(
      data = summary_counts,
      aes(x = group_lbl, y = y_max, label = label),
      vjust = 0,
      size  = 3,
      fontface = "italic",
      inherit.aes = FALSE
    ) +
    facet_wrap(~ pqtl_group, nrow = 1) +
    scale_y_continuous(
        limits = c(0, y_max * 1.02),
        breaks = c(0, 25, 50, 75, 100),
        expand = expansion(mult = c(0, 0.15))) +
    scale_fill_brewer(palette = "Set2", name = NULL) +
    labs(x = NULL, y = "Percentage of proteins") +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom")
}

############################################################
## 5.2 Blood concentration panel (g), cis-stratified
############################################################

df_blood2 <- df_uniprot2 %>%
  filter(!is.na(Blood_conc_pgL), Blood_conc_pgL > 0)

summary_counts_g2 <- df_uniprot2 %>%
  group_by(pqtl_group, group_lbl) %>%
  summarise(
    proteins_total   = n_distinct(UniProt_ID),
    genes            = n_distinct(HGNC_Symbol),
    proteins_nonmiss = n_distinct(
      UniProt_ID[!is.na(Blood_conc_pgL) & Blood_conc_pgL > 0]
    ),
    .groups = "drop"
  ) %>%
  mutate(
    proteins_missing = proteins_total - proteins_nonmiss,
    label = paste0(
       proteins_nonmiss, " proteins\n(",
       proteins_missing, " missing)"
    )
  )
#      "N=", proteins_nonmiss, " proteins\n(",
#      genes, " genes; ", proteins_missing, " missing)"

p2_g <- ggplot(df_blood2,
               aes(x = group_lbl, y = Blood_conc_pgL, fill = group_lbl)) +
  geom_violin(alpha = 0.4, scale = "area") +
  geom_boxplot(width = 0.15, outlier.size = 0.4) +
  geom_text(
    data = summary_counts_g2,
    aes(x = group_lbl, y = Inf, label = label),
    vjust = 1.2,
    size  = 3,
    fontface = "italic",
    inherit.aes = FALSE
  ) +
  facet_wrap(~ pqtl_group, nrow = 1) +
  scale_y_log10(
    expand = expansion(mult = c(0.02, 0.20))
  ) +
  scale_fill_manual(values = protein_panel_colors, name = NULL) +
  labs(x = NULL, y = "Blood conc. (log10 pg/L)") +
  theme_bw(14) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(margin = margin(r = 30))
  )

############################################################
## 5.3 QC reshaping with cis-stratification (for d/e/f)
############################################################

reshape_QC_strat <- function(df, qc_col) {
  df %>%
    select(SeqID, UniProt_ID, group_lbl, pqtl_group, {{ qc_col }}) %>%
    rename(value = {{ qc_col }}) %>%
    filter(!is.na(value))
}

make_qc_violin_strat <- function(long_df, ylab, y_limits = NULL) {

  counts <- long_df %>%
    group_by(pqtl_group, group_lbl) %>%
    summarise(
      aptamers = n_distinct(SeqID),
      proteins = n_distinct(UniProt_ID),
      label = paste0(aptamers, " aptamers\n", proteins, " proteins"),
      .groups = "drop"
    )

  ggplot(long_df, aes(x = group_lbl, y = value, fill = group_lbl)) +
    geom_violin(alpha = 0.5, scale = "area") +
    geom_boxplot(width = 0.12, outlier.size = 0.4) +
    geom_text(
      data = counts,
      aes(x = group_lbl, y = Inf, label = label),
      vjust = 1.1, size = 3, fontface = "italic",
      inherit.aes = FALSE
    ) +
    facet_wrap(~ pqtl_group, nrow = 1) +
    scale_y_continuous(limits = y_limits,
                       expand = expansion(mult = c(0.02, 0.20))) +
    scale_fill_manual(values = protein_panel_colors, name = NULL) +
    labs(x = NULL, y = ylab) +
    theme_bw(14) +
    theme(legend.position = "none")
}

############################################################
## 5. Build Figure 2
############################################################
if (make_figure2) {

  p2_a <- make_panel2(df_uniprot2,
                      "Secretion.pathway", 110,
                      c("intracellular", "membrane", "secreted"))

  p2_b <- make_panel2(df_uniprot2,
                      "RNA.tissue.distribution", 105,
                      c("Detected in all","Detected in many",
                        "Detected in some","Detected in single"))

  p2_c <- make_panel2(df_uniprot2,
                      "RNA.tissue.specificity", 105,
                      c("Group enriched","Low tissue specificity",
                        "Tissue enhanced","Tissue enriched"))

long_LOD2 <- reshape_QC_strat(apt_master, LOD_log_mean)
p2_d <- make_qc_violin_strat(long_LOD2, "LOD (log10 RFU)")

long_CV2 <- reshape_QC_strat(apt_master, CV_mean_pct)
p2_e <- make_qc_violin_strat(long_CV2, "CV (%)", y_limits = c(0, 50))

p2_f <- apt_master %>%
  filter(!is.na(pct_underLOD_mean)) %>%   # critical fix
  mutate(high = pct_underLOD_mean > 20) %>%
  group_by(pqtl_group, group_lbl) %>%
  summarise(
    total = n_distinct(SeqID),
    n_hi  = n_distinct(SeqID[high]),
    pct   = 100 * n_hi / total,
    label = paste0(
      n_hi, "/", total, "\n(",
      sprintf("%.1f", pct), "%)"
    ),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = group_lbl, y = pct, fill = group_lbl)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = label), vjust = -0.1, size = 3) +
  facet_wrap(~ pqtl_group, nrow = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = NULL, y = "apts with >20% below LOD") +
  theme_bw(14) +
  theme(legend.position = "none")

  row1_f2 <- plot_grid(p2_a, p2_d, ncol = 2,
                       labels = c("a.", "d."),
                       label_size = 14, label_fontface = "bold")

  row2_f2 <- plot_grid(p2_b, p2_e, ncol = 2,
                       labels = c("b.", "e."),
                       label_size = 14, label_fontface = "bold")

  row3_f2 <- plot_grid(p2_c, p2_f, ncol = 2,
                       labels = c("c.", "f."),
                       label_size = 14, label_fontface = "bold")

  figure2 <- plot_grid(row1_f2, row2_f2, row3_f2,
    		       plot_grid(p2_g, labels = "g.", label_size = 14, label_fontface = "bold"),
                       ncol = 1,
                       rel_heights = c(1.05, 1.05, 1.05, 1, 0.25))

  ggsave("Figure_2_FINAL.png", figure2,
         width = 13, height = 15, dpi = 300)

  cat("Figure 2 saved as Figure_2_FINAL.png\n")
}



############################################################
############## MODELS — USING CORRECTED QC METRICS #########
############################################################

if (run_models) {

  suppressPackageStartupMessages({
    library(dplyr)
    library(broom)
    library(purrr)
  })

  message("Running models with corrected QC universe...")

  ##########################################################
  ## Modeling dataset (aligned with Figure 2)
  ##########################################################

  apt_for_model <- apt_master %>%
  left_join(
    df_uniprot1 %>%
      select(UniProt_ID, Secretion.pathway, Blood_conc_pgL),
    by = "UniProt_ID"
  ) %>%
  mutate(
    secretion_class = factor(
      Secretion.pathway,
      levels = c("intracellular", "membrane", "secreted")
    ),
    log10_blood = log10(Blood_conc_pgL + 1)
  ) %>%
  filter(
    !is.na(pct_underLOD_mean),
    !is.na(CV_mean_pct)
  )

  # sanity checks (should show both TRUE/FALSE)
  # Modeling universe:
  # QC-passing aptamers with non-missing mean QC metrics
  nrow(apt_for_model)
  table(apt_for_model$has_any_pqtl)
  table(apt_for_model$has_cis_pqtl)
  table(apt_for_model$aptamer_new)

  ##########################################################
  ## UNADJUSTED MODELS
  ##########################################################

  ## % below LOD
  m_cis_pctLOD <- glm(
    has_cis_pqtl ~ pct_underLOD_mean,
    data = apt_for_model,
    family = binomial
  )

  m_any_pctLOD <- glm(
    has_any_pqtl ~ pct_underLOD_mean,
    data = apt_for_model,
    family = binomial
  )

  ## mean LOD
  m_cis_LOD <- glm(
    has_cis_pqtl ~ LOD_log_mean,
    data = apt_for_model,
    family = binomial
  )

  m_any_LOD <- glm(
    has_any_pqtl ~ LOD_log_mean,
    data = apt_for_model,
    family = binomial
  )


## CV only
m_cis_CV <- glm(
  has_cis_pqtl ~ CV_mean_pct,
  data = apt_for_model,
  family = binomial
)

m_any_CV <- glm(
  has_any_pqtl ~ CV_mean_pct,
  data = apt_for_model,
  family = binomial
)

## Secretion class only
m_cis_secretion <- glm(
  has_cis_pqtl ~ secretion_class,
  data = apt_for_model,
  family = binomial
)

m_any_secretion <- glm(
  has_any_pqtl ~ secretion_class,
  data = apt_for_model,
  family = binomial
)

## Plasma concentration only
m_cis_blood <- glm(
  has_cis_pqtl ~ log10_blood,
  data = apt_for_model,
  family = binomial
)

m_any_blood <- glm(
  has_any_pqtl ~ log10_blood,
  data = apt_for_model,
  family = binomial
)

## Panel status
m_cis_new <- glm(
  has_cis_pqtl ~ aptamer_new,
  data = apt_for_model,
  family = binomial
)

m_any_new <- glm(
  has_any_pqtl ~ aptamer_new,
  data = apt_for_model,
  family = binomial
)


  ##########################################################
  ## ADJUSTED MODELS (biology-aware)
  ##########################################################

  ## LOD-adjusted
  m_cis_LOD_adj <- glm(
    has_cis_pqtl ~ aptamer_new +
      LOD_log_mean +
      secretion_class +
      log10_blood,
    data = apt_for_model,
    family = binomial
  )

  m_any_LOD_adj <- glm(
    has_any_pqtl ~ aptamer_new +
      LOD_log_mean +
      secretion_class +
      log10_blood,
    data = apt_for_model,
    family = binomial
  )

  ## % below LOD-adjusted
m_cis_pctLOD_adj <- glm(
  has_cis_pqtl ~ aptamer_new +
    pct_underLOD_mean +
    secretion_class +
    log10_blood,
  data = apt_for_model,
  family = binomial
)

m_any_pctLOD_adj <- glm(
  has_any_pqtl ~ aptamer_new +
    pct_underLOD_mean +
    secretion_class +
    log10_blood,
  data = apt_for_model,
  family = binomial
)


  ##########################################################
  ## Collect results
  ##########################################################

all_models <- list(
  cis_vs_new          = m_cis_new,
  any_vs_new          = m_any_new,
  cis_vs_pctLOD       = m_cis_pctLOD,
  any_vs_pctLOD       = m_any_pctLOD,
  cis_vs_CV           = m_cis_CV,
  any_vs_CV           = m_any_CV,
  cis_vs_secretion    = m_cis_secretion,
  any_vs_secretion    = m_any_secretion,
  cis_vs_blood        = m_cis_blood,
  any_vs_blood        = m_any_blood,
  cis_pctLOD_adj      = m_cis_pctLOD_adj,
  any_pctLOD_adj      = m_any_pctLOD_adj
)

  supp_models <- purrr::imap_dfr(all_models, function(mod, name) {
    broom::tidy(mod) %>%
      mutate(
        model = name,
        OR = exp(estimate)
      )
  }) %>%
    select(model, term, estimate, std.error, statistic, p.value, OR)

## after you create `supp_models` in run_models:
supp_models <- supp_models %>%
  mutate(
    module = "cis_discovery",
    outcome = ifelse(grepl("^cis_", model), "has_cis_pqtl", "has_any_pqtl"),
    analysis_unit = "SeqID (all QC-passing aptamers)",
    n = stats::nobs(all_models[[1]])  # optional; or compute per-model below
  )

  write.csv(
    supp_models,
    models_file,
    row.names = FALSE
  )

  message("Models completed:", models_file)
}

  ##########################################################
  ## Complexity models
  ##########################################################
if (add_complexity_models) {

  message("Running complexity module using Supplementary Table 4 (ST4)...")

  st4 <- readxl::read_excel("supplementary_table_4.xlsx", sheet = "ST4", skip = 2)
  names(st4) <- make.names(names(st4))   # <-- critical

  message("ST4 columns:")
  print(names(st4))

  # Keep only cis (already defined by your pipeline: ±500kb around TSS)
  st4_cis <- st4 %>%
    dplyr::filter(cis_or_trans == "cis")

  message("Rows in ST4 (cis only): ", nrow(st4_cis))

# For each SeqID, count how many independent cis variants exist (i.e., how many rows in ST4_cis)
  cis_complexity_by_seqid <- st4_cis %>%
    dplyr::group_by(SeqID) %>%
    dplyr::summarise(
      UniProt_ID = dplyr::first(UniProt_ID),
      uniprot_new_in_somascan7k_vs5k = dplyr::first(uniprot_new_in_somascan7k_vs5k),
      n_indep_cis = dplyr::n(),                    # number of conditionally independent cis variants
      has_secondary_cis = (n_indep_cis > 1),       # your “secondary signal” definition
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      group_lbl = factor(
        ifelse(uniprot_new_in_somascan7k_vs5k, "New in 7k", "Present in 5k"),
        levels = c("Present in 5k", "New in 7k")
      )
    )

  # Quick summary table (proportions)
  cis_complexity_summary <- cis_complexity_by_seqid %>%
    dplyr::group_by(group_lbl) %>%
    dplyr::summarise(
      n_with_cis = dplyr::n(),
      n_secondary = sum(has_secondary_cis),
      pct_secondary = 100 * mean(has_secondary_cis),
      .groups = "drop"
    )

  print(cis_complexity_summary)

  p_complexity_bar <- cis_complexity_by_seqid %>%
    dplyr::count(group_lbl, has_secondary_cis) %>%
    dplyr::group_by(group_lbl) %>%
    dplyr::mutate(pct = 100 * n / sum(n)) %>%
    ggplot2::ggplot(ggplot2::aes(x = group_lbl, y = pct, fill = has_secondary_cis)) +
    ggplot2::geom_col() +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(
      x = NULL,
      y = "Percentage of proteins with ≥1 cis-pQTL",
      fill = "Has ≥2 independent\ncis variants"
    ) +
    ggplot2::theme_bw(14)

  print(p_complexity_bar)

m_secondary_new <- glm(
    has_secondary_cis ~ uniprot_new_in_somascan7k_vs5k,
    data = cis_complexity_by_seqid,
    family = binomial
  )
  print(summary(m_secondary_new))
  message("OR (new vs old) = ", exp(coef(m_secondary_new)[2]))


## Complexity model results in same schema
complexity_models <- list(
  cis_secondary_vs_new = m_secondary_new
  # add more here later, e.g. adjusted models
)

supp_complexity <- purrr::imap_dfr(complexity_models, function(mod, name) {
  broom::tidy(mod) %>%
    mutate(
      model = name,
      OR = exp(estimate),
      module = "cis_complexity",
      outcome = "has_secondary_cis",
      analysis_unit = "SeqID (targets with ≥1 cis-pQTL)",
      n = stats::nobs(mod)
    )
}) %>%
  select(module, outcome, analysis_unit, n, model, term, estimate, std.error, statistic, p.value, OR)

## Append to main table object (if it exists)
if (exists("supp_models")) {
  supp_models <- bind_rows(
    supp_models %>% select(module, outcome, analysis_unit, n, model, term, estimate, std.error, statistic, p.value, OR),
    supp_complexity
  )
} else {
  supp_models <- supp_complexity
}

## Write a combined table (new filename so you don’t silently overwrite old meaning)
write.csv(
  supp_models,
  models_file,
  row.names = FALSE
)

  message("Models completed:", models_file)

}


