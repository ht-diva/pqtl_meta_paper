## =========================================================
## 0. Libraries, helpers, paths
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(tidyr)
  library(cowplot)
  library(RColorBrewer)
})

data_path <- "/Users/c.giambartolomei/Downloads/annotations_of_proteins_oct222025/manuscript_tables/"
setwd(data_path)

## Mode function for categorical annotations
mode_function <- function(x) {
  ux <- na.omit(x)
  if (length(ux) == 0) return(NA_character_)
  names(sort(table(ux), decreasing = TRUE))[1]
}

## =========================================================
## 1. Load data (raw, unfiltered tables)
## =========================================================

table_annot <- read_excel("supplementary_table_2.xlsx", sheet = 3)
table_pqtl  <- read_excel("supplementary_table_3.xlsx", sheet = "ST3", skip = 1)

names(table_annot) <- make.names(names(table_annot))
names(table_pqtl)  <- make.names(names(table_pqtl))

## Harmonise SeqID name
table_annot <- table_annot %>% rename(SeqID = SeqId)

## Aptamer-level "new vs old" flag
table_annot <- table_annot %>%
  mutate(aptamer_new_in_somascan7k_vs5k = !aptamer_in_somascan5k)

## Save list of 5k SeqIDs
seqids_5k <- table_annot %>%
  filter(aptamer_in_somascan5k) %>%
  distinct(SeqID)

cat("=========================================================\n")
cat("BASIC ANNOTATION TABLE COUNTS (RAW TABLE_ANNOT)\n")
cat("=========================================================\n")
cat("Rows in table_annot (aptamer–UniProt–gene combos):", nrow(table_annot), "\n")
cat("Unique aptamers (SeqID):", n_distinct(table_annot$SeqID), "\n")
cat("Unique UniProt_ID:",      n_distinct(table_annot$UniProt_ID), "\n")
cat("Unique HGNC_Symbol:",     n_distinct(table_annot$HGNC_Symbol), "\n\n")

cat("SeqIDs present in SomaScan 5k:", nrow(seqids_5k), "\n\n")

## How many aptamers map to >1 UniProt / >1 Entrez (within this table)?
map_multi_uniprot <- table_annot %>%
  group_by(SeqID) %>%
  summarise(n_uniprot = n_distinct(UniProt_ID), .groups = "drop")

map_multi_entrez <- table_annot %>%
  group_by(SeqID) %>%
  summarise(n_entrez = n_distinct(Entrez_Gene_ID), .groups = "drop")

cat("Aptamers mapping to >1 UniProt_ID:", sum(map_multi_uniprot$n_uniprot > 1), "\n")
cat("Aptamers mapping to >1 Entrez_Gene_ID:", sum(map_multi_entrez$n_entrez > 1), "\n\n")

cat("uniprot_new_in_somascan7k_vs5k (row-level):\n")
print(table(table_annot$uniprot_new_in_somascan7k_vs5k, useNA = "ifany"))
cat("\n")

## =========================================================
## 2. FIGURE 1 – v5 vs v7 (NO pQTL stratification)
##    LEFT: protein-level (UniProt_ID)
##    RIGHT: aptamer-level QC (SeqID)
##    Data: table_annot only (raw)
## =========================================================

## 2.1 Collapse to unique UniProt_ID, keeping dominant annotations (protein-level)
df_uniprot1 <- table_annot %>%
  group_by(UniProt_ID) %>%
  summarise(
    HGNC_Symbol             = mode_function(HGNC_Symbol),
    Secretion.pathway       = mode_function(Secretion.pathway),
    RNA.tissue.distribution = mode_function(RNA.tissue.distribution),
    RNA.tissue.specificity  = mode_function(RNA.tissue.specificity),
    Blood_conc_pgL          = median(as.numeric(Blood_conc_pgL), na.rm = TRUE),
    uniprot_new_in_somascan7k_vs5k = any(uniprot_new_in_somascan7k_vs5k, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    group_lbl = factor(
      ifelse(uniprot_new_in_somascan7k_vs5k, "New in 7k", "Present in 5k"),
      levels = c("Present in 5k", "New in 7k")
    )
  )

cat("=========================================================\n")
cat("FIGURE 1 – PROTEIN-LEVEL GROUP COUNTS (RAW)\n")
cat("=========================================================\n")
prot_counts_f1 <- df_uniprot1 %>%
  group_by(group_lbl) %>%
  summarise(
    n_proteins = n_distinct(UniProt_ID),
    n_genes    = n_distinct(HGNC_Symbol),
    .groups    = "drop"
  )
print(prot_counts_f1)
cat("\n")

## 2.2 Generic stacked bar panel for A–C (no pQTL facets here)
make_panel1 <- function(df, feature, y_max, valid_levels) {

  df_filtered <- df %>%
    group_by(UniProt_ID, HGNC_Symbol, group_lbl) %>%
    summarise(
      feature_val = mode_function(.data[[feature]]),
      .groups = "drop"
    ) %>%
    filter(!is.na(feature_val), feature_val %in% valid_levels)

  df_filtered <- df_filtered %>%
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
      proteins = n_distinct(UniProt_ID),
      genes    = n_distinct(HGNC_Symbol),
      label    = paste0("N=", proteins, " proteins\n(", genes, " genes)"),
      .groups  = "drop"
    )

  ggplot(panel_data, aes(x = group_lbl, y = percent, fill = feature_val)) +
    geom_col(position = "stack", alpha = 0.9) +
    geom_text(aes(label = label),
              position = position_stack(vjust = 0.5),
              size = 3.5, fontface = "bold", color = "white") +
    geom_text(data = summary_counts,
              aes(x = group_lbl, y = y_max, label = label),
              vjust = 0, inherit.aes = FALSE,
              size = 3, fontface = "italic") +
    scale_y_continuous(limits = c(0, y_max + 10)) +
    scale_fill_brewer(palette = "Set2", name = NULL) +
    labs(x = NULL, y = "Percentage of proteins") +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom")
}

## 2.3 Left panels A–C

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
  valid_levels = c("Detected in all",
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

## 2.4 Left panel G: Blood concentration (log10), non-stratified

p1_g_blood <- df_uniprot1 %>%
  filter(!is.na(Blood_conc_pgL), Blood_conc_pgL > 0) %>%
  ggplot(aes(x = group_lbl, y = Blood_conc_pgL, fill = group_lbl)) +
  geom_violin(alpha = 0.4, position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.15, position = position_dodge(width = 0.8),
               outlier.size = 0.4) +
  scale_y_log10() +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  labs(x = NULL, y = "Expected concentration [pg/L] (log10)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

## Blood concentration Wilcoxon test (log10) for manuscript
bc_df <- df_uniprot1 %>%
  filter(!is.na(Blood_conc_pgL), Blood_conc_pgL > 0) %>%
  mutate(log_Blood = log10(Blood_conc_pgL))

cat("=========================================================\n")
cat("FIGURE 1 – BLOOD CONCENTRATION STATS (PROTEIN-LEVEL, RAW)\n")
cat("=========================================================\n")
bc_summary <- bc_df %>%
  group_by(group_lbl) %>%
  summarise(
    median   = median(Blood_conc_pgL, na.rm = TRUE),
    q1       = quantile(Blood_conc_pgL, 0.25, na.rm = TRUE),
    q3       = quantile(Blood_conc_pgL, 0.75, na.rm = TRUE),
    median_log = median(log_Blood, na.rm = TRUE),
    .groups  = "drop"
  )
print(bc_summary)

if (n_distinct(bc_df$group_lbl) == 2) {
  w_bc <- wilcox.test(log_Blood ~ group_lbl, data = bc_df, exact = FALSE)
  cat("\nWilcoxon test (log10 Blood_conc_pgL, new vs old):\n")
  cat("  W =", w_bc$statistic, "  p =", signif(w_bc$p.value, 3), "\n\n")
}

## 2.5 Aptamer-level QC (RIGHT panels, non-stratified, unique SeqID)

aptamer_qc1 <- table_annot %>%
  group_by(SeqID) %>%
  summarise(
    UniProt_ID                     = mode_function(UniProt_ID),
    HGNC_Symbol                    = mode_function(HGNC_Symbol),
    uniprot_new_in_somascan7k_vs5k = any(uniprot_new_in_somascan7k_vs5k, na.rm = TRUE),
    aptamer_in_somascan5k          = any(aptamer_in_somascan5k, na.rm = TRUE),
    aptamer_new_in_somascan7k_vs5k = !aptamer_in_somascan5k,
    INTERVAL_LOD_batch1            = mean(INTERVAL_LOD_batch1, na.rm = TRUE),
    INTERVAL_LOD_batch2            = mean(INTERVAL_LOD_batch2, na.rm = TRUE),
    INTERVAL_Percent_under_LOD_batch1 = mean(INTERVAL_Percent_under_LOD_batch1, na.rm = TRUE),
    INTERVAL_Percent_under_LOD_batch2 = mean(INTERVAL_Percent_under_LOD_batch2, na.rm = TRUE),
    INTERVAL_CV_batch1             = mean(INTERVAL_CV_batch1, na.rm = TRUE),
    INTERVAL_CV_batch2             = mean(INTERVAL_CV_batch2, na.rm = TRUE),
    CHRIS_LOD                      = mean(CHRIS_LOD, na.rm = TRUE),
    CHRIS_Percent_under_LOD        = mean(CHRIS_Percent_under_LOD, na.rm = TRUE),
    CHRIS_CV                       = mean(CHRIS_CV, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    group_lbl = factor(
      ifelse(aptamer_new_in_somascan7k_vs5k, "New in 7k", "Present in 5k"),
      levels = c("Present in 5k", "New in 7k")
    )
  )

batch_colors <- c(
  "INTERVAL Batch 1" = "#1f78b4",
  "INTERVAL Batch 2" = "#a6cee3",
  "CHRIS"            = "#ff7f00"
)

## Fix of previous error: y_value now carried into counts as y_pos
make_qc_plot1 <- function(df, value_cols, ylab, log_offset = 0) {
  long_df <- df %>%
    select(SeqID, UniProt_ID, group_lbl, all_of(value_cols)) %>%
    pivot_longer(cols = all_of(value_cols),
                 names_to = "Batch",
                 values_to = "value") %>%
    filter(!is.na(value), value >= 0) %>%
    mutate(
      Batch = recode(
        Batch,
        "INTERVAL_LOD_batch1"                 = "INTERVAL Batch 1",
        "INTERVAL_LOD_batch2"                 = "INTERVAL Batch 2",
        "CHRIS_LOD"                           = "CHRIS",
        "INTERVAL_CV_batch1"                  = "INTERVAL Batch 1",
        "INTERVAL_CV_batch2"                  = "INTERVAL Batch 2",
        "CHRIS_CV"                            = "CHRIS",
        "INTERVAL_Percent_under_LOD_batch1"   = "INTERVAL Batch 1",
        "INTERVAL_Percent_under_LOD_batch2"   = "INTERVAL Batch 2",
        "CHRIS_Percent_under_LOD"             = "CHRIS"
      ),
      y_value = value + log_offset
    )

  if (nrow(long_df) == 0) {
    return(
      ggplot() + geom_blank() +
        labs(x = NULL, y = ylab) +
        theme_bw(14)
    )
  }

  counts <- long_df %>%
    group_by(group_lbl) %>%
    summarise(
      proteins = n_distinct(UniProt_ID),
      aptamers = n_distinct(SeqID),
      label    = paste0("N=", proteins, " proteins\n(", aptamers, " aptamers)"),
      .groups  = "drop"
    ) %>%
    mutate(y_pos = max(long_df$y_value, na.rm = TRUE))

  ggplot(long_df, aes(x = group_lbl, y = y_value, fill = Batch)) +
    geom_violin(scale = "area", alpha = 0.5,
                position = position_dodge(width = 0.8)) +
    geom_boxplot(width = 0.15,
                 position = position_dodge(width = 0.8),
                 outlier.size = 0.4) +
    geom_text(
      data = counts,
      aes(x = group_lbl, y = y_pos, label = label),
      inherit.aes = FALSE,
      size = 3, fontface = "italic", vjust = 1.2
    ) +
    scale_fill_manual(values = batch_colors, name = NULL) +
    labs(x = NULL, y = ylab) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(margin = margin(b = 8))
    ) +
    scale_y_log10()
}

## Panels D–F for Figure 1 (non-stratified)

p1_d <- make_qc_plot1(
  aptamer_qc1,
  value_cols = c("INTERVAL_LOD_batch1", "INTERVAL_LOD_batch2", "CHRIS_LOD"),
  ylab       = "LOD (log10 RFU)",
  log_offset = 0
)

p1_e <- make_qc_plot1(
  aptamer_qc1,
  value_cols = c("INTERVAL_CV_batch1", "INTERVAL_CV_batch2", "CHRIS_CV"),
  ylab       = "CV (%) (log10)",
  log_offset = 0
)

p1_f <- make_qc_plot1(
  aptamer_qc1,
  value_cols = c("INTERVAL_Percent_under_LOD_batch1",
                 "INTERVAL_Percent_under_LOD_batch2",
                 "CHRIS_Percent_under_LOD"),
  ylab       = "% under LOD (log10, +0.01 offset)",
  log_offset = 0.01
)

## 2.6 QC stats + Wilcoxon tests for Figure 1 (aptamer level, raw)

qc_vars <- c(
  "INTERVAL_LOD_batch1", "INTERVAL_LOD_batch2", "CHRIS_LOD",
  "INTERVAL_CV_batch1", "INTERVAL_CV_batch2", "CHRIS_CV",
  "INTERVAL_Percent_under_LOD_batch1",
  "INTERVAL_Percent_under_LOD_batch2",
  "CHRIS_Percent_under_LOD"
)

cat("=========================================================\n")
cat("FIGURE 1 – APTAMER-LEVEL QC STATS (RAW, NON-STRATIFIED)\n")
cat("=========================================================\n")

for (v in qc_vars) {
  cat("\nVariable:", v, "\n")
  df_v <- aptamer_qc1 %>%
    select(group_lbl, !!sym(v)) %>%
    rename(value = !!sym(v)) %>%
    filter(!is.na(value))
  if (nrow(df_v) == 0) {
    cat("  No non-missing values.\n")
  } else {
    summary_v <- df_v %>%
      group_by(group_lbl) %>%
      summarise(
        median = median(value, na.rm = TRUE),
        q1     = quantile(value, 0.25, na.rm = TRUE),
        q3     = quantile(value, 0.75, na.rm = TRUE),
        n      = n(),
        .groups = "drop"
      )
    print(summary_v)
    if (n_distinct(df_v$group_lbl) == 2) {
      w <- wilcox.test(value ~ group_lbl, data = df_v, exact = FALSE)
      cat("  Wilcoxon p =", signif(w$p.value, 3), "\n")
    }
  }
}

## 2.7 Layout & save Figure 1

p1_a        <- p1_a + theme(plot.margin = margin(8, 10, 12, 10))
p1_b        <- p1_b + theme(plot.margin = margin(8, 10, 12, 10))
p1_c        <- p1_c + theme(plot.margin = margin(8, 10, 12, 10))
p1_d        <- p1_d + theme(plot.margin = margin(8, 10, 12, 10))
p1_e        <- p1_e + theme(plot.margin = margin(8, 10, 12, 10))
p1_f        <- p1_f + theme(plot.margin = margin(8, 10, 12, 10))
p1_g_blood  <- p1_g_blood + theme(plot.margin = margin(8, 10, 12, 10))

row1_f1 <- plot_grid(p1_a, p1_d, ncol = 2,
                     rel_widths = c(1.4, 1),
                     labels = c("a.", "d."),
                     label_size = 14,
                     label_fontface = "bold")

row2_f1 <- plot_grid(p1_b, p1_e, ncol = 2,
                     rel_widths = c(1.4, 1),
                     labels = c("b.", "e."),
                     label_size = 14,
                     label_fontface = "bold")

row3_f1 <- plot_grid(p1_c, p1_f, ncol = 2,
                     rel_widths = c(1.4, 1),
                     labels = c("c.", "f."),
                     label_size = 14,
                     label_fontface = "bold")

row4_f1 <- plot_grid(p1_g_blood, NULL, ncol = 2,
                     rel_widths = c(1.4, 1),
                     labels = c("g.", ""),
                     label_size = 14,
                     label_fontface = "bold")

figure1 <- plot_grid(row1_f1, row2_f1, row3_f1, row4_f1,
                     ncol = 1,
                     rel_heights = c(1.05, 1.05, 1.05, 0.8))

ggsave("Figure_1_v5_vs_v7_protein_left_aptamer_QC_right_RAW.png",
       figure1, width = 13, height = 15, dpi = 300)

cat("\nSaved: Figure_1_v5_vs_v7_protein_left_aptamer_QC_right_RAW.png\n\n")

## =========================================================
## 3. FIGURE 2 – v5 vs v7, stratified by cis-pQTL
##     LEFT: protein-level
##     RIGHT: aptamer-level QC, cis vs no cis
## =========================================================

## 3.1 Expand UniProt_IDs in pQTL (handle Q9NPF7|P29460 style)
table_pqtl_expanded <- table_pqtl %>%
  tidyr::separate_rows(UniProt_ID, sep = "\\|") %>%
  mutate(UniProt_ID = trimws(UniProt_ID)) %>%
  filter(!is.na(UniProt_ID), UniProt_ID != "")

cis_by_protein <- table_pqtl_expanded %>%
  filter(cis_or_trans == "cis") %>%
  distinct(UniProt_ID) %>%
  mutate(has_cis_pqtl = TRUE)

any_by_protein <- table_pqtl_expanded %>%
  distinct(UniProt_ID) %>%
  mutate(has_any_pqtl = TRUE)

pqtl_status_protein <- full_join(any_by_protein, cis_by_protein, by = "UniProt_ID") %>%
  mutate(
    has_any_pqtl = ifelse(is.na(has_any_pqtl), FALSE, has_any_pqtl),
    has_cis_pqtl = ifelse(is.na(has_cis_pqtl), FALSE, has_cis_pqtl)
  )

## Aptamer-level cis-pQTL
cis_by_aptamer <- table_pqtl %>%
  filter(cis_or_trans == "cis") %>%
  distinct(SeqID) %>%
  mutate(has_cis_pqtl = TRUE)

cat("=========================================================\n")
cat("FIGURE 2 – BASIC pQTL COUNTS (RAW)\n")
cat("=========================================================\n")
n_total_pqtls <- nrow(table_pqtl)
n_cis         <- sum(table_pqtl$cis_or_trans == "cis", na.rm = TRUE)
n_trans       <- sum(table_pqtl$cis_or_trans == "trans", na.rm = TRUE)

cat("Total aptamer-locus pairs (pQTLs):", n_total_pqtls, "\n")
cat("  cis:",   n_cis, "\n")
cat("  trans:", n_trans, "\n\n")

## 3.2 Protein-level df for Fig 2 (LEFT)
df_uniprot2 <- table_annot %>%
  group_by(UniProt_ID) %>%
  summarise(
    HGNC_Symbol             = mode_function(HGNC_Symbol),
    Secretion.pathway       = mode_function(Secretion.pathway),
    RNA.tissue.distribution = mode_function(RNA.tissue.distribution),
    RNA.tissue.specificity  = mode_function(RNA.tissue.specificity),
    Blood_conc_pgL          = median(as.numeric(Blood_conc_pgL), na.rm = TRUE),
    uniprot_new_in_somascan7k_vs5k =
      any(uniprot_new_in_somascan7k_vs5k, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(pqtl_status_protein, by = "UniProt_ID") %>%
  mutate(
    has_any_pqtl = ifelse(is.na(has_any_pqtl), FALSE, has_any_pqtl),
    has_cis_pqtl = ifelse(is.na(has_cis_pqtl), FALSE, has_cis_pqtl),
    pqtl_group   = factor(
      ifelse(has_cis_pqtl, "Has ≥1 cis-pQTL", "No cis-pQTL"),
      levels = c("No cis-pQTL", "Has ≥1 cis-pQTL")
    ),
    group_lbl    = factor(
      ifelse(uniprot_new_in_somascan7k_vs5k, "New in 7k", "Present in 5k"),
      levels = c("Present in 5k", "New in 7k")
    )
  )

cat("Protein-level cis-pQTL counts by new/old (Fig 2):\n")
prot_counts_f2 <- df_uniprot2 %>%
  group_by(pqtl_group, group_lbl) %>%
  summarise(
    n_proteins = n_distinct(UniProt_ID),
    n_genes    = n_distinct(HGNC_Symbol),
    .groups    = "drop"
  )
print(prot_counts_f2)
cat("\n")

## Fisher tests for cis and any pQTL (protein level)
protein_summary_f2 <- df_uniprot2 %>%
  distinct(UniProt_ID, uniprot_new_in_somascan7k_vs5k, has_cis_pqtl, has_any_pqtl)

f_any <- fisher.test(
  table(protein_summary_f2$uniprot_new_in_somascan7k_vs5k,
        protein_summary_f2$has_any_pqtl)
)
f_cis <- fisher.test(
  table(protein_summary_f2$uniprot_new_in_somascan7k_vs5k,
        protein_summary_f2$has_cis_pqtl)
)

cat("Fisher test (ANY pQTL vs new/old): p =", signif(f_any$p.value, 3), "\n")
cat("Fisher test (CIS pQTL vs new/old): p =", signif(f_cis$p.value, 3), "\n\n")

## 3.3 Left stacked panels A–C for Fig 2 (with cis/no cis facets)

make_panel2 <- function(df, feature, y_max, valid_levels) {

  df_filtered <- df %>%
    group_by(UniProt_ID, HGNC_Symbol, group_lbl, pqtl_group) %>%
    summarise(
      feature_val = mode_function(.data[[feature]]),
      .groups = "drop"
    ) %>%
    filter(!is.na(feature_val), feature_val %in% valid_levels)

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

  df_filtered <- df_filtered %>%
    mutate(feature_val = factor(feature_val, levels = valid_levels))

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
      proteins = n_distinct(UniProt_ID),
      genes    = n_distinct(HGNC_Symbol),
      label    = paste0("N=", proteins, " proteins\n(", genes, " genes)"),
      .groups  = "drop"
    )

  ggplot(panel_data, aes(x = group_lbl, y = percent, fill = feature_val)) +
    geom_col(alpha = 0.9) +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 3.5,
      fontface = "bold",
      colour = "white"
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
    scale_y_continuous(limits = c(0, y_max + 10)) +
    scale_fill_brewer(palette = "Set2", name = NULL) +
    labs(x = NULL, y = "Percentage of proteins") +
    theme_bw(14) +
    theme(legend.position = "bottom")
}

p2_a <- make_panel2(
  df_uniprot2,
  "Secretion.pathway",
  110,
  c("intracellular", "membrane", "secreted")
)

p2_b <- make_panel2(
  df_uniprot2,
  "RNA.tissue.distribution",
  105,
  c("Detected in all", "Detected in many", "Detected in some", "Detected in single")
)

p2_c <- make_panel2(
  df_uniprot2,
  "RNA.tissue.specificity",
  105,
  c("Group enriched", "Low tissue specificity", "Tissue enhanced", "Tissue enriched")
)

## 3.4 Left blood concentration (log10) stratified by cis/no cis

p2_d_blood <- df_uniprot2 %>%
  filter(!is.na(Blood_conc_pgL), Blood_conc_pgL > 0) %>%
  ggplot(aes(x = group_lbl, y = Blood_conc_pgL, fill = group_lbl)) +
  geom_violin(alpha = 0.4, position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.15,
               position = position_dodge(width = 0.8),
               outlier.size = 0.4) +
  facet_wrap(~ pqtl_group, nrow = 1) +
  scale_y_log10() +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  labs(x = NULL, y = "Expected concentration [pg/L] (log10)") +
  theme_bw(14) +
  theme(legend.position = "none")

## 3.5 Aptamer-level base df for Fig 2 (RIGHT, unique SeqID + cis)

qc_base <- table_annot %>%
  distinct(SeqID, .keep_all = TRUE) %>%
  left_join(cis_by_aptamer, by = "SeqID") %>%
  mutate(
    has_cis_pqtl = ifelse(is.na(has_cis_pqtl), FALSE, has_cis_pqtl),
    pqtl_group   = factor(
      ifelse(has_cis_pqtl, "Has ≥1 cis-pQTL", "No cis-pQTL"),
      levels = c("No cis-pQTL", "Has ≥1 cis-pQTL")
    ),
    group_lbl    = factor(
      ifelse(aptamer_new_in_somascan7k_vs5k, "New in 7k", "Present in 5k"),
      levels = c("Present in 5k", "New in 7k")
    )
  )

cat("Aptamer-level QC availability (unique SeqID, Fig 2):\n")
qc_availability <- qc_base %>%
  group_by(pqtl_group, group_lbl) %>%
  summarise(
    n_aptamers_total            = n_distinct(SeqID),
    n_LOD_nonmiss               = n_distinct(SeqID[!is.na(INTERVAL_LOD_batch1) |
                                                   !is.na(INTERVAL_LOD_batch2) |
                                                   !is.na(CHRIS_LOD)]),
    n_Percent_under_LOD_nonmiss = n_distinct(SeqID[!is.na(INTERVAL_Percent_under_LOD_batch1) |
                                                   !is.na(INTERVAL_Percent_under_LOD_batch2) |
                                                   !is.na(CHRIS_Percent_under_LOD)]),
    n_CV_nonmiss                = n_distinct(SeqID[!is.na(INTERVAL_CV_batch1) |
                                                   !is.na(INTERVAL_CV_batch2) |
                                                   !is.na(CHRIS_CV)]),
    .groups = "drop"
  )
print(qc_availability)
cat("\n")

## 3.6 QC violin+box plots (RIGHT, cis-stratified)

make_qc_plot2 <- function(df, value_cols, ylab) {

  long_df <- df %>%
    select(SeqID, UniProt_ID, HGNC_Symbol, group_lbl, pqtl_group, all_of(value_cols)) %>%
    pivot_longer(
      cols      = all_of(value_cols),
      names_to  = "Batch",
      values_to = "value"
    ) %>%
    mutate(
      Batch = recode(
        Batch,
        "INTERVAL_LOD_batch1"               = "INTERVAL Batch 1",
        "INTERVAL_LOD_batch2"               = "INTERVAL Batch 2",
        "CHRIS_LOD"                         = "CHRIS",
        "INTERVAL_CV_batch1"                = "INTERVAL Batch 1",
        "INTERVAL_CV_batch2"                = "INTERVAL Batch 2",
        "CHRIS_CV"                          = "CHRIS",
        "INTERVAL_Percent_under_LOD_batch1" = "INTERVAL Batch 1",
        "INTERVAL_Percent_under_LOD_batch2" = "INTERVAL Batch 2",
        "CHRIS_Percent_under_LOD"           = "CHRIS"
      )
    ) %>%
    filter(!is.na(value)) %>%
    mutate(value_log = log10(value + 0.01))

  if (nrow(long_df) == 0) {
    return(
      ggplot() + geom_blank() +
        labs(x = NULL, y = paste0(ylab, " (log10, +0.01 offset)")) +
        theme_bw(14)
    )
  }

  counts <- df %>%
    group_by(pqtl_group, group_lbl) %>%
    summarise(
      proteins = n_distinct(UniProt_ID),
      aptamers = n_distinct(SeqID),
      label    = paste0("N=", proteins, " proteins\n(", aptamers, " aptamers)"),
      .groups  = "drop"
    )

  max_y <- max(long_df$value_log, na.rm = TRUE)

  ggplot(long_df, aes(x = group_lbl, y = value_log, fill = Batch)) +
    geom_violin(scale = "area", alpha = 0.5, position = position_dodge(0.9)) +
    geom_boxplot(width = 0.1,
                 position = position_dodge(0.9),
                 outlier.size = 0.3) +
    facet_wrap(~ pqtl_group, nrow = 1) +
    scale_fill_manual(values = batch_colors, name = NULL) +
    geom_text(
      data = counts,
      aes(x = group_lbl, y = max_y, label = label),
      inherit.aes = FALSE,
      size        = 2.8,
      fontface    = "italic",
      vjust       = 1.2
    ) +
    labs(x = NULL, y = paste0(ylab, " (log10, +0.01 offset)")) +
    theme_bw(14) +
    theme(legend.position = "bottom")
}

p2_d_qc <- make_qc_plot2(
  qc_base,
  c("INTERVAL_LOD_batch1", "INTERVAL_LOD_batch2", "CHRIS_LOD"),
  "LOD"
)

p2_e_qc <- make_qc_plot2(
  qc_base,
  c("INTERVAL_CV_batch1", "INTERVAL_CV_batch2", "CHRIS_CV"),
  "CV (%)"
)

p2_f_qc <- make_qc_plot2(
  qc_base,
  c("INTERVAL_Percent_under_LOD_batch1",
    "INTERVAL_Percent_under_LOD_batch2",
    "CHRIS_Percent_under_LOD"),
  "% under LOD"
)

## 3.7 Wilcoxon tests new vs old within each cis stratum (aptamer-level)

cat("=========================================================\n")
cat("FIGURE 2 – WILCOXON TESTS (APTAMER-LEVEL, BY cis STRATUM)\n")
cat("=========================================================\n")

wilcox_qc <- function(var, label) {
  cat("\nVariable:", label, "\n")
  for (pg in levels(qc_base$pqtl_group)) {
    cat("  Stratum:", pg, "\n")
    df_sub <- qc_base %>%
      filter(pqtl_group == pg) %>%
      select(group_lbl, all_of(var)) %>%
      filter(!is.na(.data[[var]]))
    if (n_distinct(df_sub$group_lbl) < 2) {
      cat("    Not enough groups with data.\n")
    } else {
      w <- wilcox.test(
        df_sub[[var]] ~ df_sub$group_lbl,
        exact = FALSE
      )
      cat("    n =", nrow(df_sub),
          "| p =", signif(w$p.value, 3), "\n")
    }
  }
}

wilcox_qc("INTERVAL_LOD_batch1",        "INTERVAL LOD batch 1")
wilcox_qc("INTERVAL_LOD_batch2",        "INTERVAL LOD batch 2")
wilcox_qc("CHRIS_LOD",                  "CHRIS LOD")
wilcox_qc("INTERVAL_CV_batch1",         "INTERVAL CV batch 1")
wilcox_qc("INTERVAL_CV_batch2",         "INTERVAL CV batch 2")
wilcox_qc("CHRIS_CV",                   "CHRIS CV")
wilcox_qc("INTERVAL_Percent_under_LOD_batch1", "INTERVAL % under LOD batch 1")
wilcox_qc("INTERVAL_Percent_under_LOD_batch2", "INTERVAL % under LOD batch 2")
wilcox_qc("CHRIS_Percent_under_LOD",    "CHRIS % under LOD")

## 3.8 Layout & save Figure 2

p2_a       <- p2_a + theme(plot.margin = margin(8, 10, 12, 10))
p2_b       <- p2_b + theme(plot.margin = margin(8, 10, 12, 10))
p2_c       <- p2_c + theme(plot.margin = margin(8, 10, 12, 10))
p2_d_blood <- p2_d_blood + theme(plot.margin = margin(8, 10, 12, 10))

row1_f2 <- plot_grid(p2_a, p2_d_qc,
                     ncol = 2, rel_widths = c(1.4, 1),
                     labels = c("a.", "d."),
                     label_size = 14, label_fontface = "bold")
row2_f2 <- plot_grid(p2_b, p2_e_qc,
                     ncol = 2, rel_widths = c(1.4, 1),
                     labels = c("b.", "e."),
                     label_size = 14, label_fontface = "bold")
row3_f2 <- plot_grid(p2_c, p2_f_qc,
                     ncol = 2, rel_widths = c(1.4, 1),
                     labels = c("c.", "f."),
                     label_size = 14, label_fontface = "bold")
row4_f2 <- plot_grid(p2_d_blood, NULL,
                     ncol = 2, rel_widths = c(1.4, 1),
                     labels = c("g.", ""),
                     label_size = 14, label_fontface = "bold")

figure2 <- plot_grid(row1_f2, row2_f2, row3_f2, row4_f2,
                     ncol = 1, rel_heights = c(1.05, 1.05, 1.05, 0.8))

ggsave("Figure_2_cis_pQTL_stratified_FINAL_RAW.png",
       figure2, width = 13, height = 15, dpi = 300)

cat("\nSaved: Figure_2_cis_pQTL_stratified_FINAL_RAW.png\n")
cat("=========================================================\n")
cat("Script finished. All figures and console stats use RAW tables (no LOD filtering).\n")
cat("You can now compare these numbers with the manuscript text and adjust wording/QC description.\n")
cat("=========================================================\n")

