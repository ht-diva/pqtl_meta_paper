
# Date: February 12, 2026
# Author: Dariush Ghasemi S.


library(data.table)
library(tidyverse)


# This script is written to inspect the output of pqtl_conditional pipeline (issue #40)
# Specifically we intend to ensure that:
#   1. lifted summary stats does not have any NAs (except conditional columns)
#   2. we don't lose a big portion of the variants in the region after lifting to hg38


#-------------------------------#
# Paths
input_giulia <- "/exchange/healthds/pQTL/results/META_CHRIS_INTERVAL/bonferroni_threshold/subset_for_coloc.txt"
path_base    <- "/scratch/dariush.ghasemi/projects/pqtl_conditional/results/issue_40"
path_sumstat <- glue(path_base, "/gwas/seq.10025.1_7_107477440_107586518.tsv")
path_pos38   <- glue(path_base, "/VCF_lifted/seq.10025.1_7_107477440_107586518.txt")
path_freq    <- glue(path_base, "/freq/seq.10025.1_7_107477440_107586518.afreq")
path_example <- glue(path_base, "/toMVP/marginal_data_seq.10490.3_locus_3_127956148_128826577_target_3_128338600_A_C.tsv")
path_marginals <- glue(path_base, "/toMVP")


# Separate data file name useful for extracting SNP from RDS --- NOT USED HERE
tsv_name <- basename(path_sumstat)
seq_name <- tsv_name %>% str_extract("seq.(\\d)+.(\\d)+")
loc_name <- tsv_name %>% str_remove("_(\\d)+.(\\d)+.([ATCG])+.([ATCG])+.*$") %>% str_remove("seq.(\\d)+.(\\d)+_")
snp_name <- tsv_name %>% str_remove("seq.(\\d)+.(\\d)+_(\\d)+_(\\d)+_(\\d)+_") %>% str_remove(".tsv")


#-----------------#
# Check if input loci from Giulia are unique
fread(input_giulia) %>% distinct(SeqID, locus_START_END_37)

# Grep a random SNP from Giulia's input
fread(input_giulia) %>% filter(SNPID == "3:128338600:A:C")

# and compare values in lifted sumtats
fread(path_example) %>% filter(SNP_37 == "3:128338600:A:C")


#-----------------#
# List of marginal sumatst to transfer
marginal_list <- list.files(
  path = path_marginals, 
  pattern = ".tsv", 
  full.names = T
  )

# 1. Check if there is any missing values in the lifted sumstats
res <- map_dfr(
  marginal_list, function(df_path){
  df <- fread(df_path) %>% select(- ends_with("_cond"))
  table(is.na(df))
})


#-----------------#
# Inspect the combined report of the pipeline
reports <- fread(glue(path_marginals, "/../combined_report.tsv"))

# 2. Bar plot: Percentage of lost variants in 35 region affected by lift-over
reports %>% mutate(
    check = n_sumstat != n_liftover,
    difference =  round((n_sumstat - n_liftover) / n_sumstat, 4)*100,
    chrom = str_remove(index, ":.*$") %>% as.numeric(),
    region = fct_reorder(locus, chrom)
    ) %>%
  filter(check) %>%
  ggplot(aes(x = region, y = difference)) +
  geom_col(position = position_dodge(.9), fill = "#334567") + 
  geom_text(aes(label = difference), vjust = -.3, size = 2) +
  labs(y = "Variants lost in each region\n after liftover (%)") +
  theme_light() +
  theme(
    axis.text.x = element_text(angle = -90, vjust = .5, hjust = 0)
  )

# Save bar plot
ggsave(
  filename = "12-Feb-2026_variants_lost_in_liftover_percentage.png",
  width = 9, height = 5.5, dpi = 200
  )


#-----------------#
# Compare allele frequencies of the variants across 586 regions sumstat

read_gwas <- function(path_mvp){
  
  locuseq <- basename(path_mvp) %>% str_remove_all("marginal_data_|_locus|_target_.*$")
  path_gwas <- str_c(path_base, "/gwas/", locuseq, ".tsv")
  
  df_mvp  <- fread(path_mvp)
  df_gwas <- fread(path_gwas)
  
  # merge two datasets
  df_mvp %>%
    select(
      seqid = SeqID,
      locus = LOCUS_37,
      SNPID = SNP_37,
      EAF_MVP = EAF
      ) %>%
    left_join(
      df_gwas %>% select(SNPID, EAF_GWAS = EAF),
      join_by(SNPID)
      )
  
}

# Combine EAFs in meta-analysis and EAFs in INTERVAL for 3,011,054 variants
gwas_merged <- map_dfr(marginal_list, read_gwas)


# Scatter plot
plt_eaf <- gwas_merged %>%
  ggplot(aes(x = EAF_GWAS, y = EAF_MVP)) +
  geom_point(size = 2, alpha = .7, color = "steelblue") +
  geom_abline(slope = 1, lty = 2, size = .5, color = "grey30") +
  theme_light() + 
  labs(
    x = "Effect Allele Frequency in Meta-Analysis",
    y = "Effect Allele Frequency in INTERVAL"
  )


# Save scatter plot
ggsave(
  filename = "12-Feb-2026_allele_frequency_comparison.png",
  plot = plt_eaf, width = 8.5, height = 7, dpi = 200
)

