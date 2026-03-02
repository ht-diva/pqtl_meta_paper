rm(list=ls())
library(data.table)
library(dplyr)
library(tidyr)
library(Rmpfr)
library(stringr)
library(readxl)
library(ggplot2)

file1 <- read.delim("/group/diangelantonio/users/giulia.pontali/file_to_share/unconditional_IVs_from_MVP_to_be_merge_w_new_results_pvalNA.txt")

no_phecode <- "/exchange/healthds/pQTL/results/META_CHRIS_INTERVAL/MVP_MR/new_results_20250515/INTERVAL_CHRIS_MR_COJO_NONPHECODE_07MAY2025.txt.gz"
phecode <- "/exchange/healthds/pQTL/results/META_CHRIS_INTERVAL/MVP_MR/new_results_20250515/INTERVAL_CHRIS_MR_COJO_PHECODE_07MAY2025.txt.gz"
data_no_phecode <- fread(no_phecode, sep = "\t")
data_phecode <- fread(phecode, sep = "\t")

data <- rbind(data_no_phecode, data_phecode)
data <- as.data.frame(data)
new_ivs <- read.delim("/group/diangelantonio/users/giulia.pontali/file_to_share/unconditional_IVs_to_MVP_to_rerun.txt")
new_ivs$CHRPOS_ID <- paste0(new_ivs$CHR,"_",new_ivs$POS_38,"_",new_ivs$SeqID)
data <- merge(data, new_ivs[, c("CHRPOS_ID", "locus_START_END_37", "SNP", "POS_37")], by = "CHRPOS_ID", all.x = TRUE)
data$locus_START_END_37 <- gsub("-", "_", data$locus_START_END_37)
data$locus_START_END_37 <- paste0("chr", data$CHR, "_", data$locus_START_END_37)
class(data)


colnames(file1)[colnames(file1) == "UniProt_ID"] <- "UNIPROT"
file1 <- file1[,c(-2,-3)]
colnames(file1)[colnames(file1) == "symbol"] <- "HARMONIZED_GENE_NAME"
colnames(file1)
file1 <- file1[,c(-26)]
colnames(data)
data <- data[,c(-3,-6,-28)]

file1 <- file1[, sort(colnames(file1))]
data <- data[, sort(colnames(data))]

# Combine
final <- rbind(file1, data)

dim(final)
head(final)

final <- final[!grepl("Height", final$PHENOTYPE), ]
final$seqID <- sub(".*(seq\\.\\d+\\.\\d+).*", "\\1", final$CHRPOS_ID)
unique_triplets <- final %>% distinct(seqID, UNIPROT, PHENOTYPE) %>% nrow()
final$sign <- final$PVAL < (0.05/3497303)

table2 <- read_xlsx("/home/giulia.pontali/supplementary_table_2 (7).xlsx", sheet=2)
ann_uniprot <- table2[,c(5,37)]
ann_uniprot <- unique(ann_uniprot, by = "UniProt_ID")

ann_seqid <- table2[,c(1,33)]
ann_seqid <- unique(ann_seqid, by = "SeqID")

# Left-join UniProt annotation onto final
final2 <- merge(
  final,
  ann_uniprot,
  by.x = "UNIPROT",
  by.y = "UniProt_ID",
  all.x = TRUE,
  sort = FALSE
)

# Left-join SeqID annotation onto final2
final3 <- merge(
  final2,
  ann_seqid,
  by.x = "seqID",
  by.y = "SeqID",
  all.x = TRUE,
  sort = FALSE
)

final <- final3
table4 <- read_xlsx("/home/giulia.pontali/supplementary_table_4 (14).xlsx", sheet=2, skip=2)
table4 <- table4[,c("SeqID", "SNPID","tissues_with_encoding_gene_coloc", "encoding_gene_coloc_in_region")]

merged <- merge(
  x = final,
  y = table4,
  by.x = c("seqID", "SNP"),
  by.y = c("SeqID", "SNPID"),
  all.x = TRUE
)

final <- merged

setDT(final)
collapsed_final<- final[, lapply(.SD, function(x) paste(unique(x), collapse = "|")), 
                                    by = .(seqID, locus_START_END_37, BETA, PROTEIN_NAME, PHENOTYPE, CHR, PVAL, SE)]

write.table(collapsed_final, "/scratch/giulia.pontali/meta_file/collapsed_final.txt", sep="\t", quote=F, row.names = F)

#---
combined_data_epitope_coloc <- fread("/scratch/giulia.pontali/meta_file/collapsed_final.txt")
combined_data_epitope_coloc <- combined_data_epitope_coloc[which(combined_data_epitope_coloc$sign == TRUE), ]

combined_data_epitope_coloc <- combined_data_epitope_coloc %>%
  mutate(
    HAS_INT = str_detect(PHENOTYPE, "(?i)_INT$"),
    CORE = str_remove(PHENOTYPE, "(?i)_INT$"),
    
    CORE = CORE %>%
      str_replace("(?i)^mean_([^_]+)$",   "\\1_Mean") %>%
      str_replace("(?i)^min_([^_]+)$",    "\\1_Min") %>%
      str_replace("(?i)^max_([^_]+)$",    "\\1_Max") %>%
      str_replace("(?i)^median_([^_]+)$", "\\1_Median") %>%
      
      str_replace("(?i)_mean$",   "_Mean") %>%
      str_replace("(?i)_min$",    "_Min") %>%
      str_replace("(?i)_max$",    "_Max") %>%
      str_replace("(?i)_median$", "_Median"),
    
    PHENOTYPE = if_else(HAS_INT, paste0(CORE, "_INT"), CORE),
    PHENOTYPE = str_replace_all(PHENOTYPE, "__", "_")
  ) %>%
  select(-HAS_INT, -CORE)

# ----------------------------
# 1) Join EFO annotations
# ----------------------------
efo <- read_excel("/scratch/giulia.pontali/pheno_updates.xlsx")
efo <- efo[, c(1,2,17:20)]

combined_data_epitope_coloc <- left_join(
  combined_data_epitope_coloc,
  efo %>% select(Phenotype, MVP_description, EFO_ID, EFO_Term, EFO_Parent_ID, EFO_Parent_Term),
  by = c("PHENOTYPE" = "Phenotype")
)

# ----------------------------
# 2) Load + normalize drug file
# ----------------------------
drug_indication <- read.delim("/exchange/healthds/MVP/061624_Drug_Indication_clean_hv1.2_11JUL2024.tsv")

drug_indication <- drug_indication %>%
  mutate(across(
    c(Genetic_phenotype_phecode, EFO_ID, EFO_Term, EFO_Parent_Term, EFO_Parent_ID, Matched_genetic_phenotype),
    ~ str_replace_na(.x, "")
  )) %>%
  separate_rows(Genetic_phenotype_phecode, sep="\\|", convert=FALSE) %>%
  separate_rows(EFO_ID, sep="\\|", convert=FALSE) %>%
  separate_rows(EFO_Term, sep="\\|", convert=FALSE) %>%
  separate_rows(EFO_Parent_ID, sep="\\|", convert=FALSE) %>%
  separate_rows(Matched_genetic_phenotype, sep="\\|", convert=FALSE) %>%
  mutate(across(
    c(Genetic_phenotype_phecode, EFO_ID, EFO_Term, EFO_Parent_Term, EFO_Parent_ID, Matched_genetic_phenotype),
    str_trim
  ))

# ----------------------------
# 3) drug_summary used for classification (>=3)
# ----------------------------
drug_summary <- drug_indication %>%
  select(
    HARMONIZED_GENE_NAME, drug_indication, drug_name, moa, drug_max_phase, max_phase_for_ind,
    parent_drug_name, Indication_disease_name, year_of_first_approval,
    Genetic_phenotype_phecode, EFO_ID, EFO_Term, EFO_Parent_ID, EFO_Parent_Term,
    Matched_genetic_phenotype, action_type
  ) %>%
  distinct() %>%
  filter(max_phase_for_ind >= 3) %>%  # <-- classification rule
  separate_rows(EFO_Parent_Term, sep="\\|", convert=FALSE) %>%
  mutate(
    EFO_ID_lower     = tolower(EFO_ID),
    EFO_Term_lower   = tolower(EFO_Term),
    EFO_Parent_Term_lower = tolower(EFO_Parent_Term),
    EFO_Parent_ID_lower = tolower(EFO_Parent_ID),
    Matched_genetic_phenotype_lower = tolower(Matched_genetic_phenotype),
    Genetic_phenotype_phecode_lower = tolower(Genetic_phenotype_phecode)
  ) %>%
  select(-EFO_ID, -EFO_Term, -EFO_Parent_Term, -EFO_Parent_ID,
         -Matched_genetic_phenotype, -Genetic_phenotype_phecode)

drug_summary <- fread("/scratch/giulia.pontali/meta_file/drug_summary.txt")
# ----------------------------
# 4) Lookup all drugs/indications for the gene (same as your structure)
#    (NOTE: based on drug_summary => only >=3 drugs are included, consistent with your rule)
# ----------------------------
gene_drug_info <- drug_summary %>%
  group_by(HARMONIZED_GENE_NAME) %>%
  summarise(
    all_drug_names = paste(unique(drug_name), collapse = " | "),
    all_indications = paste(unique(Indication_disease_name), collapse = " | "),
    .groups = "drop"
  )

# ----------------------------
# 5) Prepare MR table
# ----------------------------
combined_data_epitope_coloc <- combined_data_epitope_coloc %>%
  mutate(row_id = row_number()) %>%
  mutate(across(c(HARMONIZED_GENE_NAME, PHENOTYPE, EFO_Term, EFO_ID, EFO_Parent_Term, EFO_Parent_ID),
                ~ str_replace_na(.x, ""))) %>%
  mutate(across(c(HARMONIZED_GENE_NAME, PHENOTYPE, EFO_Term, EFO_ID, EFO_Parent_Term, EFO_Parent_ID),
                str_trim)) %>%
  separate_rows(EFO_ID, sep="\\|") %>%
  separate_rows(EFO_Term, sep="\\|") %>%
  separate_rows(EFO_Parent_Term, sep="\\|") %>%
  separate_rows(EFO_Parent_ID, sep="\\|") %>%
  mutate(
    EFO_ID_lower     = str_trim(tolower(EFO_ID)),
    EFO_Term_lower   = str_trim(tolower(EFO_Term)),
    EFO_Parent_Term_lower = str_trim(tolower(EFO_Parent_Term)),
    EFO_Parent_ID_lower = str_trim(tolower(EFO_Parent_ID)),
    PHENOTYPE_lower  = str_trim(tolower(PHENOTYPE))
  ) %>%
  select(-EFO_ID, -EFO_Term, -EFO_Parent_Term, -EFO_Parent_ID)

# ----------------------------
# 6) Matching steps (unchanged logic)
# ----------------------------

### 1) Exact EFO ID
exact_id <- combined_data_epitope_coloc %>%
  inner_join(
    drug_summary %>%
      select(HARMONIZED_GENE_NAME, EFO_ID_lower, drug_name, action_type,
             Indication_disease_name, max_phase_for_ind),
    by = c("HARMONIZED_GENE_NAME", "EFO_ID_lower")
  ) %>%
  left_join(gene_drug_info, by = "HARMONIZED_GENE_NAME") %>%
  mutate(
    drug_rep = "Drug exists (Exact EFO_ID)",
    match_type = "exact_efo_id",
    drug_name = all_drug_names,
    Indication_disease_name = all_indications
  ) %>%
  select(-all_drug_names, -all_indications)

### 2) Exact EFO Term
exact_term <- combined_data_epitope_coloc %>%
  anti_join(exact_id %>% select(row_id), by = "row_id") %>%
  inner_join(
    drug_summary %>%
      select(HARMONIZED_GENE_NAME, EFO_Term_lower, drug_name, action_type,
             Indication_disease_name, max_phase_for_ind),
    by = c("HARMONIZED_GENE_NAME", "EFO_Term_lower")
  ) %>%
  left_join(gene_drug_info, by = "HARMONIZED_GENE_NAME") %>%
  mutate(
    drug_rep = "Drug exists (Exact EFO_Term)",
    match_type = "exact_efo_term",
    drug_name = all_drug_names,
    Indication_disease_name = all_indications
  ) %>%
  select(-all_drug_names, -all_indications)

### 3) Parent EFO
parent_match <- combined_data_epitope_coloc %>%
  anti_join(bind_rows(exact_id, exact_term) %>% select(row_id), by = "row_id") %>%
  inner_join(
    drug_summary %>%
      select(HARMONIZED_GENE_NAME, EFO_Parent_ID_lower, drug_name, action_type,
             Indication_disease_name, max_phase_for_ind),
    by = c("HARMONIZED_GENE_NAME", "EFO_Parent_ID_lower")
  ) %>%
  left_join(gene_drug_info, by = "HARMONIZED_GENE_NAME") %>%
  mutate(
    drug_rep = "Drug exists (Parent EFO)",
    match_type = "parent_efo",
    drug_name = all_drug_names,
    Indication_disease_name = all_indications
  ) %>%
  select(-all_drug_names, -all_indications)

### 4) Genetic phenotype phecode
phecode_match <- combined_data_epitope_coloc %>%
  anti_join(bind_rows(exact_id, exact_term, parent_match) %>% select(row_id), by = "row_id") %>%
  inner_join(
    drug_summary %>%
      select(HARMONIZED_GENE_NAME, Genetic_phenotype_phecode_lower, drug_name,
             action_type, Indication_disease_name, max_phase_for_ind),
    by = c("HARMONIZED_GENE_NAME", "PHENOTYPE_lower" = "Genetic_phenotype_phecode_lower")
  ) %>%
  left_join(gene_drug_info, by = "HARMONIZED_GENE_NAME") %>%
  mutate(
    drug_rep = "Drug exists (Genetic_phenotype_phecode)",
    match_type = "phecode",
    drug_name = all_drug_names,
    Indication_disease_name = all_indications
  ) %>%
  select(-all_drug_names, -all_indications)

### 5) Matched genetic phenotype
matched_pheno <- combined_data_epitope_coloc %>%
  anti_join(bind_rows(exact_id, exact_term, parent_match, phecode_match) %>% select(row_id), by = "row_id") %>%
  inner_join(
    drug_summary %>%
      select(HARMONIZED_GENE_NAME, Matched_genetic_phenotype_lower, drug_name,
             action_type, Indication_disease_name, max_phase_for_ind),
    by = "HARMONIZED_GENE_NAME"
  ) %>%
  filter(PHENOTYPE_lower == Matched_genetic_phenotype_lower) %>%
  left_join(gene_drug_info, by = "HARMONIZED_GENE_NAME") %>%
  mutate(
    drug_rep = "Drug exists (Matched_genetic_phenotype)",
    match_type = "matched_pheno",
    drug_name = all_drug_names,
    Indication_disease_name = all_indications
  ) %>%
  select(-all_drug_names, -all_indications)

### 6) Drug exists for gene (repurposing)
has_drug_gene <- combined_data_epitope_coloc %>%
  anti_join(bind_rows(exact_id, exact_term, parent_match, phecode_match, matched_pheno) %>% select(row_id), by = "row_id") %>%
  semi_join(drug_summary %>% select(HARMONIZED_GENE_NAME) %>% distinct(), by = "HARMONIZED_GENE_NAME") %>%
  left_join(
    drug_summary %>%
      group_by(HARMONIZED_GENE_NAME) %>%
      summarise(
        drug_name = paste(unique(drug_name), collapse = " | "),
        action_type = paste(unique(action_type), collapse = " | "),
        Indication_disease_name = paste(unique(Indication_disease_name), collapse = " | "),
        max_phase_for_ind = paste(unique(max_phase_for_ind), collapse = " | "),
        .groups = "drop"
      ),
    by = "HARMONIZED_GENE_NAME"
  ) %>%
  mutate(
    drug_rep = "Drug repurposing",
    match_type = "none"
  )

### 7) No drug for gene
### 7) No drug for gene (fix type mismatch by dropping max_phase_for_ind)
no_drug_gene <- combined_data_epitope_coloc %>%
  anti_join(
    bind_rows(
      exact_id      %>% select(-any_of("max_phase_for_ind")),
      exact_term    %>% select(-any_of("max_phase_for_ind")),
      parent_match  %>% select(-any_of("max_phase_for_ind")),
      phecode_match %>% select(-any_of("max_phase_for_ind")),
      matched_pheno %>% select(-any_of("max_phase_for_ind")),
      has_drug_gene %>% select(-any_of("max_phase_for_ind"))
    ) %>% select(row_id),
    by = "row_id"
  ) %>%
  mutate(
    drug_rep = "No drug-gene match",
    match_type = "no_drug_gene",
    drug_name = NA_character_,
    Indication_disease_name = NA_character_,
    action_type = NA_character_
  )
# ----------------------------
# 7) Final table + IMPORTANT FIXES:
#    - fill drug_rep always
#    - remove max_phase_for_ind
# ----------------------------
drop_phase <- function(df) df %>% dplyr::select(-dplyr::any_of("max_phase_for_ind"))

classified <- dplyr::bind_rows(
  drop_phase(exact_id),
  drop_phase(exact_term),
  drop_phase(parent_match),
  drop_phase(phecode_match),
  drop_phase(matched_pheno),
  drop_phase(has_drug_gene),
  drop_phase(no_drug_gene)
) %>%
  dplyr::mutate(
    drug_rep = dplyr::coalesce(drug_rep, "No drug-gene match"),
    drug_rep = ifelse(trimws(drug_rep) == "", "No drug-gene match", drug_rep)
  )

# Ensure row_id exists (ideally it should already exist upstream)
if (!"row_id" %in% names(classified)) {
  classified <- classified %>% dplyr::mutate(row_id = dplyr::row_number())
}

# Define static columns (must exist or will be ignored by any_of)
cols_statiche <- c("row_id", "HARMONIZED_GENE_NAME", "PHENOTYPE", "PHENOTYPE_lower", "seqID", "SNP", "UNIPROT")

classified_collapsed <- classified %>%
  dplyr::group_by(row_id) %>%
  dplyr::summarise(
    dplyr::across(dplyr::any_of(cols_statiche), ~ dplyr::first(.x)),
    dplyr::across(
      setdiff(names(.), cols_statiche),   # collapse everything else
      ~ paste(unique(.x[!is.na(.x) & .x != ""]), collapse = " | ")
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    dplyr::across(dplyr::where(is.character), ~ na_if(.x, "")),
    drug_rep = dplyr::coalesce(drug_rep, "No drug-gene match")
  )


write.table(
  classified_collapsed,
  file = "/scratch/giulia.pontali/meta_file/classified_collapsed.txt",
  sep="\t", quote=F, row.names = F
)

write.xlsx(
  classified_collapsed,
  "/scratch/giulia.pontali/meta_file/classified_collapsed.xlsx"
)
