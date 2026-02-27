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
combined_data_epitope_coloc <- combined_data_epitope_coloc[which(combined_data_epitope_coloc$sign==TRUE),]

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



efo <- read_excel("/scratch/giulia.pontali/pheno_updates.xlsx")
efo <- efo[,c(1,2,17:20)]

combined_data_epitope_coloc <- left_join(
  combined_data_epitope_coloc,
  efo %>% select(Phenotype, MVP_description, EFO_ID, EFO_Term, EFO_Parent_ID, EFO_Parent_Term),
  by = c("PHENOTYPE" = "Phenotype")
)


# 1) Load + normalize drug file
drug_indication <- read.delim("/exchange/healthds/MVP/061624_Drug_Indication_clean_hv1.2_11JUL2024.tsv")


drug_indication <- drug_indication %>%
  mutate(across(c(Genetic_phenotype_phecode, EFO_ID, EFO_Term, EFO_Parent_Term, EFO_Parent_ID, Matched_genetic_phenotype), 
                ~str_replace_na(.x, ""))) %>%
  separate_rows(Genetic_phenotype_phecode, sep="\\|", convert=FALSE) %>%
  separate_rows(EFO_ID, sep="\\|", convert=FALSE) %>%
  separate_rows(EFO_Term, sep="\\|", convert=FALSE) %>%
  #separate_rows(EFO_Parent_Term, sep="\\|", convert=FALSE) %>%
  separate_rows(EFO_Parent_ID, sep="\\|", convert=FALSE) %>%
  separate_rows(Matched_genetic_phenotype, sep="\\|", convert=FALSE) %>%
  mutate(across(c(Genetic_phenotype_phecode, EFO_ID, EFO_Term, EFO_Parent_Term, EFO_Parent_ID, Matched_genetic_phenotype), str_trim))


drug_summary <- drug_indication %>%
  select(
    HARMONIZED_GENE_NAME, drug_indication, drug_name, moa,
    drug_max_phase, max_phase_for_ind, year_of_first_approval,
    parent_drug_name, Indication_disease_name,
    Genetic_phenotype_phecode, EFO_ID, EFO_Term, EFO_Parent_ID, EFO_Parent_Term,
    Matched_genetic_phenotype, action_type
  ) %>%
  distinct() %>%
  filter(max_phase_for_ind >= 3) %>%
  separate_rows(EFO_Parent_Term, sep = "\\|") %>%
  mutate(
    EFO_ID_lower     = tolower(EFO_ID),
    EFO_Term_lower   = tolower(EFO_Term),
    EFO_Parent_Term_lower = tolower(EFO_Parent_Term),
    EFO_Parent_ID_lower = tolower(EFO_Parent_ID),
    Matched_genetic_phenotype_lower = tolower(Matched_genetic_phenotype),
    Genetic_phenotype_phecode_lower = tolower(Genetic_phenotype_phecode)
  ) %>%
  select(
    -EFO_ID, -EFO_Term, -EFO_Parent_Term, -EFO_Parent_ID,
    -Matched_genetic_phenotype, -Genetic_phenotype_phecode
  )

gene_drug_info <- drug_summary %>%
  group_by(HARMONIZED_GENE_NAME) %>%
  summarise(
    all_drug_names = paste(unique(drug_name), collapse = " | "),
    all_indications = paste(unique(Indication_disease_name), collapse = " | "),
    .groups = "drop"
  )

combined_data_epitope_coloc <- combined_data_epitope_coloc %>%
  mutate(row_id = row_number()) %>%
  mutate(across(c(HARMONIZED_GENE_NAME, PHENOTYPE, EFO_Term, EFO_ID, EFO_Parent_Term, EFO_Parent_ID),
                ~ str_replace_na(.x, ""))) %>%
  mutate(across(c(HARMONIZED_GENE_NAME, PHENOTYPE, EFO_Term, EFO_ID, EFO_Parent_Term, EFO_Parent_ID),
                str_trim)) %>%
  separate_rows(EFO_ID, sep = "\\|") %>%
  separate_rows(EFO_Term, sep = "\\|") %>%
  separate_rows(EFO_Parent_Term, sep = "\\|") %>%
  separate_rows(EFO_Parent_ID, sep = "\\|") %>%
  mutate(
    EFO_ID_lower     = str_trim(tolower(EFO_ID)),
    EFO_Term_lower   = str_trim(tolower(EFO_Term)),
    EFO_Parent_Term_lower = str_trim(tolower(EFO_Parent_Term)),
    EFO_Parent_ID_lower = str_trim(tolower(EFO_Parent_ID)),
    PHENOTYPE_lower  = str_trim(tolower(PHENOTYPE))
  ) %>%
  select(-EFO_ID, -EFO_Term, -EFO_Parent_Term, -EFO_Parent_ID)

exact_id <- combined_data_epitope_coloc %>%
  inner_join(
    drug_summary %>%
      select(
        HARMONIZED_GENE_NAME, EFO_ID_lower,
        drug_name, action_type, Indication_disease_name,
        drug_max_phase, max_phase_for_ind, year_of_first_approval
      ),
    by = c("HARMONIZED_GENE_NAME", "EFO_ID_lower")
  ) %>%
  left_join(gene_drug_info, by = "HARMONIZED_GENE_NAME") %>%
  mutate(
    drug_rep = "Drug exists (Exact EFO_ID)",
    match_type = "exact_efo_id",
    drug_name = all_drug_names,
    Indication_disease_name = all_indications
  ) %>%
  select(-all_drug_names, -all_indications) %>%
  distinct(row_id, .keep_all = TRUE)


exact_term <- combined_data_epitope_coloc %>%
  anti_join(exact_id %>% select(row_id), by = "row_id") %>%
  inner_join(
    drug_summary %>%
      select(
        HARMONIZED_GENE_NAME, EFO_Term_lower,
        drug_name, action_type, Indication_disease_name,
        drug_max_phase, max_phase_for_ind, year_of_first_approval
      ),
    by = c("HARMONIZED_GENE_NAME", "EFO_Term_lower")
  ) %>%
  left_join(gene_drug_info, by = "HARMONIZED_GENE_NAME") %>%
  mutate(
    drug_rep = "Drug exists (Exact EFO_Term)",
    match_type = "exact_efo_term",
    drug_name = all_drug_names,
    Indication_disease_name = all_indications
  ) %>%
  select(-all_drug_names, -all_indications) %>%
  distinct(row_id, .keep_all = TRUE)


parent_match <- combined_data_epitope_coloc %>%
  anti_join(bind_rows(exact_id, exact_term) %>% select(row_id), by = "row_id") %>%
  inner_join(
    drug_summary %>%
      select(
        HARMONIZED_GENE_NAME, EFO_Parent_ID_lower,
        drug_name, action_type, Indication_disease_name,
        drug_max_phase, max_phase_for_ind, year_of_first_approval
      ),
    by = c("HARMONIZED_GENE_NAME", "EFO_Parent_ID_lower")
  ) %>%
  left_join(gene_drug_info, by = "HARMONIZED_GENE_NAME") %>%
  mutate(
    drug_rep = "Drug exists (Parent EFO)",
    match_type = "parent_efo",
    drug_name = all_drug_names,
    Indication_disease_name = all_indications
  ) %>%
  select(-all_drug_names, -all_indications) %>%
  distinct(row_id, .keep_all = TRUE)

phecode_match <- combined_data_epitope_coloc %>%
  anti_join(bind_rows(exact_id, exact_term, parent_match) %>% select(row_id), by = "row_id") %>%
  inner_join(
    drug_summary %>%
      select(
        HARMONIZED_GENE_NAME, Genetic_phenotype_phecode_lower,
        drug_name, action_type, Indication_disease_name,
        drug_max_phase, max_phase_for_ind, year_of_first_approval
      ),
    by = c("HARMONIZED_GENE_NAME", "PHENOTYPE_lower" = "Genetic_phenotype_phecode_lower")
  ) %>%
  left_join(gene_drug_info, by = "HARMONIZED_GENE_NAME") %>%
  mutate(
    drug_rep = "Drug exists (Genetic_phenotype_phecode)",
    match_type = "phecode",
    drug_name = all_drug_names,
    Indication_disease_name = all_indications
  ) %>%
  select(-all_drug_names, -all_indications) %>%
  distinct(row_id, .keep_all = TRUE)

matched_pheno <- combined_data_epitope_coloc %>%
  anti_join(bind_rows(exact_id, exact_term, parent_match, phecode_match) %>% select(row_id), by = "row_id") %>%
  inner_join(
    drug_summary %>%
      select(
        HARMONIZED_GENE_NAME, Matched_genetic_phenotype_lower,
        drug_name, action_type, Indication_disease_name,
        drug_max_phase, max_phase_for_ind, year_of_first_approval
      ),
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
  select(-all_drug_names, -all_indications) %>%
  distinct(row_id, .keep_all = TRUE)


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
        drug_max_phase = paste(unique(drug_max_phase), collapse = " | "),
        max_phase_for_ind = paste(unique(max_phase_for_ind), collapse = " | "),
        year_of_first_approval = paste(unique(year_of_first_approval), collapse = " | "),
        .groups = "drop"
      ),
    by = "HARMONIZED_GENE_NAME"
  ) %>%
  mutate(
    drug_rep = "Drug repurposing",
    match_type = "none"
  ) %>%
  distinct(row_id, .keep_all = TRUE)

dfs <- list(exact_id, exact_term, parent_match, phecode_match, matched_pheno, has_drug_gene)

# Trova le colonne comuni
common_cols <- Reduce(intersect, lapply(dfs, names))

# Escludi row_id dalla conversione
cols_to_convert <- setdiff(common_cols, "row_id")

dfs_norm <- lapply(dfs, function(df) {
  df %>% mutate(across(all_of(cols_to_convert), as.character))
})

no_drug_gene <- combined_data_epitope_coloc %>%
  anti_join(
    bind_rows(dfs_norm) %>% select(row_id),
    by = "row_id"
  ) %>%
  mutate(
    drug_rep = "No drug-gene match",
    match_type = "no_drug_gene"
  ) %>%
  distinct(row_id, .keep_all = TRUE)


classified <- bind_rows(
  exact_id,
  exact_term,
  parent_match,
  phecode_match,
  matched_pheno,
  has_drug_gene,
  no_drug_gene
) %>%
  select(any_of(c(
    colnames(combined_data_epitope_coloc),
    "drug_rep", "drug_name", "action_type",
    "drug_indication", "match_type",
    "Indication_disease_name",
    "drug_max_phase", "max_phase_for_ind", "year_of_first_approval"
  ))) %>%
  distinct(row_id, .keep_all = TRUE)

no_drug_gene_norm <- no_drug_gene %>% mutate(across(any_of(cols_to_convert), as.character))
classified <- bind_rows(
  dfs_norm,
  no_drug_gene_norm
) %>%
  select(any_of(c(
    colnames(combined_data_epitope_coloc),
    "drug_rep", "drug_name", "action_type",
    "drug_indication", "match_type",
    "Indication_disease_name",
    "drug_max_phase", "max_phase_for_ind", "year_of_first_approval"
  ))) %>%
  distinct(row_id, .keep_all = TRUE)


library(dplyr)

# Create a lookup table collapsing duplicate gene matches
drug_lookup <- drug_summary %>%
  group_by(HARMONIZED_GENE_NAME) %>%
  summarise(
    Indication_disease_name = paste(unique(Indication_disease_name), collapse = " | "),
    .groups = "drop"
  )


# --- Step: Add approved_indications_all per gene–drug ---

# Create lookup with all approved indications per gene–drug
approved_lookup <- drug_summary %>%
  group_by(HARMONIZED_GENE_NAME, drug_name) %>%
  summarise(
    approved_indications_all = paste(unique(Indication_disease_name), collapse = " | "),
    .groups = "drop"
  )

# Join the lookup to add approved_indications_all
classified <- classified %>%
  left_join(approved_lookup, by = c("HARMONIZED_GENE_NAME", "drug_name"))

# For repurposing rows, clear Indication_disease_name (set to NA)
classified_a <- classified %>%
  mutate(
    Indication_disease_name = if_else(
      drug_rep == "Drug repurposing",
      NA_character_,
      Indication_disease_name
    )
  )

table2 <- read_xlsx("/home/giulia.pontali/supplementary_table_2 (7).xlsx", sheet=2)

final2 <- merge(
  classified,
  table2,
  by.x = c("seqID", "HARMONIZED_GENE_NAME"),
  by.y = c("SeqID", "HGNC_Symbol"),
  all.x = TRUE,
  sort = FALSE
)

final2 <- final2 %>%
  group_by(SNP) %>%
  mutate(n_other_phenotypes_for_SNP = n() - 1) %>%
  ungroup()


write_xlsx(final2, "/home/giulia.pontali/final_drug.xlsx")



