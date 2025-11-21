###############################################################################
# Claudia Giambartolomei
# Date: Oct 24 2025
# Decisions after meeting of Tuesday 21 October for matching historical assays: 
# Source: SomaLogic menu mapping file, SomaScan.db (v4.0), Olink HT and Olink 3K 
###############################################################################
# load libraries
library(dplyr)
library(data.table)
library(stringr)
library(AnnotationDbi)
library(SomaScan.db)
# needed for section 3 and 4
library(readr)
library(stringr)
library(tidyr)
library(openxlsx)
library(purrr)

# DATA
## --- [1] Load mapping file ---
# mapping_file_from_eva <- "/Users/c.giambartolomei/Downloads/annotations_of_proteins_oct222025/Fig1_Solene_Meta_paper-main/Results_files_for_plots/somascan_tss_ncbi_grch37_version_20250326.txt"
# mapping_file_from_eva <- "/exchange/healthds/pQTL/Reference_datasets_for_QC_proteomics/Cis_trans_mapping/somascan_tss_ncbi_grch37_version_20250326.txt"

# Giulia did this part so I do not need to replicate:
add_HGNC_Symbol = FALSE

# mapping_file <- "/exchange/healthds/pQTL/Reference_datasets_for_QC_proteomics/Cis_trans_mapping/somascan_tss_ncbi_grch37_version_20251023.txt"
mapping_file <- "/Users/c.giambartolomei/Downloads/annotations_of_proteins_oct222025/Fig1_Solene_Meta_paper-main/Results_files_for_plots/somascan_tss_ncbi_grch37_version_20251023.txt"


## data from olink downloaded from https://insight.olink.com/plan-study/panel-selection/explore
dir_olink <- "/Users/c.giambartolomei/Downloads/annotations_of_proteins_oct222025/olink_downloadOct202015/"

## PATH WHERE INTERVAL AND CHRIS LOD AND CV ARE
data_path <- "/Users/c.giambartolomei/Downloads/annotations_of_proteins_oct222025/Fig1_Solene_Meta_paper-main/Results_files_for_plots"

# PATH WHERE WE HAVE BLOOD CONCENTRATIONS AND RNA TISSUE INFO
hpa_path = "/Users/c.giambartolomei/Downloads/annotations_of_proteins_oct222025/uniprot_tissue/"


################################################################################
################################################################################
# START FROM MAPPING FILE
################################################################################
################################################################################
## --- [1] Load mapping file ---

mapping <- fread(mapping_file)
mapping$target = paste("seq.", gsub("-", ".", mapping$SeqId), sep="")
# Rename to avoid conflicts
mapping <- mapping %>%
  dplyr::rename(SeqId_original = SeqId, SeqId = target)

################################################################################
################################################################################
### SECTION 0: Add CV and LOD of each study
################################################################################
################################################################################
setwd(data_path)
# ---- 1. Load INTERVAL and CHRIS CV/LOD ----
LOD_int <- readRDS("LOD_in_interval.Rds")
CV_int  <- readRDS("CV_in_interval.Rds")
LOD_chr <- readRDS("LOd_in_chris.Rds")
CV_chr  <- readRDS("CV_in_chris.Rds")


# ---- 2. Filter to human proteins ----
LOD_int_HP <- LOD_int %>% filter(Type == "Protein", Organism %in% c("Human"))
CV_int_HP  <- CV_int  %>% filter(Type.y == "Protein", Organism.y %in% c("Human"))

LOD_chr_HP <- LOD_chr %>% filter(Type == "Protein", Organism %in% c("Human"))
CV_chr_HP  <- CV_chr  %>% filter(Type.y == "Protein", Organism.y %in% c("Human"))

# ---- Join and Clean INTERVAL ----
#    Samples_under_LOD_batch_1 = Samples_under_LOD_batch_1,
#    Samples_under_LOD_batch_2 = Samples_under_LOD_batch_2,

interval_merged <- LOD_int_HP %>%
  dplyr::select(
    SeqId = seq,
    LOD_batch1 = LOD1,
    LOD_batch2 = LOD2,
    Percent_under_LOD_batch1 = Percent_under_LOD_batch_1,
    Percent_under_LOD_batch2 = Percent_under_LOD_batch_2,
    Dilution
  ) %>%
  left_join(
    CV_int_HP %>% dplyr::select(SeqId = targets, CV_batch1 = CV_QC1, CV_batch2 = CV_QC2),
    by = "SeqId"
  ) %>%
  rename_with(~ paste0("INTERVAL_", .x), .cols = -SeqId)

# ---- Join and Clean  CHRIS ----
#     Samples_under_LOD = Samples_under_LOD,
chris_merged <- LOD_chr_HP %>%
  dplyr::select(
    SeqId = seq,
    LOD = LOD1,
    Percent_under_LOD = Percent_under_LOD,
    Dilution
  ) %>%
  left_join(
    CV_chr_HP %>% dplyr::select(SeqId = targets, CV = CV_QC),
    by = "SeqId"
  ) %>%
  rename_with(~ paste0("CHRIS_", .x), .cols = -SeqId)

# ---- 3. Merge INTERVAL and CHRIS ----
lod_cv_combined <- full_join(interval_merged, chris_merged, by = "SeqId")


# ---- Add QC to mapping ----
mapping <- mapping %>%
  left_join(lod_cv_combined, by = "SeqId")

n_seqid <- dplyr::n_distinct(mapping$SeqId)
n_uniprot <- dplyr::n_distinct(mapping$UniProt_ID)

cat("Unique SeqIds in final file:", n_seqid, "\n")
cat("Unique UniProt IDs in final file:", n_uniprot, "\n")

# Unique SeqIds in final file: 7289 
# Unique UniProt IDs in final file: 6381 


# print("Whatch for NA values interpreted as passing thresholds: convert to 0")
print("Number of rows in original mapping file:")
print(nrow(mapping))
# ---- 1. Replace NA values with 0 in the relevant percent LOD columns ----
#mapping <- mapping %>%
#  mutate(
#    INTERVAL_Samples_under_LOD_batch_1 = ifelse(is.na(INTERVAL_Samples_under_LOD_batch_1), 0, INTERVAL_Samples_under_LOD_batch_1),
#    INTERVAL_Samples_under_LOD_batch_2 = ifelse(is.na(INTERVAL_Samples_under_LOD_batch_2), 0, INTERVAL_Samples_under_LOD_batch_2),
#    CHRIS_Samples_under_LOD            = ifelse(is.na(CHRIS_Samples_under_LOD), 0, CHRIS_Samples_under_LOD)
#  )
# ---- 2. Replace NA values in sample-under-LOD columns with 0
replace = FALSE # I set this as TRUE in the version of the Supp Table 1 befoe Nov 7
if (replace) {
mapping <- mapping %>%
                mutate(
                INTERVAL_Percent_under_LOD_batch1 = ifelse(is.na(INTERVAL_Percent_under_LOD_batch1), 0, INTERVAL_Percent_under_LOD_batch1),
                INTERVAL_Percent_under_LOD_batch2 = ifelse(is.na(INTERVAL_Percent_under_LOD_batch2), 0, INTERVAL_Percent_under_LOD_batch2),
                CHRIS_Percent_under_LOD = ifelse(is.na(CHRIS_Percent_under_LOD), 0, CHRIS_Percent_under_LOD)
        )
}

# ---- 3. Create the pass_qc column ----
# TRUE if the SeqId passes both INTERVAL and CHRIS thresholds:
# - INTERVAL: both batch1 and batch2 must be < 50%
# - CHRIS: must be < 80%
mapping <- mapping %>%
  mutate(
    pass_qc = (INTERVAL_Percent_under_LOD_batch1 < 50) &
              (INTERVAL_Percent_under_LOD_batch2 < 50) &
              (CHRIS_Percent_under_LOD < 80)
  )

if (!replace) { 
mapping <- mapping %>% 
  mutate( 
    pass_qc = ( 
      (coalesce(INTERVAL_Percent_under_LOD_batch1, 0) < 50) & 
      (coalesce(INTERVAL_Percent_under_LOD_batch2, 0) < 50) & 
      (coalesce(CHRIS_Percent_under_LOD, 0) < 80) 
    ) 
  ) 
} 

# ---- 4. Count how many pass ----
num_pass <- sum(mapping$pass_qc, na.rm = TRUE)
num_unique_pass <- length(unique(mapping$SeqId[mapping$pass_qc == TRUE]))

print(paste("Number of rows passing QC:", num_pass))
print(paste("Number of unique SeqIds passing QC:", num_unique_pass))

print("Number of rows in final mapping file:")
print(nrow(mapping))

################################################################################
################################################################################
### SECTION 1: Add historical assays
################################################################################
################################################################################

## --- [2] Get SomaScan v4.0 panel aptamers (≈5k) ---
seqids_v40 <- somascan_menu[["v4.0"]]
seqids_v40_fmt <- paste0("seq.", gsub("-", ".", seqids_v40))
mapping$aptamer_in_somascan5k <- mapping$SeqId %in% seqids_v40_fmt

## --- [3] Load Olink HT + Explore 3072 UniProt sets ---
ht_file <- fread(file.path(dir_olink, "5416-proteins_HT.csv"), sep = ";")
ht_file_expanded <- ht_file %>%
  filter(!is.na(`UniProt ID`)) %>%
  separate_rows(`UniProt ID`, sep = "_") 

# Eight constituent Explore 3072 panel files
panel_files <- list.files(dir_olink, pattern = "proteins_.*\\.csv", full.names = TRUE)
panel_files <- panel_files[!str_detect(panel_files, "HT|Reveal")]

# Load and combine Explore 3072 panels
explore_3072 <- rbindlist(lapply(panel_files, fread, sep = ";"))
explore_3072_expanded <- explore_3072 %>%
  filter(!is.na(`UniProt ID`)) %>%
  separate_rows(`UniProt ID`, sep = "_")

# Unique protein sets
u_3072 <- unique(explore_3072_expanded$`UniProt ID`)
u_ht <- unique(ht_file_expanded$`UniProt ID`)

## --- [4] Add columns: UniProt seen in Olink panels? ---
mapping$uniprot_in_ExploreHT <- mapping$UniProt_ID %in% u_ht
mapping$uniprot_in_Explore3072 <- mapping$UniProt_ID %in% u_3072


# But caution if I use this: If a protein has multiple aptamers and only one was newly added, it may be flagged as “new” even if it was measured before (less strict definition).
#mapping <- mapping %>%
#  group_by(UniProt_ID) %>%
#  mutate(previously_in_somascan_v4 = all(aptamer_in_somascan5k == TRUE),
#         new_in_somascan_v7 = !previously_in_somascan_v4) %>%
#  ungroup()

# define as: protein as “present in v5 (5k)” if any of its aptamers were present in 5k, “new to v7” if none of its aptamers were present in v5.
# avoids mixed classification due to multiple aptamers mapping to one target.
# uniprot_present_in_v5	any(aptamer_in_somascan5k == TRUE)	Protein was measured in the 5k assay
# uniprot_new_in_somascan7k_vs5k	all(aptamer_in_somascan5k == FALSE)	Protein was not measured in 5k; newly added in v7
mapping <- mapping %>%
  group_by(UniProt_ID) %>%
  mutate(
    # Strict definition:
    # A protein is considered present in v5 if ANY of its aptamers were in the v4/v5 panel
    uniprot_in_somascan5k = any(aptamer_in_somascan5k == TRUE, na.rm = TRUE),
    # A protein is considered truly new in v7 only if NONE of its aptamers were in v5
    uniprot_new_in_somascan7k_vs5k = !uniprot_in_somascan5k
  ) %>%
  ungroup()


## --- [5] Flag proteins previously seen in *any* platform
# Any platform = SomaScan v4.0 (by SeqId) OR Olink HT/3072 (by UniProt)
mapping <- mapping %>%
  group_by(UniProt_ID) %>%
  mutate(uniprot_in_any_platform = any(
    uniprot_in_somascan5k == TRUE |
    uniprot_in_ExploreHT == TRUE |
    uniprot_in_Explore3072 == TRUE,
    na.rm = TRUE
  )) %>%
  ungroup()


## --- [6] Optional: write output
# fwrite(mapping, "somascan_mapping_annotated_Oct22_2025.txt", sep = "\t")



################################################################################
################################################################################
## SECTION 2: Add Annotation of plasma proteins with secretion pathway: files and snippet taken from Solene Cadiou 
# Solene script: Fig1_Solene_Meta_paper-main/202411_figure_panel_old_vs_new.R
# Solene data as listed and under Results_files_for_plots
################################################################################
################################################################################

##prepare location attribution
# df_uniprot <- fread("/exchange/healthds/pQTL/pQTL_workplace/protein_annotation/uniprotkb_AND_reviewed_true_AND_model_o_2024_09_09.tsv", header = TRUE, sep = "\t")
df_uniprot <- fread(paste(data_path, "/uniprotkb_AND_reviewed_true_AND_model_o_2024_09_09.tsv", sep=""), header = TRUE, sep = "\t")

# Classify proteins based on UniprotDB keywords for the blood proteome:
# 1. Secreted proteins (KW-0964)
# 2. Membrane (KW-0472)
# 3. Intracellular (the rest)

# df_uniprot<-df_uniprot[,c( "Entry",      "Entry Name"     ,  "Protein names"    ,"Gene Names"   ,       "Keyword ID"  ,     "Protein families")]
df_uniprot<-df_uniprot[,c( "Entry", "Keyword ID")]
# Function to annotate secretion type
location_annotation <- function(mapping_df, uniprot_df) {
  # Merge only on UniProt_ID ↔ Entry
  annotated <- mapping_df %>%
    left_join(uniprot_df, by = c("UniProt_ID" = "Entry")) %>%
    mutate(
      label_cel_loc = case_when(
        is.na(`Keyword ID`) ~ NA_character_,
        str_detect(`Keyword ID`, "KW-0964") ~ "secreted",
        str_detect(`Keyword ID`, "KW-0472") ~ "membrane",
        TRUE ~ "intracellular"
      )
    ) %>%
    dplyr::select(-`Keyword ID`)  # Drop UniProt column after classification
  return(annotated)
}

mapping = location_annotation(mapping, df_uniprot)
# Check how many proteins are unclassified (no UniProt match)
unmatched <- mapping[is.na(mapping$label_cel_loc), .N]
cat("Number of UniProt IDs not found in UniProt file:", nrow(unmatched), "\n")

################################################################################
################################################################################
## SECTION 3: RNA tissue annotations (Human Protein Atlas).
## SECTION 4: Blood concentrations (Human Protein Atlas). Script in /Users/c.giambartolomei/Downloads/annotations_of_proteins_oct222025/uniprot_tissue/hpa_tissue_v3.R 
################################################################################
################################################################################
# Giulia did this part so I do not need to replicate: 
add_HGNC_Symbol = FALSE
if (add_HGNC_Symbol) {
	# This section needs gene symbol
	library(org.Hs.eg.db) # Version should be 3.20.0 per Methods
	library(AnnotationDbi) # Version should be 1.68.0 per Methods
	packageVersion("AnnotationDbi") # ‘1.66.0’
	packageVersion("org.Hs.eg.db") # ‘3.19.1’
 
	# Ensure Entrez_Gene_ID is character (required for mapIds)
	mapping$Entrez_Gene_ID <- as.character(mapping$Entrez_Gene_ID)
 
	# Map Entrez IDs to gene symbols
	gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = mapping$Entrez_Gene_ID,
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")

	# Add gene symbol to mapping
	mapping$HGNC_Symbol <- gene_symbols

}

###################################################################################
setwd(hpa_path)

# ---- 0) Load input files ----
mydata <- mapping
cat_raw   <- read_tsv("tissue_category_rna_Tissue.tsv", show_col_types = FALSE)
consensus <- read_tsv("rna_tissue_consensus.tsv", show_col_types = FALSE)
blood     <- read_tsv("blood_ms_concentration.tsv", show_col_types = FALSE)
map       <- read_tsv("proteinatlas.tsv", show_col_types = FALSE)

# ---- 1) Prepare HPA RNA tissue annotation ----
cat_clean <- cat_raw %>%
  mutate(
    Uniprot = gsub("\\s","", ifelse(is.na(Uniprot),"",Uniprot)),
    Gene    = toupper(ifelse(is.na(Gene),"",Gene))
  ) %>%
  separate_rows(Uniprot, sep = ",") %>%
  filter(Uniprot != "") %>%
  mutate(Gene_is_ENSG = str_detect(Gene, "^ENSG"),
         Gene_symbol  = if_else(Gene_is_ENSG, NA_character_, Gene)) %>%
  filter(!is.na(Gene_symbol)) %>%
  dplyr::select(Uniprot, Gene_symbol,
         `RNA tissue specificity`,
         `RNA tissue distribution`,
         `RNA tissue specific nTPM`) %>%
  distinct()

# Gene-level mapping
gene_to_cat <- cat_clean %>%
  distinct(Gene_symbol,
           `RNA tissue specificity`,
           `RNA tissue distribution`,
           `RNA tissue specific nTPM`)

# UniProt → all HPA genes
uniprot_to_gene <- cat_clean %>%
  group_by(Uniprot) %>%
  summarise(HPA_genes = paste(sort(unique(Gene_symbol)), collapse = "; "), .groups = "drop")

# ---- 2) RNA consensus: tissues detected at nTPM ≥1 ----
detected_by_gene <- consensus %>%
  mutate(GENE_UP = toupper(`Gene name`)) %>%
  filter(nTPM >= 1) %>%
  group_by(GENE_UP) %>%
  summarise(`RNA tissues detected (nTPM>=1)` = paste(sort(unique(Tissue)), collapse = "; "), .groups = "drop")

# ---- 3) Blood concentration mapping ----
submap <- map %>% dplyr::select(Gene, Uniprot, Ensembl)

blood_map_full <- blood %>%
  left_join(submap, by = c("ENSG ID" = "Ensembl"))

# Check that each UniProt has unique concentration
dup_check <- blood_map_full %>%
  filter(!is.na(Uniprot)) %>%
  group_by(Uniprot) %>%
  summarise(n_conc = n_distinct(`Conc [pg/L]`, na.rm = TRUE),
            conc_vals = paste(sort(unique(`Conc [pg/L]`)), collapse = "; "),
            .groups = "drop")

if (any(dup_check$n_conc > 1)) {
  stop("Some UniProt IDs have conflicting concentrations")
}

blood_map <- blood_map_full %>%
  arrange(Uniprot) %>%
  distinct(Uniprot, .keep_all = TRUE) %>%
  dplyr::select(
    Uniprot,
    Blood_gene_symbol = Gene.x,
    Blood_conc_pgL    = `Conc [pg/L]`
  )

# ---- 4) Combine everything ----
mydata_up <- mydata %>%
  mutate(HGN_UP = toupper(HGNC_Symbol)) %>%
  left_join(gene_to_cat,      by = c("HGN_UP" = "Gene_symbol")) %>%
  left_join(uniprot_to_gene,  by = c("UniProt_ID" = "Uniprot")) %>%
  left_join(detected_by_gene, by = c("HGN_UP" = "GENE_UP")) %>%
  left_join(blood_map,        by = c("UniProt_ID" = "Uniprot")) %>%
  mutate(
    `RNA tissue specific NX` = `RNA tissue specific nTPM`,

    # fixed: per-row check
    HPA_mapping_OK = if_else(
      is.na(HPA_genes), NA,
      mapply(function(g, hs) g %in% strsplit(hs, ";\\s*")[[1]], HGN_UP, HPA_genes)
    ),

    Blood_mapping_OK = case_when(
      is.na(Blood_gene_symbol) ~ NA,
      HGN_UP == toupper(Blood_gene_symbol) ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  dplyr::select(
    everything(),
    `RNA tissue specificity`,
    `RNA tissue distribution`,
    `RNA tissues detected (nTPM>=1)`,
    `RNA tissue specific NX`,
    Blood_conc_pgL,
    HPA_mapping_OK,
    Blood_mapping_OK
  ) %>%
  dplyr::select(-c(`RNA tissue specific nTPM`, HGN_UP)) %>%
  mutate(across(where(is.character), ~ na_if(., "")))




# ---- 5) Summary counts ----
n_total_uni  <- n_distinct(mydata_up$UniProt_ID)
n_with_blood <- n_distinct(mydata_up$UniProt_ID[!is.na(mydata_up$Blood_conc_pgL)])

n_hpa_true  <- n_distinct(mydata_up$UniProt_ID[mydata_up$HPA_mapping_OK %in% TRUE])
n_hpa_false <- n_distinct(mydata_up$UniProt_ID[mydata_up$HPA_mapping_OK %in% FALSE])
n_hpa_na    <- n_distinct(mydata_up$UniProt_ID[is.na(mydata_up$HPA_mapping_OK)])

n_blood_true  <- n_distinct(mydata_up$UniProt_ID[mydata_up$Blood_mapping_OK %in% TRUE])
n_blood_false <- n_distinct(mydata_up$UniProt_ID[mydata_up$Blood_mapping_OK %in% FALSE])
n_blood_na    <- n_distinct(mydata_up$UniProt_ID[is.na(mydata_up$Blood_mapping_OK)])

summary_counts <- tibble::tibble(
  metric = c(
    "Unique UniProt_IDs",
    "With blood concentration (pg/L)",
    "HPA mapping OK (TRUE)",
    "HPA mapping mismatch (FALSE)",
    "HPA mapping missing (NA)",
    "Blood mapping OK (TRUE)",
    "Blood mapping mismatch (FALSE)",
    "Blood mapping missing (NA)"
  ),
  count = c(
    n_total_uni,
    n_with_blood,
    n_hpa_true,
    n_hpa_false,
    n_hpa_na,
    n_blood_true,
    n_blood_false,
    n_blood_na
  )
)

# for the full mapping fle with 7,289 seqids:
#                           metric count
#              Unique UniProt_IDs  6381
# With blood concentration (pg/L)  2664
#           HPA mapping OK (TRUE)  6308
#    HPA mapping mismatch (FALSE)    16
#        HPA mapping missing (NA)    60
#         Blood mapping OK (TRUE)  2663
#  Blood mapping mismatch (FALSE)     4
#      Blood mapping missing (NA)  3717


mapping = mydata_up 

################################################################################
################################################################################
# END and write
################################################################################
################################################################################
rearrange = TRUE
if (rearrange) {
	# Move the HGNC_Symbol column near the front (e.g., right after Entrez_Gene_ID)
	mapping <- mapping %>%
  		dplyr::relocate(HGNC_Symbol, .after = Entrez_Gene_ID)
	mapping <- mapping %>%
  		dplyr::relocate(Ensembl_ID, .after = Entrez_Gene_ID)
	mapping <- mapping %>%
  		dplyr::relocate(SeqId, .after = SeqId_original)
}


cols_to_remove <- c("SeqId_original")
mapping <- mapping %>%
  dplyr::select(-dplyr::any_of(cols_to_remove))

remove = FALSE 
if (remove) {
# Remove unwanted columns if present

 cols_to_remove <- c(
  "SomaScan_UniProt_ID",
  "SomaScan_Entrez_Gene_ID",
  "feature",
  "class",
  "assembly_unit",
  "UniProt_EntrezID_match",
  "HPA_mapping_OK", "Blood_mapping_OK"

)

mapping <- mapping %>%
  dplyr::select(-dplyr::any_of(cols_to_remove))

}

mapping <- mapping %>% dplyr::rename(Secretion.pathway = label_cel_loc)


write.xlsx(mapping, file = paste(data_path, "/somascan_tss_ncbi_grch37_version_20250326_eva_version_20251023_addedgenesbygiulia_annotated_v2.xlsx", sep=""))


n_seqid <- dplyr::n_distinct(mydata_up$SeqId)
n_uniprot <- dplyr::n_distinct(mydata_up$UniProt_ID)

cat("Unique SeqIds in final file:", n_seqid, "\n")
cat("Unique UniProt IDs in final file:", n_uniprot, "\n")

# Unique SeqIds in final file: 7289
# Unique UniProt IDs in final file: 6381




# now continue to plot: /Users/c.giambartolomei/Downloads/annotations_of_proteins_oct222025/uniprot_tissue/new_figure1_v5_2025_oct24.R
# need to:
# figure out what needs changes by putting it in bold with the new set-up? Also, I need to add a new panel F boxplot to the rigth demonstarting demonstarting this statement "more samples were below the LOD for aptamers targeting new proteins (median= ; IQR= ; P-value < 10-16 in all three batches) compared to previously assessed aptamers (median = ; IQR = ;)."

