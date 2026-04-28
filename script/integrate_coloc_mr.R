# Date: March 30 2026
# Claudia Giambartolomei

library(data.table)
library(readxl)

### set WD
setwd("/Users/c.giambartolomei/Downloads/")

### ------------------------------------------------------------
### READ DATA 
### ------------------------------------------------------------

### read coloc
# coloc1 <- fread("INTERVAL_CHRIS_coloc_final_result.txt")
coloc <- fread("/Users/c.giambartolomei/Downloads/INTERVAL_CHRIS_coloc_rerun_final_result.txt")
cat("Original coloc rows:", nrow(coloc), "\n")

### save exact duplicate rows
n_full_dups <- nrow(coloc) - uniqueN(coloc)
cat("Number of fully duplicated rows:", n_full_dups, "\n")

coloc_full_dups <- coloc[duplicated(coloc) | duplicated(coloc, fromLast = TRUE)]
# fwrite(coloc_full_dups, "coloc_full_duplicate_rows.txt", sep = "\t")

### remove exact duplicate rows
coloc <- unique(coloc)
cat("After removing full duplicates:", nrow(coloc), "\n")

### read mr (sheet ST7)
mr <- as.data.table(read_excel("supplementary_table_7.xlsx", sheet = "ST7"))

### rename coloc columns to match mr
setnames(coloc, "gene_name", "HARMONIZED_GENE_NAME")
setnames(coloc, "seqID", "SeqID")
setnames(coloc, "phenotype", "PHENOTYPE")
setnames(coloc, "protein", "UniProt_ID")

### ------------------------------------------------------------
### columns are not clean e.g. coloc data contains  whitespaces like  "KIR2DL5A " so need to clean first
### ------------------------------------------------------------

### cleaning function for character matching
clean_chr <- function(x) {
  x <- as.character(x)
  x <- gsub("\u00A0", " ", x, fixed = TRUE)   # non-breaking spaces
  x <- gsub("[[:space:]]+", " ", x)           # collapse repeated spaces/tabs
  x <- trimws(x)                              # trim ends
  x[x == ""] <- NA_character_
  x
}

cols_to_clean_coloc <- c("HARMONIZED_GENE_NAME", "SeqID", "PHENOTYPE", "UniProt_ID")

coloc[, (cols_to_clean_coloc) := lapply(.SD, clean_chr), .SDcols = cols_to_clean_coloc]

### make sure numeric coloc columns are numeric
num_cols_coloc <- intersect(
  c("GWAS_beta", "GWAS_se", "GWAS_Pvalue",
    "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf", "nsnps"),
  names(coloc)
)

coloc[, (num_cols_coloc) := lapply(.SD, as.numeric), .SDcols = num_cols_coloc]

cols_to_clean_mr    <- c("HARMONIZED_GENE_NAME", "SeqID", "PHENOTYPE", "UniProt_ID", "SNPID")
mr[,    (cols_to_clean_mr)    := lapply(.SD, clean_chr), .SDcols = cols_to_clean_mr]

### ------------------------------------------------------------
### check relevant coloc entries are present in MR entries
### ------------------------------------------------------------

## GENE
table(coloc$HARMONIZED_GENE_NAME %in% mr$HARMONIZED_GENE_NAME)
table(coloc$SeqID %in% mr$SeqID)
table(coloc$UniProt_ID %in% mr$UniProt_ID)
table(coloc$PHENOTYPE %in% mr$PHENOTYPE)

# phenotype needs some modifications to match...
## mismatching coloc phenotypes
coloc_not_in_mr <- setdiff(unique(coloc$PHENOTYPE),
				unique(mr$PHENOTYPE))

# coloc_not_in_mr
# [1] "max_BMI"            "max_BMI_INT"        "max_Diastolic_INT" 
# [4] "max_Pulse"          "max_Pulse_INT"      "max_Systolic"      
# [7] "max_Systolic_INT"   "max_Weight"         "max_Weight_INT"    
#[10] "mean_BMI"           "mean_BMI_INT"       "mean_Diastolic"    
#[13] "mean_Diastolic_INT" "mean_Pulse"         "mean_Pulse_INT"    
#[16] "mean_Systolic"      "mean_Systolic_INT"  "mean_Weight"       
#[19] "mean_Weight_INT"    "min_BMI"            "min_BMI_INT"       
#[22] "min_Diastolic"      "min_Diastolic_INT"  "min_Pulse"         
#[25] "min_Pulse_INT"      "min_Systolic"       "min_Weight"        
#[28] "min_Weight_INT"    

# recode the few mismatching PHENOTYPE values
pheno_map <- c(
  "max_BMI" = "BMI_Max",
  "max_BMI_INT" = "BMI_Max_INT",
  "max_Diastolic_INT" = "Diastolic_Max_INT",
  "max_Pulse" = "Pulse_Max",
  "max_Pulse_INT" = "Pulse_Max_INT",
  "max_Systolic" = "Systolic_Max",
  "max_Systolic_INT" = "Systolic_Max_INT",
  "max_Weight" = "Weight_Max",
  "max_Weight_INT" = "Weight_Max_INT",
  "mean_BMI" = "BMI_Mean",
  "mean_BMI_INT" = "BMI_Mean_INT",
  "mean_Diastolic" = "Diastolic_Mean",
  "mean_Diastolic_INT" = "Diastolic_Mean_INT",
  "mean_Pulse" = "Pulse_Mean",
  "mean_Pulse_INT" = "Pulse_Mean_INT",
  "mean_Systolic" = "Systolic_Mean",
  "mean_Systolic_INT" = "Systolic_Mean_INT",
  "mean_Weight" = "Weight_Mean",
  "mean_Weight_INT" = "Weight_Mean_INT",
  "min_BMI" = "BMI_Min",
  "min_BMI_INT" = "BMI_Min_INT",
  "min_Diastolic" = "Diastolic_Min",
  "min_Diastolic_INT" = "Diastolic_Min_INT",
  "min_Pulse" = "Pulse_Min",
  "min_Pulse_INT" = "Pulse_Min_INT",
  "min_Systolic" = "Systolic_Min",
  "min_Weight" = "Weight_Min",
  "min_Weight_INT" = "Weight_Min_INT"
)

coloc[PHENOTYPE %in% names(pheno_map), PHENOTYPE := pheno_map[PHENOTYPE]]


### ------------------------------------------------------------
# Create mr_long with instrument order
### ------------------------------------------------------------
# now for the coloc it is by instrument, so need to make a long data frame for MR after splitting by instrument

## build MR long table from CHR + split POS_38
## create a stable MR row ID: I need this to build the coloc into the collapsed MR version again
mr[, mr_row_id := .I]

## split POS_38 in the original order
mr_long <- mr[
  ,
  {
    pos_vec <- strsplit(POS_38, "|", fixed = TRUE)[[1]]
    .(
      instrument = paste(CHR, pos_vec, sep = "_"),
      instrument_order = seq_along(pos_vec)
    )
  },
  by = .(mr_row_id)
]

## add all original MR columns back
mr_long <- mr[mr_long, on = "mr_row_id"]

table(coloc$instrument %in% mr_long$instrument)

## add all original MR columns back
mr_long <- mr[mr_long, on = "mr_row_id"]

### ------------------------------------------------------------
### merge long MR with coloc
### ------------------------------------------------------------
# merge long MR with coloc keeping all columns from both datasets
merge_cols <- c("HARMONIZED_GENE_NAME", "SeqID", "PHENOTYPE", "instrument")

mr_coloc_long <- merge(
  mr_long,
  coloc,
  by = merge_cols,
  all.x = TRUE,
  suffixes = c("_mr", "_coloc")
)

### ------------------------------------------------------------
### back to MR level summaries
### ------------------------------------------------------------

# Summarize back to MR-row level
pp_h4_threshold <- 0.8

collapse_in_order <- function(x) {
  paste(ifelse(is.na(x), NA_character_, as.character(x)), collapse = "|")
}

mr_summary <- mr_coloc_long[
  order(mr_row_id, instrument_order),
  {
    n_instr <- .N
    n_with_coloc <- sum(!is.na(PP.H4.abf))
    any_coloc <- n_with_coloc > 0

    cis_epi <- unique(cis_epitope_effect)

    max_PP.H4.abf <- if (any_coloc) max(PP.H4.abf, na.rm = TRUE) else NA_real_
    best_instrument <- if (any_coloc) instrument[which.max(fifelse(is.na(PP.H4.abf), -Inf, PP.H4.abf))] else NA_character_
    best_instrument_order <- if (any_coloc) which.max(fifelse(is.na(PP.H4.abf), -Inf, PP.H4.abf)) else NA_integer_

    GWAS_beta_by_instrument   <- collapse_in_order(GWAS_beta)
    GWAS_se_by_instrument     <- collapse_in_order(GWAS_se)
    GWAS_Pvalue_by_instrument <- collapse_in_order(GWAS_Pvalue)

    coloc_status <- if (!any_coloc && cis_epi == "TRUE") {
      "excluded_epitope_effect"
    } else if (!any_coloc) {
      "no_coloc_available"
    } else if (max_PP.H4.abf > pp_h4_threshold) {
      "colocalizing"
    } else {
      "not_colocalizing"
    }

    coloc_coverage <- if (n_with_coloc == 0) {
      "none"
    } else if (n_with_coloc == n_instr) {
      "complete"
    } else {
      "partial"
    }

    .(
      coloc_status = coloc_status,
      any_coloc_result = any_coloc,
      n_instruments = n_instr,
      n_instruments_with_coloc = n_with_coloc,
      coloc_coverage = coloc_coverage,

      max_PP.H4.abf = max_PP.H4.abf,
      best_instrument = best_instrument,
      best_instrument_order = best_instrument_order,

      GWAS_beta_by_instrument   = GWAS_beta_by_instrument,
      GWAS_se_by_instrument     = GWAS_se_by_instrument,
      GWAS_Pvalue_by_instrument = GWAS_Pvalue_by_instrument,

      PP.H0.abf_by_instrument = collapse_in_order(PP.H0.abf),
      PP.H1.abf_by_instrument = collapse_in_order(PP.H1.abf),
      PP.H2.abf_by_instrument = collapse_in_order(PP.H2.abf),
      PP.H3.abf_by_instrument = collapse_in_order(PP.H3.abf),
      PP.H4.abf_by_instrument = collapse_in_order(PP.H4.abf)
    )
  },
  by = mr_row_id
]
### -----------------------------
### attach summary back to original MR table
### -----------------------------
mr_final <- merge(
  mr,
  mr_summary,
  by = "mr_row_id",
  all.x = TRUE
)

### If there are MR rows completely absent from mr_coloc_long summary,
### assign them explicitly as no coloc available
mr_final[is.na(coloc_status) & cis_epitope_effect == "TRUE", coloc_status := "excluded_epitope_effect"]
mr_final[is.na(coloc_status) & cis_epitope_effect != "TRUE", coloc_status := "no_coloc_available"]
mr_final[is.na(any_coloc_result), any_coloc_result := FALSE]
mr_final[is.na(n_instruments), n_instruments := 0L]
mr_final[is.na(n_instruments_with_coloc), n_instruments_with_coloc := 0L]
mr_final[is.na(coloc_coverage), coloc_coverage := "none"]

mr_final[is.na(PP.H0.abf_by_instrument), PP.H0.abf_by_instrument := "NA"]
mr_final[is.na(PP.H1.abf_by_instrument), PP.H1.abf_by_instrument := "NA"]
mr_final[is.na(PP.H2.abf_by_instrument), PP.H2.abf_by_instrument := "NA"]
mr_final[is.na(PP.H3.abf_by_instrument), PP.H3.abf_by_instrument := "NA"]
mr_final[is.na(PP.H4.abf_by_instrument), PP.H4.abf_by_instrument := "NA"]


### -----------------------------
### Report how many MR rows have no corresponding coloc result
### -----------------------------
n_no_coloc_available <- mr_final[coloc_status == "no_coloc_available", .N]
cat("MR rows with no coloc available:", n_no_coloc_available, "\n")

### -----------------------------
## Rows with partial coloc availability, e.g. PP.H4.abf_by_instrument = "0.91|NA|0.12|0.85" -> These have some instruments with coloc, some without (n_instruments_with_coloc != n_instruments)
### -----------------------------
n_partial_coloc <- mr_final[coloc_coverage == "partial", .N]
cat("MR rows with partial coloc availability:", n_partial_coloc, "\n")

### -----------------------------
### save
### -----------------------------
# fwrite(mr_final, "mr_with_coloc_summary.txt", sep = "\t")
library(openxlsx)

wb <- createWorkbook()
addWorksheet(wb, "ST7_with_coloc")
writeData(wb, "ST7_with_coloc", mr_final)
saveWorkbook(wb, "supplementary_table_7_with_coloc_summary.xlsx", overwrite = TRUE)


cat("Total MR rows:", nrow(mr_final), "\n")
cat("Colocalizing rows:", mr_final[coloc_status == "colocalizing", .N], "\n")
cat("Non-colocalizing rows:", mr_final[coloc_status == "not_colocalizing", .N], "\n")
cat("No coloc available rows:", mr_final[coloc_status == "no_coloc_available", .N], "\n")
cat("Rows with complete coloc coverage:", mr_final[coloc_coverage == "complete", .N], "\n")
cat("Rows with partial coloc coverage:", mr_final[coloc_coverage == "partial", .N], "\n")
cat("Rows with no coloc coverage:", mr_final[coloc_coverage == "none", .N], "\n")



# 
#Total MR rows: 7,591
#3,454 (45.5%) colocalizing
#3,363 (44.3%) not colocalizing
#774 (10.2%) with no coloc available
#Coverage across instruments:
#6,456 (85.0%) complete
#361 (4.8%) partial
#774 (10.2%) none



# Quantify how similar PP values are across instruments within multi-instrument MR rows
library(data.table)

## helper: turn "0.91|NA|0.12|0.85" into numeric vector
split_num <- function(x) {
  vals <- strsplit(x, "|", fixed = TRUE)[[1]]
  vals[vals == "NA"] <- NA_character_
  as.numeric(vals)
}

## look at PP4 similarity across instruments
pp4_similarity <- mr_final[
  n_instruments > 1 & coloc_status != "excluded_epitope_effect",
  {
    pp4 <- split_num(PP.H4.abf_by_instrument)
    pp4 <- pp4[!is.na(pp4)]
    
    .(
      n_pp4 = length(pp4),
      pp4_mean = if (length(pp4) > 0) mean(pp4) else NA_real_,
      pp4_sd = if (length(pp4) > 1) sd(pp4) else NA_real_,
      pp4_range = if (length(pp4) > 0) max(pp4) - min(pp4) else NA_real_
    )
  },
  by = mr_row_id
]

summary(pp4_similarity$pp4_sd)
summary(pp4_similarity$pp4_range)

## e.g. proportion of multi-instrument rows with very similar PP4
mean(pp4_similarity$pp4_range < 0.05, na.rm = TRUE)
mean(pp4_similarity$pp4_range < 0.10, na.rm = TRUE)
mean(pp4_similarity$pp4_sd < 0.05, na.rm = TRUE)




#############
### Are multi-SNP rows are less likely to colocalize?
#############
## exclude rows not tested for coloc
tmp <- mr_final[coloc_status %in% c("colocalizing", "not_colocalizing")]

## define single vs multi SNP
tmp[, snp_group := ifelse(NUM_SNPS == 1, "single_SNP", "multi_SNP")]

## counts and fractions
tab <- tmp[
  ,
  .(
    n = .N,
    colocalizing = sum(coloc_status == "colocalizing"),
    not_colocalizing = sum(coloc_status == "not_colocalizing"),
    prop_colocalizing = round(mean(coloc_status == "colocalizing") * 100, 1)
  ),
  by = snp_group
]

print(tab)







# A couple of notes:
# 1. The posterior columns are returned as strings like: PP.H4.abf_by_instrument = "0.91|NA|0.12|0.85" in the same order as the original POS_37 instruments in the MR row.
# For each original MR row:
# 2. coloc_status == "colocalizing" if at least one instrument has coloc and at least one has PP.H4.abf > 0.8
#   coloc_status == "not_colocalizing" if there is at least one coloc result but none pass the threshold
#   coloc_status == "no_coloc_available" if no instrument in that MR row has any coloc result
# 3. any_coloc_result tells you whether there was any coloc row at all for that MR row.
# 4. n_instruments_with_coloc tells you how many instruments in that MR row found a coloc match.
