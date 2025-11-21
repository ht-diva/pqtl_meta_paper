# Add column uniprot_new_in_somascan7k_vs5k to Supp Table 2
# Column legend: uniprot_new_in_somascan7k_vs5k
# Indicates whether the protein(s) represented by the UniProt ID(s) in each row are newly introduced in the SomaScan 7k (v4.1) assay compared with the SomaScan 5k (v4.0) panel, defined strictly at the UniProt level.
# TRUE – All UniProt ID(s) in the row are newly introduced in the SomaScan 7k assay; none of their corresponding aptamers were present in the SomaScan 5k (v4.0) panel.
# FALSE – At least one UniProt ID in the row was previously represented in the SomaScan 5k (v4.0) panel (i.e., the protein was already covered in the prior assay).
# For rows containing multiple UniProt IDs (collapsed entries separated by “|”), the value is assigned as FALSE if any UniProt ID in the group was previously present in the 5k panel, and TRUE only if all UniProt IDs were newly introduced in the 7k panel

# Load necessary packages
library(readxl)
library(dplyr)
library(stringr)
library(openxlsx)

# -----------------------------
# Step 1: Read in the data
# -----------------------------

# Read Supp Table 1 (sheet 2)
supp1 <- readxl::read_excel("supplementary_table_1.xlsx", sheet = 2)

# Read Supp Table 2 (sheet 2)
# Header in row 3, data starts at row 2 (skip first 2 rows)
supp2 <- readxl::read_excel("supplementary_table_2.xlsx", sheet = 2, skip = 2)


# Split UniProt_IDs into lists
supp2_split <- stringr::str_split(supp2$UniProt_ID, "\\|")

# Flatten and extract all unique UniProt IDs from Table 2
uniprot2_unique <- unique(unlist(supp2_split))

# Unique UniProt IDs from Table 1
uniprot1_unique <- unique(supp1$UniProt_ID)

table(uniprot2_unique %in% uniprot1_unique)

#TRUE 
#3930 

# -----------------------------
# Step 2: Prepare lookup table
# -----------------------------

if(!("UniProt_ID" %in% colnames(supp1))) stop("Missing 'UniProt_ID' in supp1.")
if(!("uniprot_new_in_somascan7k_vs5k" %in% colnames(supp1))) stop("Missing 'uniprot_new_in_somascan7k_vs5k' in supp1.")

lookup <- supp1 %>%
  dplyr::select(UniProt_ID, uniprot_new_in_somascan7k_vs5k) %>%
  dplyr::distinct()

# -----------------------------
# Step 3: Handle UniProt_IDs
# -----------------------------

# Count rows before
n_before <- nrow(supp2)

# Split UniProt_IDs into lists (temporary)
supp2_temp <- supp2 %>%
  dplyr::mutate(UniProt_ID_split = stringr::str_split(UniProt_ID, "\\|"))

# Count rows with multiple UniProt IDs
multi_count <- sum(sapply(supp2_temp$UniProt_ID_split, length) > 1)
cat("Number of rows with multiple UniProt_IDs:", multi_count, "\n")

# -----------------------------
# Step 4: Compute 7k vs 5k status
# -----------------------------

get_coverage_status <- function(ids, lookup_df) {
  ids <- ids[ids != "" & !is.na(ids)]  # clean
  matches <- lookup_df %>%
    dplyr::filter(UniProt_ID %in% ids)
  
  if (nrow(matches) == 0) return(NA)
  
  # If ANY UniProt_ID was previously present (FALSE) → FALSE
  # Only TRUE if all UniProt_IDs were new
  if (any(matches$uniprot_new_in_somascan7k_vs5k == FALSE, na.rm = TRUE)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

# Apply row-wise to temporary df
supp2_temp <- supp2_temp %>%
  dplyr::rowwise() %>%
  dplyr::mutate(uniprot_new_in_somascan7k_vs5k = get_coverage_status(UniProt_ID_split, lookup)) %>%
  dplyr::ungroup()

# -----------------------------
# Step 5: Clean up (remove temporary column)
# -----------------------------

supp2_final <- supp2_temp %>%
  dplyr::select(-UniProt_ID_split)  # remove temporary list column

# -----------------------------
# Step 6: Count and verify
# -----------------------------

n_after <- nrow(supp2_final)
cat("Row count before:", n_before, "\n")
cat("Row count after :", n_after, "\n")

if (n_before == n_after) {
  cat("✅ Row counts match — no rows were lost or added.\n")
} else {
  cat("⚠️ Warning: row count mismatch!\n")
}

table(is.na(supp2_final$uniprot_new_in_somascan7k_vs5k))

#FALSE 
# 7870 
# Summary of TRUE/FALSE
cat("Summary of new column:\n")
print(table(supp2_final$uniprot_new_in_somascan7k_vs5k, useNA = "ifany"))

# -----------------------------
# Step 7: Write output
# -----------------------------

openxlsx::write.xlsx(supp2_final, "supplementary_table_2_with_uniprot_status.xlsx")

cat("✅ Done! Added 'uniprot_new_in_somascan7k_vs5k' column and saved as 'supplementary_table_2_with_uniprot_status.xlsx'.\n")

