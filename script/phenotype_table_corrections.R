###############################################
## Setup
###############################################
 
# Load required packages
library(readxl)
library(dplyr)
library(stringr)
library(writexl)   # only needed if you want to write out updated tables

###############################################
## 1. Read Excel file and basic checks
###############################################
 
input_file  <- "supplementary_tables_natmed_v1.xlsx"
sheet_index <- 2
 
# Read second sheet, skipping first 2 rows so that the 3rd row becomes the header
t <- read_excel(input_file, sheet = sheet_index, skip = 2)
 
# Basic checks
cat("First rows of phenotype table (after skipping first 2 rows):\n")
print(head(t))
 
cat("\nNumber of rows in table: ", nrow(t), "\n", sep = "")
cat("Number of unique Phenotype: ", dplyr::n_distinct(t$Phenotype), "\n\n", sep = "")
 
###############################################
## 2. Apply text substitutions to MVP_description immediately
###############################################
 
# Make sure MVP_description exists
if (!"MVP_description" %in% names(t)) {
  t <- t %>% mutate(MVP_description = NA_character_)
}
 
# Text substitutions
t <- t %>%
  mutate(
    MVP_description =
      MVP_description %>%
      # high density lipoprotein cholestrol measurement → cholesterol
      str_replace_all(
        regex("high density lipoprotein cholestrol measurement", ignore_case = TRUE),
        "high density lipoprotein cholesterol measurement"
      ) %>%
      # Creat → Creatinine
      str_replace_all(
        regex("\\bCreat\\b"),
        "Creatinine"
      ) %>%
      # Althete's foot → Athlete's foot
      str_replace_all(
        regex("Althete['’]s foot", ignore_case = TRUE),
        "Athlete's foot"
      ) %>%
      # NEC → not elsewhere classified
      str_replace_all(
        regex("\\bNEC\\b", ignore_case = TRUE),
        "not elsewhere classified"
      ) %>%
      # NOS → not otherwise specified
      str_replace_all(
        regex("\\bNOS\\b", ignore_case = TRUE),
        "not otherwise specified"
      ) %>%
      # glomerular filtration rate → estimated glomerular filtration rate
      str_replace_all(
        regex("(?<!estimated )glomerular filtration rate", ignore_case = TRUE),
        "estimated glomerular filtration rate"
      ) %>%
      # Migrain with aura → Migraine with aura
      str_replace_all(
        regex("Migrain with aura", ignore_case = TRUE),
        "Migraine with aura"
      ) %>%
      # SHBG → sex hormone-binding globulin measurement
      str_replace_all(
        regex("\\bSHBG\\b", ignore_case = TRUE),
        "sex hormone-binding globulin measurement"
      ) %>%
      # Sepsis and SIRS → full expansion
      str_replace_all(
        regex("Sepsis and SIRS", ignore_case = TRUE),
        "Sepsis and Systemic Inflammatory Response Syndrome (SIRS)"
      ) %>%
      # BSP → bone sialoprotein measurement
      str_replace_all(
        regex("\\bBSP\\b", ignore_case = TRUE),
        "bone sialoprotein measurement"
      )
  )
 
###############################################
## 3. Add snippet-based columns (Data_source, PHENOTYPE_CLASS)
###############################################
 
t_summary <- t %>%
  mutate(
  Data_source = case_when(
    (MVP_Cases + MVP_Controls > 0) & (UKBB_Cases + UKBB_Controls > 0) & (FinnGen_Cases + FinnGen_Controls > 0) ~ "MVP+UKB+FinnGen",
    (MVP_Cases + MVP_Controls > 0) & (UKBB_Cases + UKBB_Controls > 0) ~ "MVP+UKB",
    (MVP_Cases + MVP_Controls > 0) & (FinnGen_Cases + FinnGen_Controls > 0) ~ "MVP+FinnGen",
    (UKBB_Cases + UKBB_Controls > 0) & (FinnGen_Cases + FinnGen_Controls > 0) ~ "UKB+FinnGen",
    (MVP_Cases + MVP_Controls > 0) ~ "MVP",
    (UKBB_Cases + UKBB_Controls > 0) ~ "UKB",
    (FinnGen_Cases + FinnGen_Controls > 0) ~ "FinnGen",
    TRUE ~ NA_character_
),
 
    # Initial phenotype class based on string matching
    PHENOTYPE_CLASS = case_when(
      str_detect(Phenotype, "^phe") ~ "diseases",
      str_detect(Phenotype, "Med")  ~ "medications",
      TRUE                          ~ "traits"
    )
  )
 

###############################################
## 4. Define PHENOTYPE_CLASS_new from case-control logic
###############################################
 
t_summary <- t_summary %>%
  mutate(
    # Identifies true disease phenotypes
    has_mvp_case_control = MVP_Cases > 0 & MVP_Controls > 0,
    has_ukb_case_control = UKBB_Cases > 0 & UKBB_Controls > 0,
    has_finngen_case_control = FinnGen_Cases > 0 & FinnGen_Controls > 0,
 
    PHENOTYPE_CLASS_new = case_when(
      str_detect(Phenotype, "Med") ~ "medications",
      has_mvp_case_control | has_ukb_case_control | has_finngen_case_control ~ "diseases",
      TRUE ~ "traits"
    ),
 
    PHENOTYPE_CLASS_mismatch = if_else(
      dplyr::coalesce(PHENOTYPE_CLASS, "NA") != dplyr::coalesce(PHENOTYPE_CLASS_new, "NA"),
      1L, 0L
    )
  ) %>%
  select(-has_mvp_case_control, -has_ukb_case_control, -has_finngen_case_control)
 
###############################################
## 5. Check mismatches
###############################################

mismatches <- t_summary %>%
  filter(PHENOTYPE_CLASS_mismatch == 1L)
 
cat("Number of mismatches between PHENOTYPE_CLASS and PHENOTYPE_CLASS_new: ",
    nrow(mismatches), "\n", sep = "")
 
if (nrow(mismatches) > 0) {
  cat("First mismatches:\n")
  print(head(mismatches %>%
               select(Phenotype, PHENOTYPE_CLASS, PHENOTYPE_CLASS_new, MVP_Cases, MVP_Controls,
                      UKBB_Cases, UKBB_Controls, FinnGen_Cases, FinnGen_Controls)))
}
 
###############################################

## 6. Save updated table
###############################################
 
output_file <- "MVP_phenotypes_updated_with_substitutions.xlsx"
 
write_xlsx(
  list(
    "Phenotypes_cleaned" = t_summary
  ),
  path = output_file
)
 
cat("Updated table written to:", output_file, "\n")

