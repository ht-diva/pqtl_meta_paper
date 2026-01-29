# October 6, 2025
# 

# We aim to check if the cis signal on SDF2 (seq.15333.11) 
# is instead mediated by VTN, and it is actually a trans (correcting for VTN 
# value in the genomic association)

# Among the 463 ‘trans’ is an apparent cis-pQTL for SDF2, 
# which is nearby but seems more likely to be a trans-pQTL 
# driven by the rs704 variant - we could check this by doing 
# a mediation analysis in INTERVAL of the SDF2 ‘trans-pQTL’ 
# adjusting for VTN levels -> Dariush and Claudia: we need 
# to adjust for this variable in the GWAS model. 

library(pgenlibr)
library(data.table)
library(tidyverse)
library(Rmpfr) # to encounter 0 P-value


# Point to the .pgen/.pvar/.psam trio
path_pgen <- "/scratch/dariush.ghasemi/projects/pqtl_susie/results/test/fm/tmp/seq.13125.45_17_25283026_29975846_dosage.pgen"
path_pvar <- gsub(".pgen", ".pvar", path_pgen)
path_psam <- gsub(".pgen", ".psam", path_pgen)

# Phenotypes: proteins_residuals
path_pheno <- "/exchange/healthds/pQTL/INTERVAL/residuals/INTERVAL_NonImp_residuals_final.txt"


# ---------------------------------- #

# Function to recompute p-value from beta and se
safe_pnorm <- function(b, se, p=FALSE) {
  
  # Ensure the vectors are of the same length
  if(length(b) != length(se)) {
    stop("Beta and SE must be of the same length")
  }
  
  k  <- length(b)
  b  <- as.numeric(b)
  se <- as.numeric(se)
  
  # Initialize result vector with NA values
  result <- rep(NA, k)
  
  # Identify non-missing and non-zero indices
  i <- which(!is.na(b) & !is.na(se) & se != 0)
  
  # compute z-score and take absolute, raise digits with mpfr, apply pnorm for non-missing values
  z_score <- b[i] / se[i]
  z_mpfr <- Rmpfr::mpfr(- abs(z_score), 120)
  p_mpfr <- 2 * pnorm(z_mpfr)
  mlog10p <- - log10(p_mpfr)
  
  # print p-value in character format and mlog10p in numeric
  if(p==TRUE){
    # reformat to mpfr character, then to numeric (don't set digits for MLOG10P)
    mlog10p_mpfr <- Rmpfr::formatMpfr(p_mpfr, scientific = TRUE, digits = 6)
    result[i] <- mlog10p_mpfr
  } else {
    mlog10p_mpfr <- Rmpfr::formatMpfr(mlog10p, scientific = TRUE)
    result[i] <- as.numeric(mlog10p_mpfr)
  }
  
  return(result)
}

# ---------------------------------- #

# gene, seqid, locus, cis/trans, pQTL
#=========================#
# SDF2, seq.15333.11, chr17_25854347_28884532, cis, 17:26694861:A:G
# VTN,  seq.8280.238,	chr17_26592946_26739127, cis,	17:26694861:A:G
# Vitronectin, seq.13125.45, chr17_25283026_29975846, cis, 17:26694861:A:G

# proteins entailing VTN locus
seqid_list <- c("seq.15333.11", "seq.8280.238", "seq.13125.45")
snp_list <- c("17:26694861:A:G")

# Read psam and pvar
psam_df <- read.delim(path_psam, header = TRUE, comment.char = "")
pvar_df <- read.delim(path_pvar, header = TRUE, comment.char = "")
pvar_df <- dplyr::rename(pvar_df, CHROM = X.CHROM)

# read pgen
pvar <- pgenlibr::NewPvar(path_pvar)
pgen <- pgenlibr::NewPgen(path_pgen, pvar=pvar)


# Check the number of variants and samples.
n_variants <- pgenlibr::GetVariantCt(pgen)
n_samples  <- pgenlibr::GetRawSampleCt(pgen)

# Extract genotypes for all of variants
geno <- pgenlibr::ReadList(pgen, 1:n_variants, meanimpute = FALSE)

# Add variant IDs as column names
colnames(geno) <- pvar_df$ID

# Add sample IDs as row names
rownames(geno) <- psam_df$IID


geno_subset <- geno[, snp_list] %>% as.data.frame()
colnames(geno_subset) <- snp_list
rownames(geno_subset) <- psam_df$IID

geno_subset <- geno_subset %>% tibble::rownames_to_column(var = "IID")

# Read phenotype and merge it with genotype matrix
pheno <- fread(path_pheno)

# mediation input data
lm_df <- pheno %>%
  mutate(IID = as.character(IID)) %>%
  dplyr::select(IID, all_of(seqid_list)) %>%
  inner_join(geno_subset, join_by(IID))



# ---------------------------------- #
# ------  Mediation Analysis  ------ #
# ---------------------------------- #

# SDF2, seq.15333.11
# VTN,  seq.8280.238 (week pQTL)
# VTN,  seq.13125.45 (super signif. pQTL)

# model formula
my_fml <- "seq.13125.45 ~ `17:26694861:A:G`" # Step 1
my_fml <- "seq.13125.45 ~ `17:26694861:A:G` + seq.15333.11" # Step 2
my_fml <- "seq.15333.11 ~ `17:26694861:A:G` + seq.13125.45" # Step 4
my_int <- "seq.13125.45 ~ `17:26694861:A:G` * seq.15333.11" # interaction

my_fml <- "seq.8280.238 ~ `17:26694861:A:G`"
my_fml <- "seq.8280.238 ~ `17:26694861:A:G` + seq.15333.11"
my_int <- "seq.8280.238 ~ `17:26694861:A:G` * seq.15333.11"
my_fml <- "seq.15333.11 ~ `17:26694861:A:G` + seq.8280.238"

# run the model
model_fit <- lm(data = lm_df, formula = my_int)

# see model coefficients 
broom::tidy(model_fit) %>%
  mutate(Pval = safe_pnorm(estimate, std.error, p = TRUE)) %>%
  DT::datatable(caption = "model: VTN ~ SNP + SDF2 + SNP*SDF2") %>%
  DT::formatSignif(columns = c('estimate', 'std.error', 'statistic'), digits = 3)






library(ggplot2)

# Fit interaction model, then

# Prepare new data for prediction
newdat <- expand.grid(
  `17:26694861:A:G` = sort(unique(lm_df$`17:26694861:A:G`)),  # genotypes (0,1,2)
  seq.15333.11 = seq(min(lm_df$seq.15333.11, na.rm = TRUE),
                     max(lm_df$seq.15333.11, na.rm = TRUE),
                     length.out = 100)
)

# Predict fitted values and confidence intervals
newdat$pred <- predict(model_fit, newdata = newdat)
newdat$SNP_cat <- cut(newdat$`17:26694861:A:G`, breaks = c(0, 0.499,1.499, 2), labels = c("0", "1", "2"), include.lowest = TRUE)
newdat$genotype <- factor(newdat$SNP_cat, labels = c("0 alleles", "1 allele", "2 alleles"))

# Plot
ggplot(newdat, aes(x = seq.15333.11, y = pred, color = genotype)) +
  geom_line(size = 1.2) +





interactions::interact_plot(
  model = model_fit,
  pred = `17:26694861:A:G`,
  modx = seq.15333.11,
  modx.labels = c("Low", "Average", "High"),
  centered = "seq.8280.238", 
  linearity.check = F,
  line.thickness = 2,
  plot.points = F,
  #colors = "seagreen",
  x.label = "Dosage of 17:26694861:A:G",
  y.label = "VTN protein level (seq.8280.238)"
  ) +
  labs(
    #fill = "SDF2 protein level (seq.15333.11)",
    title = "Genotype modifies the SDF2–VTN relationship",
    subtitle = "Significant SNP × SDF2 interaction (p = 2×10⁻⁷)"
  ) +
  theme_classic(base_size = 14)

  # theme(
  #   axis.text = element_text(size = 12),
  #   axis.ticks.length = unit(2, "mm")
  #   )+

ggsave("08-Oct-25_SNP_SDF2_interaction_on_VTN_seq.8280.238.png", width = 7.5, height = 5.5, dpi = 300)
