
library(susieR)


#============================#
#  SuSiE Parameters
#============================#

# Load parameters for susieR model
susie_min_abs_cor <- 0.5
susie_iter <- 5000
susie_L <- 10
susie_est_resvar <- FALSE

#============================#
#  SuSiE Inputs
#============================#

# Set TRUE to compute correlation from X, FALSE to load pre-computed LD
path_ld_plink1 <- "/scratch/dariush.ghasemi/projects/pqtl_susie/plink_ld/ld/prrx1_eqtls_gtex_v8_ld.ld"
#path_ld_matrix <- "/scratch/dariush.ghasemi/projects/pqtl_susie/plink_ld/ld/prrx1_eqtls_gtex_v8_ld.matrix"
path_ld_header <- gsub(".matrix", ".headers", path_ld_matrix)


# How many variants (eQTLs) overlap in GTEx and INTERVAL
intersect(sumstat$SNPID, ld_headers$SNP) %>% length()
setdiff(sumstat$SNPID, ld_headers$SNP) %>% length()

# Read LD Correlation Matrix
ld_headers <- fread(path_ld_header, header = FALSE, col.names = "SNP")

# Add positions in hg38 from GTEx genotype
ld_headers <- ld_headers %>%
  left_join(prrx1_snps, join_by(SNP == rsid)) %>%
  mutate(SNPID = str_c(chr, pos, ALT, REF, sep = ":"))
  

R <- fread(path_ld_plink1) %>% as.matrix()
rownames(R) <- colnames(R) <- ld_headers$SNPID

# Once computing LD, only 6569 of 8327 GTEx eQTLs existed in INTERVAL 
#sumstat <- prrx1_eqtl_harm %>% filter(SNPID %in% ld_headers$SNP)

# Computing LD in GTEx resulted
R <- R[sumstat$SNPID, sumstat$SNPID]

dim(R)
dim(sumstat)

#============================#
#  SuSiE RSS
#============================#

betas    <- sumstat$BETA
se_betas <- sumstat$SE
n        <- median(sumstat$N, na.rm = TRUE)

# The estimated λ is
lambda <- estimate_s_rss(betas/se_betas, R=R, n=n)

err_handling <- function(e) { stop("❌ SuSiE failed: ", e$message) }

res_rss <- tryCatch(
  susie_rss(
    bhat = betas,
    shat = se_betas,
    n = n,
    R = R,
    L = susie_L,
    max_iter = susie_iter,
    min_abs_corr = susie_min_abs_cor,
    estimate_residual_variance = susie_est_resvar # TRUE if using in-sample LD
  ),
  error = err_handling
)



# plot credible sets
susie_plot(
  res_rss,
  y = "PIP",
  b = betas,
  xlab = "Variants",
  add_bar = FALSE,
  add_legend = TRUE,
  #main = paste("SeqID:", tag_seqid, "\nRegion:", tag_locus)
)


# full model summary
full_res <- summary(res_rss)
cs   <- full_res$cs    # containing CS impurity indices
vars <- full_res$vars  # containing CS Posterior Inclusion Probabilities

# list of the entire SNPs with PIP
snps_pip <- vars %>%
  transmute(
    cs_id = cs,
    SNPID = sumstat$SNPID[variable],
    PIP = variable_prob
  )

sumstat %>%
  left_join(snps_pip, by = "SNPID") %>%
  filter(cs_id > 0)


# Annotate and save full model fitness with LD matrix for coloc
res_rss_annot_eqtls <- coloc::annotate_susie(res_rss, sumstat$SNPID, R)

saveRDS(res_rss_annot_eqtls,
        file = "susie_annotated_res_for_prrx1_gtex_eqtls.RDS")
