
# PRRX1 pQTLs in Interval -- lifted to GRCh38
path_pqtl_b38 <- "/scratch/dariush.ghasemi/projects/pqtl_liftover/results/vcf/input_positions_interval_GRCh37_GRCH38_converted.txt"
path_pqtl_ld_matrix <- "/scratch/dariush.ghasemi/projects/pqtl_susie/plink_ld/ld/seq.24215.8_1_169258178_171220925_ld.matrix"
path_pqtl_ld_header <- gsub(".matrix", ".headers", path_pqtl_ld_matrix)


# Save list of pQTLs in INTERVAL to subset INTERVAL genotype 
pqtl_interval %>% select(SNPID) %>%
  write.table(
    "/scratch/dariush.ghasemi/projects/pqtl_susie/plink_ld/seq.24215.8_1_169258178_171220925_snps.list",
    quote = F, row.names = F, col.names = F)

# Read lifted pQTLs
pqtl_b38 <- fread(
  path_pqtl_b38,
  header = F,
  col.names = c("chr", "pos_b38", "variant_id")
  )

# Read LD Correlation Matrix
ld_headers <- fread(path_pqtl_ld_header, header = FALSE, col.names = "SNP")
R <- fread(path_pqtl_ld_matrix) %>% as.matrix()
rownames(R) <- colnames(R) <- ld_headers$SNP

# common pQTLs between two builds 
common_pqtls <- pqtl_b38 %>%
  mutate(
    EA  = str_split_fixed(variant_id, ":", 4)[,3], # in Interval, SNP ids are in chr_pos_ea_nef allele.
    NEA = str_split_fixed(variant_id, ":", 4)[,4],
    SNPID_b38 = str_c(chr, pos_b38, EA, NEA, sep = ":") %>% str_remove("chr")
    ) %>%
  inner_join(ld_headers, join_by(variant_id == SNP))

# subset LD to common SNPs
R <- R[common_pqtls$variant_id, common_pqtls$variant_id]
rownames(R) <- colnames(R) <- common_pqtls$SNPID_b38

# Append positions in build GRCh38 to INTERVAL pQTLs
sumstat <- pqtl_interval %>%
  inner_join(common_pqtls %>% select(variant_id, SNPID_b38, pos_b38), join_by(SNPID == variant_id)) %>% 
  mutate(SNPID = SNPID_b38, POS = pos_b38) %>%
  select(- SNPID_b38,  - pos_b38)

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

# Annotate and save full model fitness with LD matrix for coloc
res_rss_annot_pqtls <- coloc::annotate_susie(res_rss, sumstat$SNPID, R)

saveRDS(res_rss_annot_pqtls,
        file = "susie_annotated_res_for_prrx1_interval_pqtls.RDS")


#-------------------------------#
# -----  Coloc with SuSiE  -----
#-------------------------------#

# run coloc with susie
library(coloc)

res_coloc <- coloc.susie(
  dataset1 = res_rss_annot_eqtls,
  dataset2 = res_rss_annot_pqtls,
  back_calculate_lbf = FALSE,
  susie.args = list()    # only needed if you want coloc to re-run susie via runsusie
)

res_coloc$summary

