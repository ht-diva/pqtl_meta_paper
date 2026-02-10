

path_eqtl_prrx1 <- "/project/cdh/gtex/data/Heart_Left_Ventricle/Heart_Left_Ventricle.v10.ENSG00000116132.12.csv.gz"
path_eqtl_interval <- "/exchange/healthds/pQTL/INTERVAL/summary_stats/eQTL/cis/INTERVAL_eQTL_nominal_chr1.tsv"
path_eqtl_trans <- "/exchange/healthds/pQTL/INTERVAL/summary_stats/eQTL/trans/INTERVAL_trans_eQTL_summary_statistics_1e5.tsv"
path_eqtl_b37 <- "/scratch/dariush.ghasemi/projects/pqtl_liftover/results/vcf/output_positions_GRCh37.txt"
path_geno_gtex <- "/exchange/healthds/private/GTEXv8_FUSION/GTEx.WGS.838.passOnly.geno0.05.hwe0.00001.dbsnp.SNPsOnly.NoAmbig_1.bim"



# pQTL meta-analysis for PRRX1 belonged to seq.24215.8
# pQTL region in meta-analysis (hg37): chr1:169636068-171152683
# ENSG00000116132 -- PRRX1
# eQTL locus in GTEx:
#   - GRCh37: chr1:169632033-171631688
#   - GRCh38: chr1:169662795-171662548


#----------------------------------------#
#-----     Store GTEx eQTLs RSID    -----
#----------------------------------------#

# subset region (GTEx genotype in build 38 -- bed/bim/fam)
gtex_geno <- fread(
  path_geno_gtex,
  col.names = c("chr", "rsid", "distance", "pos", "ALT", "REF"),
  colClasses = c(rep("character", 3), "integer", rep("character", 2))
)

prrx1_snps <- gtex_geno %>%
  filter(
    between(pos, 169662795, 171662548)
    #chr == 1, pos >= 169662795 & pos <= 171662548
  )

# Store list of SNP ids without multi-allelics
prrx1_snps %>%
  select(rsid) %>%
  write.table(
    "/scratch/dariush.ghasemi/projects/pqtl_susie/plink_ld/prrx1_eqtls_gtex_v8.list",
    quote = F, row.names = F, col.names = F)


#----------------------------------------#
#-----      Wrangle Input Data      -----
#----------------------------------------#

# https://gtexportal.org/home/faq
# The variant IDs from GTEx are in the format:
#       chromosome_position_ref_alt_build.


# Read data
prrx1_eqtl <- fread(path_eqtl_prrx1)
# eqtl_interval <- fread(path_eqtl_interval, fill = TRUE)
# eqtl_trans <- fread(path_eqtl_trans)

# Read GTEx lifted position to hg37
eqtl_b37 <- fread(path_eqtl_b37, header = F,
                  col.names = c("chr", "pos_b37", "variant_id"))

# Drop junk columns
eqtl_interval[, c("V14", "V16") := NULL]


# check if PRRX1 is among eQTL gene in Interval study --- No, it was NOT
eqtl_trans %>% 
  #distinct(phenotype_id) %>%
  filter(phenotype_id == "ENSG00000116132")


# Append position in b37 to GTEx eQTLs for PRRX1 after lift-over
prrx1_eqtl_hg37 <- prrx1_eqtl %>% left_join(eqtl_b37, join_by(variant_id))

message(
  "eQTL region in GTEx was (hg37): ",
  min(prrx1_eqtl_hg37$pos_b37),
  ",",
  max(prrx1_eqtl_hg37$pos_b37)
  )


#----------------------------------------#
#-----      eQTL Regional Plot      -----
#----------------------------------------#

# eQTLs plot
prrx1_eqtl_hg37 %>%
  mutate(MLOG10P = -log10(pval_nominal)) %>%
  ggplot(aes(x = pos_b37, y = MLOG10P, color = slope))+
  geom_point() +
  scale_color_gradient2(
    #low = "green3", mid = "grey50", high = "brown", midpoint = .02
    low = "#075AFF", mid = "#FFFFCC", high = "#FF0000"
    ) +
  scale_y_continuous(breaks = seq(0,10,1), limits = c(0,10), expand = c(0, 0)) +
  scale_x_continuous(
    breaks = c(seq(169612033, 171631688, 2e5)),
    labels = c(seq(169.6, 171.6, .2)),
    limits = c(169600000, 171632000),
    expand = c(0, 0)
    ) +
  labs(
    color = "Effect size",
    x = "\nGenomic position (Mb, GRCh37)",
    title = "PRRX1 eQTLs at Heart Left Ventricle in GTEx v10"
    ) +
  theme_light() +
  theme(
    legend.position = c(0.91, 0.75),
    legend.background = element_blank(),
    axis.title = element_text(size = 12, face = 2),
    axis.text = element_text(size = 11)
  )


ggsave(
  filename = "27-Jan-26_query_of_prrx1_eqtls_gtex.png", 
  height = 6.5, width = 11, dpi = 300
)
