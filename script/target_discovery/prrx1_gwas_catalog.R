#=====================#
# Query PRRX1 region 
#   in GWAS Catalog
#    21 Jan, 2026
#=====================#

# inputs: query results
path_query <- "/scratch/dariush.ghasemi/projects/pqtl_meta_paper/data/gwas_catalog_2026-01-21-locus_1_170162728_171239421_PRRX1.tsv"

# read results
query_res <- fread(path_query, sep = "\t", stringsAsFactors = FALSE)


#----------------------------------------#
#-----          PheWeb Plot         -----
#----------------------------------------#

query_res %>%
  filter(
    str_detect(MAPPED_TRAIT, "protein"),
    PVALUE_MLOG > 7.2
    ) %>%
  group_by(CHR_POS) %>%
  slice_max(PVALUE_MLOG, n = 1) %>% # pick most significant signal across identical SNPs
  ggplot(aes(y = `DISEASE/TRAIT`, x = PVALUE_MLOG, color = MAPPED_GENE)) +
  geom_point(
    position = position_dodge(.7)
    ) +
  labs(x = "MLOG10P")+
  theme_light()+
  theme(
    legend.position = c(0.8, 0.15),
    legend.background = element_blank(),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10, face = 3)
  )


ggsave(
  "27-Jan-26_query_gwas_catalog_PRRX1_assoc_proteins.jpg",
  plot = last_plot(), width = 11, height = 6.5, dpi = 150
  )

