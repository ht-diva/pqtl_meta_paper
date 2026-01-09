#=====================#
# Query C1QTNF4 region GWAS Catalog
#=====================#

# inputs: query results
path_query <- "/scratch/dariush.ghasemi/projects/pqtl_meta_paper/data/gwas-downloaded_2026-01-06_chr11_46589667_48595304_C1QTNF4.tsv"


query_res <- fread(path_query, sep = "\t")

# clean query results
clean_res <- query_res %>%
  dplyr::select(
    CHR_ID,
    CHR_POS,
    SNPS,
    reported_gene = "REPORTED GENE(S)",
    MAPPED_GENE,
    disease_trait = "DISEASE/TRAIT",
    MAPPED_TRAIT,
    MLOG10P = PVALUE_MLOG,
    effect_size = "OR or BETA"
  ) %>% 
  dplyr::filter(
    MLOG10P > 9.3,
    CHR_POS > 47089667, CHR_POS < 48095304
    ) %>%
  dplyr::mutate(
    gene = str_replace_all(MAPPED_GENE, "SLC39A13 - PSMC3", "SLC39A13"),
    gene = str_replace_all(gene, "OR4A.*$", "OR4A40P-OR4A46P"),
    gene = str_replace_all(gene, "EPIC1 x OR8D2 - OR8B7P", "OR8D2"),
    gene = str_replace_all(gene, "SLC39A13-AS1, SPI1, SLC39A13|SPI1, SLC39A13-AS1", "SLC39A13-AS1, SPI1"),
    gene = str_replace_all(gene, "AGBL2 - FNBP4", "AGBL2"),
    gene = str_replace_all(gene, "FNBP4.*$", "FNBP4"),
    gene = str_replace_all(gene, "NR1H3, ACP2|NR1H3-DDB2", "NR1H3"),
    gene = str_replace_all(gene, "ACP2, DDB2|DDB2, ACP2|ACP2|DDB2", "ACP2-DDB2"),
    gene = str_replace_all(gene, "NDUFS3.*$|KBTBD4, NDUFS3", "NDUFS3"),
    gene = str_replace_all(gene, "MADD - MYBPC3|MADD-AS1, MADD", "MADD"),
    gene = str_replace_all(gene, "MYBPC3 - SPI1.*$", "MYBPC3"),
    gene = str_replace_all(gene, "PACSIN3.*$", "PACSIN3"), 
    gene = str_replace_all(gene, "PSMC3.*$", "PSMC3"),
    gene = str_replace_all(gene, "FAM180B - C1QTNF4$|^C1QTNF4$|^FAM180B$", "FAM180B-C1QTNF4"),
    gene = str_replace_all(gene, "NUP160 - PTPRJ|NUP160|PTPRJ", "NUP160-PTPRJ"),
    gene = str_replace_all(gene, "RAPSN - CELF1|RAPSN|CELF1", "RAPSN-CELF1"),
    #gene = str_replace_all(gene, "MTCH2 - AGBL2|MTCH2|AGBL2", "AGBL2-MTCH2"),
    # gene = str_replace_all(gene, ".* DERL3|DERL3.*$", "SMARCB1"),
    # gene = str_replace_all(gene, ".* CABIN1|CABIN1.*$", "GSTT4"),
    trait = str_remove_all(MAPPED_TRAIT, "measurement| serum | serum|serum | in blood| blood serum | amount"), #| protein level ratio
    trait = str_replace_all(trait, "glomerular filtration rate", "GFR"),
    trait = str_replace_all(trait, "IgG .*$", "IgG related"),
    trait = str_replace_all(trait, "high density lipoprotein .*$", "HDL-C"),
    trait = str_replace_all(trait, "low density lipoprotein .*$", "LDL-C"),
    trait = str_replace_all(trait, "gamma-glutamyl transferase .*$", "GGT"),
    trait = str_replace_all(trait, "alkaline phosphatase .*$", "ALP"),
    trait = str_replace_all(trait, "aspartate aminotransferase .*$", "AST"),
    trait = str_replace_all(trait, "alanine aminotransferase .*$", "ALT"),
    trait = str_replace_all(trait, "body mass index", "BMI"),
    trait = str_replace_all(trait, "C-reactive protein", "CRP"),
    trait = str_replace_all(trait, "urea nitrogen", "BUN"),
    trait = fct_reorder(trait, -MLOG10P), 
  )


# pheweb plot
clean_res %>%
  group_by(trait) %>%
  mutate(
    n_signals = n_distinct(CHR_POS)
    ) %>%
  ungroup() %>%
  filter(
    n_signals > 3,
    #!str_detect(trait, "in medium HDL")
    ) %>%
  ggplot(aes(x = fct_reorder(trait, - n_signals),
             y = MLOG10P, color = gene, shape = gene)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = seq(0, 300, 50)) +
  scale_shape_manual(values = rep(c(0:2, 7:10), 3)) + #c(19:17)
  guides(shape = guide_legend(ncol = 3)) +
  #scale_color_brewer(palette = "Paired")+
  #scale_y_discrete(position = "right") +
  #paletteer::scale_colour_paletteer_d("vapoRwave::hotlineBling")+
  #ggtitle("Query of ... region in GWAS Catalog") +
  theme_classic()+
  theme(
    legend.position = c(.75, .7),
    legend.background = element_blank(),
    axis.text.x = element_text(size=8, face = 1, hjust = 0, vjust = 0.5, angle = -90),
    axis.title.x = element_blank(),
    axis.ticks.length = unit(2, "mm"),
    plot.margin = margin(r = 2, b = 2, t = 2, l = 2, unit = "mm")
  )


ggsave(
  filename = "09-Jan-26_query_of_celf1.png", 
  height = 6.5, width = 11, dpi = 300
)


