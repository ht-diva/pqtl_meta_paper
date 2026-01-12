#=====================#
# Query C1QTNF4 region GWAS Catalog
#=====================#

# inputs: query results
path_query <- "/scratch/dariush.ghasemi/projects/pqtl_meta_paper/data/gwas-downloaded_2026-01-06_chr11_46589667_48595304_C1QTNF4.tsv"


query_res <- fread(path_query, sep = "\t", stringsAsFactors = FALSE)

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
  dplyr::mutate(
    gene = str_replace_all(MAPPED_GENE, "SLC39A13 - PSMC3", "SLC39A13"),
    gene = str_replace_all(gene, "OR4A.*$", "OR4A40P-OR4A46P"),
    gene = str_replace_all(gene, "EPIC1 x OR8D2 - OR8B7P", "OR8D2"),
    gene = str_replace_all(gene, "^SLC39A13-AS1, SPI1$|^SLC39A13-AS1, SPI1, SLC39A13$|^SPI1, SLC39A13-AS1$", "SLC39A13-AS1"),
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
    trait = str_remove_all(MAPPED_TRAIT, "measurement| serum | serum|serum | in blood| blood serum | amount"), #| protein level ratio
    trait = str_replace_all(trait, "glomerular filtration rate", "GFR"),
    trait = str_replace_all(trait, "IgG .*$", "IgG related"),
    trait = str_replace_all(trait, "high density lipoprotein .*$", "HDL-C"),
    trait = str_replace_all(trait, "low density lipoprotein .*$", "LDL-C"),
    trait = str_replace_all(trait, "triglyceride $|triglycerides:total lipids ratio", "triglycerides"),
    trait = str_replace_all(trait, "cholesterol:total lipids ratio, HDL-C|(total cholesterol|free cholesterol|cholesteryl).*, HDL-C|cholesteryl ester ", "HDL-C"),
    trait = str_replace_all(trait, "lipid , HDL-C", "HDL-C"),
    trait = str_replace_all(trait, "Alzheimer disease, polygenic risk score", "Alzheimer disease"),
    trait = str_replace_all(trait, "mood disorder, major depressive disorder", "major depressive disorder"),
    trait = str_replace_all(trait, "phospholipid level|phospholipids:total lipids ratio", "phospholipids"),
    trait = str_replace_all(trait, "gamma-glutamyl transferase .*$", "GGT"),
    trait = str_replace_all(trait, "alkaline phosphatase .*$", "ALP"),
    trait = str_replace_all(trait, "aspartate aminotransferase .*$", "AST"),
    trait = str_replace_all(trait, "alanine aminotransferase", "ALT"),
    trait = str_replace_all(trait, "body mass index", "BMI"),
    trait = str_replace_all(trait, "C-reactive protein.*$", "CRP"),
    trait = str_replace_all(trait, "(blood) urea nitrogen.*$", "BUN")
  ) %>%
  group_by(CHR_POS, trait, gene) %>% # take the most significant signal
  top_n(MLOG10P, n=1) %>%
  ungroup() %>%
  distinct(CHR_POS, trait, gene, .keep_all = TRUE) %>%   # avoid double counting
  filter(
    !str_detect(trait, "health trait"),
    MLOG10P > 9.3,
    CHR_POS > 47089667, CHR_POS < 48095304
    )

#----------------------------#
# -----   pheweb plot   -----
#----------------------------#
# pheweb plot
clean_res %>%
  group_by(trait) %>%
  mutate(n_hits_trait = n_distinct(CHR_POS)) %>% # no. of gwas hits per trait
  ungroup() %>%
  filter(n_hits_trait > 3) %>% #distinct(trait)
  ggplot(aes(x = fct_reorder(trait, -n_hits_trait),
             y = MLOG10P, color = gene, shape = gene)) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = seq(0, 300, 50)) +
  scale_shape_manual(values = rep(c(0:2, 7:10), 3)) + #c(19:17)
  guides(shape = guide_legend(ncol = 3)) + #, title = "Gene symbol"
  #scale_color_brewer(palette = "Paired")+
  #paletteer::scale_colour_paletteer_d("vapoRwave::hotlineBling")+
  #ggtitle("Query of ... region in GWAS Catalog") +
  theme_classic()+
  theme(
    legend.position = c(.75, .75),
    legend.background = element_blank(),
    legend.text = element_text(face=4),
    axis.text.x = element_text(size=8, face = 1, hjust = 0, vjust = 0.5, angle = -35),
    axis.title.x = element_blank(),
    axis.ticks.length = unit(2, "mm"),
    plot.margin = margin(r = 10, b = 2, t = 2, l = 2, unit = "mm")
  )


ggsave(
  filename = "12-Jan-26_query_of_celf1.png", 
  height = 6.5, width = 11, dpi = 300
)

#----------------------------#
# -----   Bar plot   -----
#----------------------------#

clean_res %>%
  count(trait, gene, sort = TRUE, name = "n_hits") %>%
  group_by(trait) %>%
  mutate(n_hits_trait = sum(n_hits)) %>%
  group_by(gene) %>%
  mutate(n_hits_gene = sum(n_hits)) %>%
  ungroup() %>%
  filter(n_hits_trait > 3) %>% 
  #distinct(trait, n_hits_trait) %>% arrange(-n_hits_trait)
  distinct(gene, n_hits_gene) %>% arrange(-n_hits_gene)
  ggplot(aes(x = n_hits,
             y = fct_reorder(trait, n_hits_trait),
             fill = fct_reorder(gene, n_hits_gene))) +
  geom_bar(position="stack", stat="identity") +
  scale_x_continuous(breaks = c(seq(0, 55, 5)), limits = c(0,56), 
                     expand = c(0,0)) + # remove space between plot and y-axis
  guides(fill = guide_legend(ncol = 2)) +
  labs(fill = "Gene symbol") +
  theme_light() +
  theme(
    legend.position = c(.65, .20),
    legend.background = element_blank(),
    legend.text = element_text(face=4),
    axis.text.y = element_text(size = 10, face = 1),
    axis.title.y = element_blank(),
    axis.ticks.length = unit(2, "mm"),
    plot.margin = margin(r = 2, b = 2, t = 2, l = 2, unit = "mm")
  )

ggsave(
  filename = "12-Jan-26_query_of_celf1_gene_distribution.png", 
  height = 11, width = 7.5, dpi = 300
)
