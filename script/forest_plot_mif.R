
library(readxl)

path_st7_mr <- "/scratch/dariush.ghasemi/projects/pqtl_meta_paper/data/supplementary_table_7_Nov18.xlsx"

st7mr <- read_excel(path_st7_mr, sheet = 2, range = NULL)

st7mr %>%
  separate(locus_START_END_37, into = c("chrom", "beg", "end"), sep = "_") %>%
  mutate(target_endpoint = paste0(HARMONIZED_GENE_NAME, " - ", PHENOTYPE_description)) %>%
  dplyr::filter(
    CHR == 22,
    #locus_START_END_37 == "chr22_23942068_26128449"
    beg >= 23942068,
    end <= 26128449
  ) %>%
  ggplot(aes(
    y = fct_reorder(target_endpoint, -BETA),
    x = BETA,
    xmin = BETA - 1.96*SE,
    xmax = BETA + 1.96*SE,
    colour = SeqID
    )) +
  geom_pointrange(position = position_dodge(0.9)) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(breaks = seq(-0.1, 0.6, .1)) +
  labs(x = "Causal effect (95% confidence interval)", y = "Target gene-endpoint") +
  theme_light() +
  theme(
    legend.position = c(0.8, 0.8),
    legend.background = element_blank(),
    axis.title = element_text(size = 12, face = 2),
    axis.text  = element_text(size = 11),
    axis.ticks.length = unit(2, "mm")
    )


ggsave("18-Nov-25_forest_plot_mif.jpg",
       plot = last_plot(), width = 8, height = 6, dpi = 300)
