

# Supplementary figure 7 of pQTL meta-analysis paper
# Number of conditional variants in cis and trans loci

path_freez <- "/exchange/healthds/pQTL/results/META_CHRIS_INTERVAL/Locus_breaker_cojo_frozen_version_1812024/"
path_st3_epitop <- paste0(path_freez, "06-Oct-25_suppl_table_3_cojo_variants_vep_annotated_credible_sets_epitope.tsv")

#-------------------#
# read COJO results (supplementary table 3 in paper)
st3_epitop <- fread(path_st3_epitop)


# check data
st3_epitop %>%
  distinct(study_id, locus, SNP, cis_or_trans, ncojo) %>%
  arrange(-ncojo)

#-------------------#
# bar plot
st3_epitop %>%
  distinct(study_id, locus, SNP, cis_or_trans, ncojo) %>%
  # group_by(seqid, locus) %>%
  # summarise(n_cojo = n()) %>% arrange(n_cojo) %>% View
  ggplot(aes(ncojo, fill = cis_or_trans)) +
  geom_histogram(bins = 18,
                 stat = "count", 
                 position = position_dodge(0.9, preserve = 'single')) + # keep half column empty
  scale_x_continuous(breaks = c(1:18)) +
  scale_y_continuous(breaks = seq(0, 5500, 500)) +
  #scale_fill_manual(values = c("cis" = "orange2", "trans" = "steelblue2"))+
  labs(
    fill = "Locus",
    x = "\nNumber of conditionally independent SNPs per locus",
    y = "Number of loci\n"
  ) +
  theme_light() +
  theme(
    legend.position = c(0.8, 0.8),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.title = element_text(size = 13, face = 2),
    legend.text = element_text(size = 11, face = 4),
    legend.key.height = unit(2, "mm"),
    legend.key.width = unit(8, "mm"),
    legend.key.spacing = unit(4, "mm"),
    legend.key.spacing.y = unit(2, "mm"),
    axis.title = element_text(size = 12, face = 2),
    axis.text  = element_text(size = 11),
    axis.ticks.length = unit(2, "mm")
  )



ggsave("06-Nov-25_ncojo_per_locus_cis_trans_default.jpg",
       plot = last_plot(), width = 8, height = 6, dpi = 300)

#-------------------#



