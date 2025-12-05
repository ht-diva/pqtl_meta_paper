

library(tidyverse)
library(readxl)
library(data.table)
library(gtsummary)

path_base <- "/scratch/dariush.ghasemi/projects/pqtl_meta_paper/data/"
path_st2_map  <- glue::glue(path_base, "st2_mapping_02Dec.xlsx")
path_st3_loci <- glue::glue(path_base, "st3_loci_list_Nov26.xlsx")
path_st4_cojo <- glue::glue(path_base, "supplementary_table_4_Nov17_4plotting.xlsx")
path_freez <- "/exchange/healthds/pQTL/results/META_CHRIS_INTERVAL/Locus_breaker_cojo_frozen_version_1812024/"
path_lb <- "16-Dec-24_collected_loci_excluding_mhc.csv"

st2_map <- read_excel(path_st2_map)
st3_lb  <- read_excel(path_st3_loci)
st4_cojo <- read_excel(path_st4_cojo, sheet = NULL, range = NULL)


res <- st4_cojo %>%
  dplyr::rename(
    seqid = SeqID,
    locus = locus_START_END_37,
    location = cis_or_trans,
    status = uniprot_new_in_somascan7k_vs5k
  ) %>%
  group_by(seqid) %>%
  dplyr::mutate(nloci = n_distinct(locus)) %>% # Count number of loci per seqid
  group_by(seqid, locus) %>%
  dplyr::mutate(
    ncojo = n_distinct(SNPID), # No. of COJO SNPs per seqid-locus
    nover = ncojo > 1,         # loci with secondary signals
    signal = ifelse(ncojo > 1, "secondary", "primary") %>% factor(levels = c("primary", "secondary")),
    status = ifelse(status, "new", "old") %>% factor(levels = c("old", "new"))
    ) %>%
  select(seqid, locus, location, status, nloci, ncojo, signal) %>%
  ungroup() %>% 
  distinct()


# appending Z-scores of index variants
lb_cojo <- st3_lb %>%
  dplyr::mutate(
    #locus = paste0("chr", chr, "_", start, "_", end),
    locus = locus_START_END_37,
    zscore = abs(BETA/SE)
  ) %>%
  dplyr::select(seqid = SeqID, locus, SNPID, zscore) %>%
  inner_join(res, join_by(seqid, locus)) %>%
  filter(location == "cis")

# appending LOD and CV from mapping to LB summarized
map_cojo <- st2_map %>%
  dplyr::select(
    seqid = SeqID,
    INTERVAL_LOD_batch1,
    INTERVAL_LOD_batch2,
    INTERVAL_CV_batch1,
    INTERVAL_CV_batch2,
    CHRIS_LOD,
    CHRIS_CV,
    ) %>%
  distinct() %>%
  right_join(res, join_by(seqid)) %>%
  filter(location == "cis")


lb_pathway <- st2_map %>%
  dplyr::select(
    seqid = SeqID,
    pathway = `Secretion pathway`
  ) %>%
  summarise(func_prot = paste0(unique(pathway), collapse = ","), .by = seqid) %>%
  right_join(res, join_by(seqid)) %>%
  dplyr::filter(func_prot %in% c("intracellular", "membrane", "secreted")) %>%
  #count(func_prot) %>% DT::datatable()
  dplyr::filter(location == "cis")




#-------------------#
# contingency table 
res %>%
  dplyr::filter(location == "cis") %>%
  summarise(
    n_pri = sum(signal == "primary"),
    prop_pri = n_pri/n(),
    prop_pri_tot = n_pri/nrow(.),
    n_sec = sum(signal == "secondary"),
    prop_sec = n_sec/n(),
    prop_sec_tot = n_sec/nrow(.),
    .by = status
    )


table(res$status, res$location)  %>% prop.table(margin = 1)
table(res$location, res$nover) 
table(res$status, res$nover)  
prop.table(table(res$status, res$nover), margin = 1)


lb_pathway %>%
  dplyr::filter(func_prot == "membrane") %>%
  summarise(
    n_pri = sum(signal == "primary"),
    prop_pri = n_pri/n(),
    #prop_pri_tot = n_pri/nrow(.),
    n_sec = sum(signal == "secondary"),
    prop_sec = n_sec/n(),
    #prop_sec_tot = n_sec/nrow(.),
    .by = status
  )


#-------------------------------#
# -----   Logistic Model   -----
#-------------------------------#


# Logistic regression
fit_logistic <- glm(
  #signal ~ location,
  #signal ~ status,
  #signal ~ status * location,
  #signal ~ status + zscore,
  #signal ~ status + CHRIS_CV,
  signal ~ status + func_prot,
  family = binomial(link='logit'),
  data = lb_pathway
  )

gtsummary::tbl_regression(fit_logistic, exponentiate = TRUE)

broom::tidy(
  fit_logistic,
  conf.int = TRUE,
  exponentiate = TRUE
  ) #%>% DT::datatable()




#-------------------#
# forest plot
coef_df <- broom::tidy(fit_logistic, conf.int = TRUE, exponentiate = TRUE)

ggplot(coef_df[-1,], aes(
  y = term,
  x = estimate,
  xmin = conf.low,
  xmax = conf.high
)) +
  geom_pointrange() +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  labs(
    x = "Odds ratio (log scale)",
    y = "",
    title = "Effect sizes from logistic regression"
  ) +
  theme_bw(15)

#-------------------#
# bar plot
res %>%
  ggplot(aes(signal, fill = status)) +
  geom_histogram(
    stat = "count", #bins = 18,
    position = position_dodge(0.9, preserve = 'single')
    ) +
  facet_wrap(~location) +
  stat_count(
    aes(label=..count.., y=..count.. + 130), 
    geom = 'text', color = '#bd5734',
    position = position_dodge(.9),
  ) +
  scale_y_continuous(breaks = seq(0, 4600, 500)) +
  #scale_x_continuous(breaks = c(1:18)) +
  scale_fill_manual(
    name = "Protein",
    labels = c("Present in 5k", "New in 7k"),
    values = c("#50394c", "#d4ac6e") #"#d8b365", "#5ab4ac"
    )+ 
  labs(
    #fill = "",
    #x = "\nNumber of conditionally independent SNPs per locus",
    x = "\nNumber of COJO signals per locus",
    y = "Number of loci\n"
  ) +
  theme_light() +
  theme(
    legend.position = c(0.15, 0.85),
    legend.background = element_blank(),
    legend.title = element_text(size = 13, face = 2),
    legend.text = element_text(size = 11, face = 2),
    legend.key.height = unit(2, "mm"),
    legend.key.width = unit(8, "mm"),
    legend.key.spacing = unit(4, "mm"),
    legend.key.spacing.y = unit(2, "mm"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12, color = "Black", face = 4),
    panel.border = element_rect(colour = "black", fill = NA),
    #panel.grid.minor = element_blank(),
    #panel.grid.major = element_line(color = "lightgray", size = .25),
    strip.text = element_text(size = 10, face = 4),
    axis.title = element_text(size = 12, face = 2),
    axis.text  = element_text(size = 11),
    axis.ticks.length = unit(2, "mm")
  )



ggsave("20-Nov-25_ncojo_per_locus_new_old.jpg",
       plot = last_plot(), width = 8, height = 6, dpi = 300)

#-------------------#
# box plot
lb_cojo %>%
  ggplot(aes(x = status, y = zscore, fill = signal))+
  geom_boxplot() +
  scale_x_discrete(
    name = "Protein",
    labels = c("Present in 5k", "New in 7k")
    )+
  scale_fill_manual(
    name = "Signal",
    values = c("#50394c", "#d4ac6e")
  )+
  theme_light() +
  theme(
    legend.position = c(0.15, 0.85),
    legend.background = element_blank(),
    legend.title = element_text(size = 13, face = 2),
    legend.text = element_text(size = 11, face = 2),
    legend.key.height = unit(8, "mm"),
    legend.key.width = unit(8, "mm"),
    legend.key.spacing = unit(3, "mm"),
    axis.title = element_text(size = 12, face = 2),
    axis.text  = element_text(size = 11),
    axis.ticks.length = unit(2, "mm")
  )


ggsave("25-Nov-25_ncojo_per_locus_new_old_adjusted.jpg",
       plot = last_plot(), width = 8, height = 6, dpi = 300)

#-------------------#
# box plot for LOD
map_cojo %>%
  ggplot(aes(x = status, y = INTERVAL_CV_batch1, fill = signal))+
  geom_boxplot(color = "grey40")+
  scale_x_discrete(
    name = "Protein",
    labels = c("Present in 5k", "New in 7k")
  )+
  scale_fill_manual(
    name = "Signal",
    values = c("#50394c", "#d4ac6e")
  )+
  theme(legend.position = c(0.85, 0.85))



ggstatsplot::ggbetweenstats(data = map_cojo, x = status, y = INTERVAL_CV_batch2)


# bar plot for protein function/location in blood
lb_pathway %>%
  ggplot(aes(signal, fill = status)) +
  geom_histogram(
    stat = "count",
    position = position_dodge(0.9, preserve = 'single')
  ) +
  facet_wrap(~func_prot) +
  stat_count(
    aes(label=..count.., y=..count.. + 20), 
    geom = 'text', color = '#bd5734',
    position = position_dodge(.9),
  ) +
  scale_y_continuous(breaks = seq(0, 650, 50)) +
  scale_fill_manual(
    name = "Protein",
    labels = c("Present in 5k", "New in 7k"),
    values = c("#50394c", "#d4ac6e")
  )+ 
  labs(
    x = "\nNumber of COJO signals per locus",
    y = "Number of loci\n"
  ) +
  theme_light() +
  theme(
    legend.position = c(0.15, 0.85),
    legend.background = element_blank(),
    legend.title = element_text(size = 13, face = 2),
    legend.text = element_text(size = 11, face = 2),
    legend.key.height = unit(2, "mm"),
    legend.key.width = unit(8, "mm"),
    legend.key.spacing = unit(4, "mm"),
    legend.key.spacing.y = unit(2, "mm"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12, color = "Black", face = 4),
    panel.border = element_rect(colour = "black", fill = NA),
    strip.text = element_text(size = 10, face = 4),
    axis.title = element_text(size = 12, face = 2),
    axis.text  = element_text(size = 11),
    axis.ticks.length = unit(2, "mm")
  )

ggsave("02-Dec-25_barplot_new_old_proteins_pathway_cis.jpg",
       plot = last_plot(), width = 9, height = 6, dpi = 300)

