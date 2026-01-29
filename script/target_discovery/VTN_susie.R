

library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)
library(glue)
library(susieR)

#-------------------#
# read susie output
path_susie_VTN <- "/scratch/dariush.ghasemi/projects/pqtl_susie/results"

vtn_fit <- readRDS(glue(path_susie_VTN, "/VTN_interval/susierss/cs_fitness/seq.13125.45_17_25283026_29975846_fit.rds"))

susie_vtn_interval <- fread(
  glue(path_susie_VTN, "/VTN_interval/susierss/collected_credible_sets.tsv")
  )

susie_vtn_meta <- fread(
  glue(path_susie_VTN, "/VTN/susierss/collected_credible_sets.tsv")
)

#-------------------#
# verify results
vtn_fit
summary(vtn_fit)
susie_get_cs(vtn_fit, X = NULL)

susie_vtn_interval %>% summarise( n(), .by = c(seqid, locus))

#-------------------#
# number of SNPs in each credible set
tbl_interval <- susie_vtn_interval %>%
  select(seqid, locus, cs_snps) %>%
  na.omit() %>%
  mutate(ncs = str_count(cs_snps, ",") + 1)

# number of SNPs in each credible set
tbl_meta <- susie_vtn_meta %>%
  select(seqid, locus, cs_snps) %>%
  na.omit() %>%
  mutate(ncs = str_count(cs_snps, ",") + 1)

#-------------------#
# number of credible set in each dataset
tbl_interval %>%
  mutate(study = "interval") %>%
  rbind(tbl_meta %>% mutate(study = "meta")) %>%
  summarise(n_sets = n(), .by = c(seqid, locus, study)) %>%
  ggplot(aes(x = seqid, y = n_sets, fill = study)) +
  geom_col(position = position_dodge(.9)) +
  geom_text(aes(label = n_sets), position = position_dodge(.9)) + 
  scale_y_continuous(breaks = seq(0, 10, 1)) +
  scale_fill_manual(name = "dataset", values = c("steelblue2", "orange"))+
  labs(x = "", y = "Number of credible sets") +
  theme_light()+
  theme(panel.grid.minor = element_blank())

ggsave("14-Oct-25_susie_VTN_interval_vs_meta.png", width = 7.5, height = 5.5, dpi = 300)

#-------------------#
# save number of SNPs in each credible set
tbl_interval %>% select(- cs_snps) %>% write.csv()
tbl_meta %>% select(- cs_snps) %>% write.csv()




