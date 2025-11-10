setwd("/group/diangelantonio/users/alessia_mapelli/pQTL/INTERVAL/Hotspots and lonespots/New_coloc_results")
rm(list=ls(all=TRUE))

library(readr)
library(dplyr)
library(scales)
library(ggrepel)
library(data.table)
library(tidyverse)

mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann <- read_delim("/exchange/healthds/pQTL/results/META_CHRIS_INTERVAL/Locus_breaker_cojo_frozen_version_1812024/mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann.csv", 
                                                              delim = ";", escape_double = FALSE, trim_ws = TRUE)

colnames(mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann)
bar_plot_hotspot_signal_df <- data.frame(n_signals = nrow(mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann))
bar_plot_hotspot_signal_df$n_signals_cis <- sum(mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann$cis_or_trans == "cis")
bar_plot_hotspot_signal_df$n_signals_trans <- sum(mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann$cis_or_trans == "trans")
bar_plot_hotspot_signal_df$n_signals_in_hotspot <- sum(mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann$hotspot)
bar_plot_hotspot_signal_df$n_signals_in_hotspot_cis <- sum(mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann$hotspot & mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann$cis_or_trans == "cis")
bar_plot_hotspot_signal_df$n_signals_in_hotspot_trans <- sum(mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann$hotspot & mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann$cis_or_trans == "trans")
bar_plot_hotspot_signal_df$n_signals_not_in_hotspot <- sum(!mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann$hotspot)
bar_plot_hotspot_signal_df$n_signals_not_in_hotspot_cis <- sum(!mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann$hotspot & mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann$cis_or_trans == "cis")
bar_plot_hotspot_signal_df$n_signals_not_in_hotspot_trans <- sum(!mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann$hotspot & mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann$cis_or_trans == "trans")
bar_plot_hotspot_signal_df

df_long <- bar_plot_hotspot_signal_df %>%
  pivot_longer(cols = everything(),
               names_to = "category",
               values_to = "count")

# Plot
ggplot(df_long, aes(x = category, y = count, fill = category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "", y = "Number of signals") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))


df <- data.frame(
  category = c("All signals", "All signals", "Hotspot signals", "Hotspot signals",  "Non hotspot signals", "Non hotspot signals"),
  type = c("cis", "trans", "cis", "trans", "cis", "trans"),
  count = c(1784, 6086, 287, 4136, 1497, 1950 )
)

# Plot stacked bars
ggplot(df, aes(x = category, y = count, fill = type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count),
            position = position_stack(vjust = 0.5),
            color = "black", size = 4) +
  theme_minimal() +
  labs(x = "", y = "Number of signals", fill = "") +
  theme(axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "top")

#### Adam's plot following Solène's code 202411_figure_panel_old_vs_new.R
# mapping<-fread("/home/solene.cadiou/basic_GWAS_protein/meta_results/MR/MR_instruments_selection/code/mapped_gene_file_GRCh37_21052025.txt")
#mapping <- read.delim("/exchange/healthds/pQTL/Reference_datasets_for_QC_proteomics/Cis_trans_mapping/somascan_tss_ncbi_grch37_version_20241212.txt", sep=",")
mapping <- read_csv("Supp_table_1_mapping_file.csv", skip = 2)
mapping$target<-mapping$SeqID
mapping$cis_end<-(mapping$TSS+500000)
mapping$cis_start<-(mapping$TSS-500000)
mapping$symbol <- mapping$Target_Name

lit <- mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann
minimap<-mapping[,c("target","symbol","chromosome","TSS","cis_start", "uniprot_new_in_somascan7k_vs5k")]
minimap$chromosome<-ifelse(minimap$chromosome=="X","23",minimap$chromosome)
minimap$chromosome<-ifelse(minimap$chromosome=="Y","23",minimap$chromosome)
table(minimap$chromosome)
minimap$chromosome<-as.numeric(minimap$chromosome)
table(table(minimap$target))
minimap<-minimap[!duplicated(minimap$target),]
minimap$study_id<-minimap$target
assoc2<-left_join(lit,minimap,by=c("phenotype_id"="study_id"))

data_cum <- assoc2|>
  group_by(chr) |>
  dplyr::summarise(max_bp = max(start)) |>
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) |>
  select(chr, bp_add)

data_map_cum <- assoc2|>
  group_by(chromosome) |>
  dplyr::summarise(max_bp_map = max(cis_start)) |>
  mutate(bp_add_map = lag(cumsum(as.numeric(max_bp_map)), default = 0)) |>
  select(chromosome, bp_add_map)

gwas_data <- assoc2 |>
  inner_join(data_cum, by = "chr") |>
  mutate(bp_cum = start+ bp_add)
gwas_data <- gwas_data |>
  inner_join(data_map_cum, by = "chromosome") |>
  mutate(bp_cum_map = cis_start+ bp_add_map)
gwas_data$MLOG10P<-as.numeric(gwas_data$MLOG10P)


df2=data.frame(chr=1:22,center=NA,stringsAsFactors=F)

for(chr in 1:22){
  
  start=min(gwas_data$bp_cum[gwas_data$chr==chr],na.rm=T)
  
  end=max(gwas_data$bp_cum[gwas_data$chr==chr],na.rm=T)
  
  center=mean(c(start,end))
  
  df2[df2$chr==chr,"center"]=center
}  


df3=data.frame(chr=1:22,center=NA,stringsAsFactors=F)

for(chr in 1:22){
  
  start=min(gwas_data$bp_cum_map[gwas_data$chromosome==chr],na.rm=T)
  
  end=max(gwas_data$bp_cum_map[gwas_data$chromosome==chr],na.rm=T)
  
  center=mean(c(start,end))
  
  df3[df3$chr==chr,"center"]=center
  
}

p<-ggplot(gwas_data, aes(x=bp_cum,y=bp_cum_map
                         # label=outlier
))+geom_point(colour = as.factor(assoc2$chr))+
  scale_x_continuous(breaks=df2$center,labels=c(1:22))+
  scale_y_continuous(breaks=df3$center,labels=c(1:22))+
  theme(axis.text.x = element_text(angle = 40))+ xlab("pQTL positions") + ylab("Coding gene position")
p

gwas_data$new_uniprot <- ifelse(gwas_data$uniprot_new_in_somascan7k_vs5k, "New in v7", "Present in v5" )
gwas_data$new_uniprot<-factor(gwas_data$new_uniprot,levels=c("New in v7","Present in v5"))

p<-ggplot(gwas_data, aes(x=bp_cum,y=bp_cum_map,col=new_uniprot
                         # label=outlier
))+geom_point(alpha=0.5)+
  guides(size = 'none')+
  scale_x_continuous(breaks=df2$center,labels=c(1:22))+
  scale_y_continuous(breaks=df3$center,labels=c(1:22))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 40),legend.position="bottom")+ 
  xlab("pQTL positions") + 
  ylab("Coding gene position")+ 
  scale_color_manual(values = c("#69b3a2","#404080"),
                     labels=c("New in v7","Present in v5"))+
  guides(alpha="none")+labs(color="")

p


# top_symbol_per_hotspot <- gwas_data %>%
#   filter(hotspot == TRUE & cis_or_trans == "cis") %>%
#   separate_rows(symbol.x, sep = "\\|") %>%
#   mutate(symbol.x = trimws(symbol.x)) %>%
#   filter(!is.na(symbol.x), symbol.x != "") %>%
#   group_by(chr, full_hotspot_gene_window, symbol.x) %>%
#   summarise(
#     n = n(),
#     max_strength = suppressWarnings(max(MLOG10P, na.rm = TRUE)),
#     .groups = "drop"
#   ) %>% 
#   mutate(max_strength = ifelse(is.finite(max_strength), max_strength, NA_real_)) %>%
#   group_by(chr, full_hotspot_gene_window) %>%
#   arrange(desc(n), desc(max_strength), symbol.x) %>%
#   slice(1) %>%
#   ungroup() %>%
#   rename(label = symbol.x)

top_symbol_per_hotspot <- read_csv("Hotspot_Summary_Final.csv")
colnames(top_symbol_per_hotspot)[1] <- "label"
  
hotspot_label_position <- gwas_data %>%
  filter(hotspot == TRUE) %>%
  group_by(chr, full_hotspot_gene_window) %>%
  summarise(
    x = median(bp_cum, na.rm = TRUE),
    y = max(bp_cum_map, na.rm = TRUE),
    .groups = "drop"
  )

labels_df <- top_symbol_per_hotspot %>%
  left_join(hotspot_label_position, by = c("chr", "full_hotspot_gene_window")) %>%
  mutate(label = str_trim(label))



p_annot <- p +
  expand_limits(y = max(gwas_data$bp_cum_map, na.rm = TRUE) * 1.05) +
  geom_text_repel(
    data = labels_df,
    inherit.aes = FALSE,
    aes(x = x, y = y, label = label),
    size = 3,
    max.overlaps = Inf,
    ylim = c(max(gwas_data$bp_cum_map, na.rm = TRUE) * 1.02, NA),
    box.padding = 0.3,
    point.padding = 0.2,
    min.segment.length = 0
  )

p_annot

#### Rest of Figure 2 update

Supp_table_2_LB_results <- read_csv("Supp_table_2_LB_results.csv", skip = 3)

colnames(Supp_table_2_LB_results)
df <- data.frame(
  category = c(rep("All signals",2),rep( "New signals",2), rep("Rep signals", 2),
               rep("Cis", 2), rep("Trans", 2), rep("New cis", 2), rep("New trans", 2),
               rep("Signals in hotspots",2), rep("Cis signals in hotspots",2),
               rep("Trans signals in hotspots",2),rep("Heterogeneus signals",2)),
  type = rep(c("Present in v5", "New in v7"), 11),
  count = c(
    sum(!Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(is.na(Supp_table_2_LB_results$unip_matching_study) & !Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(is.na(Supp_table_2_LB_results$unip_matching_study) & Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(!is.na(Supp_table_2_LB_results$unip_matching_study) & !Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(!is.na(Supp_table_2_LB_results$unip_matching_study) & Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$cis_or_trans == "cis" &!Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$cis_or_trans == "cis" & Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$cis_or_trans == "trans" &!Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$cis_or_trans == "trans" & Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(is.na(Supp_table_2_LB_results$unip_matching_study) & Supp_table_2_LB_results$cis_or_trans == "cis" &!Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(is.na(Supp_table_2_LB_results$unip_matching_study) & Supp_table_2_LB_results$cis_or_trans == "cis" & Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(is.na(Supp_table_2_LB_results$unip_matching_study) & Supp_table_2_LB_results$cis_or_trans == "trans" &!Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(is.na(Supp_table_2_LB_results$unip_matching_study) & Supp_table_2_LB_results$cis_or_trans == "trans" & Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$hotspot &!Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$hotspot & Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$hotspot & Supp_table_2_LB_results$cis_or_trans == "cis" &!Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$hotspot & Supp_table_2_LB_results$cis_or_trans == "cis" & Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$hotspot & Supp_table_2_LB_results$cis_or_trans == "trans" &!Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$hotspot & Supp_table_2_LB_results$cis_or_trans == "trans" & Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$HETISQ>90  &!Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$HETISQ>90  & Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k)
  )
)

df <- df %>%
  mutate(
    category = factor(category, levels = unique(category))  # optional, keeps x order too
  )
df_totals <- df %>%
  group_by(category) %>%
  summarise(total = sum(count), .groups = "drop")

thr <- 0.10
labs_inside <- df %>% filter(count >= 100)


# Plot stacked bars
ggplot(df, aes(x = category, y = count, fill = type)) +
  geom_bar(stat = "identity", alpha = 0.5, aes(fill = type, color = type),width = 0.7) +
  geom_text(data= labs_inside,
            aes(label = count),
            position = position_stack(vjust = 0.5),
            color = "black", size = 4) +
  geom_text(
    data = df_totals,
    inherit.aes = FALSE,
    aes(x = category, y = total, label = total),
    vjust = -0.6, fontface = "bold", size = 4
  ) +
  theme_minimal() +
  labs(x = "", y = "Number of signals", fill = "") +
  scale_fill_manual(values = c("Present in v5" = "#404080","New in v7" = "#69b3a2")) +
  scale_color_manual(values = c("Present in v5" = "#404080","New in v7" = "#69b3a2"), , guide = "none") +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 12),
        legend.position = "top")
