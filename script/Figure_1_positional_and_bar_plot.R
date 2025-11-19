setwd("/group/diangelantonio/users/alessia_mapelli/pQTL/INTERVAL/Hotspots and lonespots/New_coloc_results")
rm(list=ls(all=TRUE))

library(readr)
library(dplyr)
library(scales)
library(ggrepel)
library(data.table)
library(tidyverse)
library(readxl)

mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann <- read_excel("supplementary_table_3.xlsx", sheet = 2)
colnames(mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann) <- mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann[3,]
lit <- mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann[-c(1,2,3),]

#### Adam's plot following Solène's code 202411_figure_panel_old_vs_new.R
mapping <- read_excel("supplementary_table_2.xlsx", sheet = 2)
mapping$target<-mapping$SeqId
mapping$cis_end<-(mapping$TSS+500000)
mapping$cis_start<-(mapping$TSS-500000)
mapping$symbol <- mapping$Target_Name

minimap<-mapping[,c("target","symbol","chromosome","TSS","cis_start")]
minimap$chromosome<-ifelse(minimap$chromosome=="X","23",minimap$chromosome)
minimap$chromosome<-ifelse(minimap$chromosome=="Y","23",minimap$chromosome)
table(minimap$chromosome)
minimap$chromosome<-as.numeric(minimap$chromosome)
table(table(minimap$target))
minimap<-minimap[!duplicated(minimap$target),]
minimap$study_id<-minimap$target
lit$START <- unlist(str_split(lit$locus_START_END_37[1], "_"))[2]
lit$START <- as.numeric(lit$START)
lit$CHR <- as.numeric(lit$CHR)
assoc2<-left_join(lit,minimap,by=c("SeqID"="study_id"))

data_cum <- assoc2|>
  group_by(CHR) |>
  dplyr::summarise(max_bp = max(START)) |>
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) |>
  select(CHR, bp_add)

data_map_cum <- assoc2|>
  group_by(chromosome) |>
  dplyr::summarise(max_bp_map = max(cis_start)) |>
  mutate(bp_add_map = lag(cumsum(as.numeric(max_bp_map)), default = 0)) |>
  select(chromosome, bp_add_map)

gwas_data <- assoc2 |>
  inner_join(data_cum, by = "CHR") |>
  mutate(bp_cum = START+ bp_add)
gwas_data <- gwas_data |>
  inner_join(data_map_cum, by = "chromosome") |>
  mutate(bp_cum_map = cis_start+ bp_add_map)
gwas_data$MLOG10P<-as.numeric(gwas_data$MLOG10P)


df2=data.frame(chr=1:22,center=NA,stringsAsFactors=F)

for(chr in 1:22){
  
  start=min(gwas_data$bp_cum[gwas_data$CHR==chr],na.rm=T)
  
  end=max(gwas_data$bp_cum[gwas_data$CHR==chr],na.rm=T)
  
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
))+geom_point(colour = as.factor(assoc2$CHR))+
  scale_x_continuous(breaks=df2$center,labels=c(1:22))+
  scale_y_continuous(breaks=df3$center,labels=c(1:22))+
  theme(axis.text.x = element_text(angle = 40))+ xlab("pQTL positions") + ylab("Coding gene position")
p

gwas_data$new_uniprot <- ifelse(gwas_data$uniprot_new_in_somascan7k_vs5k, "New in 7k", "Present in 5k" )
gwas_data$new_uniprot<-factor(gwas_data$new_uniprot,levels=c("New in 7k","Present in 5k"))

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
  group_by(CHR, full_hotspot_gene_window) %>%
  summarise(
    x = median(bp_cum, na.rm = TRUE),
    y = max(bp_cum_map, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(Gene_window = full_hotspot_gene_window)

labels_df <- top_symbol_per_hotspot %>%
  left_join(hotspot_label_position, by = c("CHR", "Gene_window")) %>%
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

p <- ggplot(gwas_data, aes(x = bp_cum, y = bp_cum_map, col = new_uniprot)) +
  geom_point(alpha = 0.5) +
  scale_x_continuous(breaks = df2$center, labels = 1:22) +
  scale_y_continuous(breaks = df3$center, labels = 1:22) +
  scale_color_manual(
    values = c("#69b3a2", "#404080"),
    labels = c("New in v7", "Present in v5")
  ) +
  theme_bw() +
  theme(
    legend.position = "none",                        # no legend
    axis.text.x  = element_text(angle = 40, size = 12),
    axis.text.y  = element_text(size = 12),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13)
  ) +
  xlab("pQTL positions") +
  ylab("Coding gene position") +
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

p

#### Rest of Figure 2 update

Supp_table_2_LB_results <- lit
str(Supp_table_2_LB_results)
Supp_table_2_LB_results$hotspot <- as.logical(Supp_table_2_LB_results$hotspot )
Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k <- as.logical(Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k )
Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k <- as.logical(Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k )

colnames(Supp_table_2_LB_results)
head(Supp_table_2_LB_results)
df <- data.frame(
  category = c(rep("All signals",2),rep( "New signals",2), rep("Rep signals", 2),
               rep("Cis", 2), rep("Trans", 2), rep("New cis", 2), rep("New trans", 2),
               rep("Signals in hotspots",2), rep("Cis signals in hotspots",2),
               rep("Trans signals in hotspots",2),rep("Heterogeneus signals",2)),
  type = rep(c("Present in 5k", "New in 7k"), 11),
  count = c(
    sum(!Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$unip_matching_study=="NA" & !Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$unip_matching_study=="NA" & Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$unip_matching_study !="NA" & !Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$unip_matching_study !="NA" & Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$cis_or_trans == "cis" &!Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$cis_or_trans == "cis" & Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$cis_or_trans == "trans" &!Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$cis_or_trans == "trans" & Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$unip_matching_study=="NA" & Supp_table_2_LB_results$cis_or_trans == "cis" &!Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$unip_matching_study=="NA" & Supp_table_2_LB_results$cis_or_trans == "cis" & Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$unip_matching_study=="NA" & Supp_table_2_LB_results$cis_or_trans == "trans" &!Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
    sum(Supp_table_2_LB_results$unip_matching_study=="NA" & Supp_table_2_LB_results$cis_or_trans == "trans" & Supp_table_2_LB_results$uniprot_new_in_somascan7k_vs5k),
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

df_reordered <- df[c(1,2,7,8,9,10,3,4,11,12,13,14,15,16,17,18,19,20,21,22),]

df_reordered <- df_reordered %>%
  mutate(
    category = factor(category, levels = unique(category))  # optional, keeps x order too
  )
df_totals <- df_reordered %>%
  group_by(category) %>%
  summarise(total = sum(count), .groups = "drop")

thr <- 0.10
labs_inside <- df_reordered %>% filter(count >= 100)

df_reordered$group <- c(rep(1,6), rep(2,6), rep(3,6), rep(4,2))

# Plot stacked bars
ggplot(df_reordered, aes(x = category, y = count, fill = type)) +
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
  scale_fill_manual(values = c("Present in 5k" = "#404080","New in 7k" = "#69b3a2")) +
  scale_color_manual(values = c("Present in 5k" = "#404080","New in 7k" = "#69b3a2"), , guide = "none") +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 12),
        legend.position = "top")





df_reordered <- df_reordered %>%
  mutate(
    category = factor(category, levels = unique(category)),
    group    = factor(group, levels = sort(unique(group)))
  )

# recompute totals & inside labels AFTER adding `group`
df_totals <- df_reordered %>%
  group_by(group, category) %>%
  summarise(total = sum(count), .groups = "drop")

labs_inside <- df_reordered %>% filter(count >= 100)

ggplot(df_reordered, aes(x = category, y = count, fill = type)) +
  geom_bar(stat = "identity", alpha = 0.5, aes(color = type), width = 0.7) +
  # inside labels
  geom_text(data = labs_inside, aes(label = count),
            position = position_stack(vjust = 0.5),
            color = "black", size = 4, inherit.aes = TRUE) +
  # totals on top
  geom_text(data = df_totals,
            aes(x = category, y = total, label = total),
            vjust = -0.6, fontface = "bold", size = 4,
            inherit.aes = FALSE) +
  # colors
  scale_fill_manual(values = c("Present in 5k" = "#404080", "New in 7k" = "#69b3a2"), name = "") +
  scale_color_manual(values = c("Present in 5k" = "#404080", "New in 7k" = "#69b3a2"), guide = "none") +
  # create SPACE between groups
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  labs(x = "", y = "Number of signals") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(size = 12, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 12),
    legend.position = "top",
    strip.text.x = element_blank(),           # hide facet strip labels
    strip.background = element_blank(),
    panel.spacing.x = unit(1, "cm")           # <-- space between groups
  )

ggplot(df_reordered, aes(x = category, y = count, fill = type)) +
  geom_bar(stat = "identity", alpha = 0.5, aes(color = type), width = 0.7) +
  # inside labels
  geom_text(
    data = labs_inside, aes(label = count),
    position = position_stack(vjust = 0.5),
    color = "black", size = 4, inherit.aes = TRUE
  ) +
  # totals on top
  geom_text(
    data = df_totals,
    aes(x = category, y = total, label = total),
    vjust = -0.6, fontface = "bold", size = 4,
    inherit.aes = FALSE
  ) +
  # colors
  scale_fill_manual(
    values = c("Present in 5k" = "#404080", "New in 7k" = "#69b3a2"),
    name = ""
  ) +
  scale_color_manual(
    values = c("Present in 5k" = "#404080", "New in 7k" = "#69b3a2"),
    guide = "none"
  ) +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  labs(x = "", y = "Number of signals") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 12),
    axis.title.y = element_text(size = 13),
    legend.position = "bottom",
    legend.text = element_text(size = 14, face = "bold"),
    strip.text.x = element_blank(),
    strip.background = element_blank(),
    panel.spacing.x = unit(1, "cm")
  )

