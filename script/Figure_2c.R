setwd("/group/diangelantonio/users/alessia_mapelli/pQTL/INTERVAL/Hotspots and lonespots/New_coloc_results/Coloc_networks_top_cond")
rm(list=ls(all=TRUE))

if (!requireNamespace("ggraph", quietly = TRUE)) {
  install.packages("ggraph")
}

library(readr)
library(dplyr)
library(igraph)
library(ggplot2)
library(scales)
library(ggrepel)
library(Matrix)
library(ggraph)
library(tidygraph)
library(data.table)


#################################################
## 0. USER DEFINED PARAMETERS (MODIFY THIS PART)
################################################
colocalization_results_filtered_path ="/group/diangelantonio/users/alessia_mapelli/pQTL/INTERVAL/Hotspots and lonespots/New_coloc_results/Colocalization_results_filtered.csv"
hotspots_df_path = "/group/diangelantonio/users/alessia_mapelli/pQTL/INTERVAL/Hotspots and lonespots/New_coloc_results/Hotspots_df.csv"
type_of_clustering = "connected" # either "greedy", "connected"
##############################################

#################################################
## Include parameter info in the logs
cat("Coloc results: ", colocalization_results_filtered_path,"\n")
cat("Hotspost df: ", hotspots_df_path,"\n")
#################################################


#############################################
###### 1. Upload all the data needed
#############################################
args <- commandArgs(trailingOnly = TRUE)
k = 17
colocalization_results_filtered <- fread(colocalization_results_filtered_path)
hotspot_df <- read_csv(hotspots_df_path)[,-1]

LB <- read_excel("/group/diangelantonio/users/alessia_mapelli/pQTL/INTERVAL/Hotspots and lonespots/New_coloc_results/supplementary_table_3.xlsx", sheet = 2)
colnames(LB) <- LB[3,]
LB <- LB[-c(1,2,3),]
LB_in_hotsposts <- LB %>%
  filter(hotspot== TRUE)

LB_in_hotsposts$CHR <- as.character(as.numeric(LB_in_hotsposts$CHR))


##############################################
##### 2.Computation
##########################################
cat("Processing row", k, "\n")
# Create a graph that has as node the seqid linked to a specific hotspot
sel_hotspot <- as.character(hotspot_df[k, 1])
chr_hotspot <- as.character(hotspot_df[k, 2])
LB_in_hotspost <- LB_in_hotsposts %>%
  filter(full_hotspot_gene_window == sel_hotspot & CHR == chr_hotspot)
colocalization_results_filtered_in_hotspot <- colocalization_results_filtered %>%
  filter(full_hotspot_gene_window_locus_a == sel_hotspot & chr_locus_a == chr_hotspot &
           full_hotspot_gene_window_locus_b == sel_hotspot & chr_locus_b == chr_hotspot)

cat("Starting computation of the adjacency matrix \n")

colocalization_results_filtered_in_hotspot <- colocalization_results_filtered_in_hotspot %>%
  mutate(
    node_name_a = paste(trait_a, top_cond_a, sep=":"),
    node_name_b = paste(trait_b, top_cond_b, sep=":")
  )
nodes_name <- unique(c(colocalization_results_filtered_in_hotspot$node_name_a, colocalization_results_filtered_in_hotspot$node_name_b))
# Keep only the nodes that at least colocalize with another signal

node_index <- setNames(seq_along(nodes_name), nodes_name)
row_idx <- node_index[colocalization_results_filtered_in_hotspot$node_name_a]
col_idx <- node_index[colocalization_results_filtered_in_hotspot$node_name_b]

adj_matrix <- sparseMatrix(
  i = c(row_idx, col_idx),
  j = c(col_idx, row_idx),
  x = 1,
  dims = c(length(nodes_name), length(nodes_name)),
  dimnames = list(nodes_name, nodes_name)
)
adj_matrix@x[] <- 1

cat("Computation of the adjacency matrix finished \n")

g_est <- graph_from_adjacency_matrix(adj_matrix , mode = "undirected", diag = F,add.colnames = NA, add.rownames = NULL )
summary_adj_matrix <- summary(adj_matrix)
summary_adj_matrix$i_Name <- rownames(adj_matrix)[summary_adj_matrix$i]
summary_adj_matrix$j_Name <- colnames(adj_matrix)[summary_adj_matrix$j]
nodes_plot_label <- match(V(g_est)$name, nodes_name)
V(g_est)$plot_label <- nodes_plot_label
nod_legend <- data.frame(nodes_names = V(g_est)$name, label= match(V(g_est)$name, nodes_name))

cat("Starting community computation \n")
cfg <- cluster_fast_greedy(g_est)

# AVOID COMPUTING cluster_edge_betweenness TO REDUCE COMPUTATIONAL TIME
#ceb <- cluster_edge_betweenness(g_est)

comp <- components(g_est)
cat("Community computation done \n")

#Annotate results
communities <-data.frame(nodes = V(g_est)$name,
                         group_conncetd_component = comp$membership,
                         #group_between_clustering = ceb$membership,
                         group_greedy_clustering = cfg$membership )

communities$phenotype_id <- rep(NA, nrow(communities))
communities$chr <- rep(NA, nrow(communities))
communities$start <- rep(NA, nrow(communities))
communities$end <- rep(NA, nrow(communities))

for(i in 1: nrow(communities)){
  node_name <- communities$nodes[i]
  temp_a <- colocalization_results_filtered_in_hotspot %>% 
    filter(node_name_a == node_name) %>%
    select(trait_a, chr_locus_a, start_locus_a, end_locus_a ) %>% 
    rename(phenotype_id = trait_a, chr=chr_locus_a,
           start= start_locus_a, end = end_locus_a)
  temp_b <- colocalization_results_filtered_in_hotspot %>% 
    filter(node_name_b == node_name) %>%
    select(trait_b, chr_locus_b, start_locus_b, end_locus_b ) %>% 
    rename(phenotype_id = trait_b, chr=chr_locus_b,
           start= start_locus_b, end = end_locus_b)
  temp <- rbind(temp_a,temp_b)
  if(nrow(temp)>0){
    communities$phenotype_id[i] <- temp$phenotype_id[1]
    communities$chr[i] <- temp$chr[1]
    communities$start[i] <- temp$start[1]
    communities$end[i] <- temp$end[1]
  }
}

if(type_of_clustering == "connected"){
  communities$group = communities$group_conncetd_component
  V(g_est)$group <- communities$group_conncetd_component
  #} else if(type_of_clustering == "between") {
  #  communities$group = communities$group_between_clustering
  #  V(g_est)$group <- communities$group_between_clustering
} else {
  communities$group = communities$group_greedy_clustering
  V(g_est)$group <- communities$group_greedy_clustering
}

communities$CHR <- as.character(communities$chr)
communities$locus_START_END_37 <- paste(paste("chr",communities$CHR, sep=""), communities$start, communities$end, sep="_" )
communities$SeqID <- communities$phenotype_id

communities <-  full_join(LB_in_hotspost, communities, by= c("CHR", "locus_START_END_37", "SeqID"))


communities$community_type <-"all_trans"
for (group_id in unique(communities$group)) {
  temp <- communities %>% filter(group == group_id)
  
  # If the community has a single node
  if (nrow(temp) == 1) {
    communities$community_type[communities$group == group_id] <- temp$cis_or_trans
  } else {
    # Count "cis" and "trans" nodes
    cis_nodes <- sum(temp$cis_or_trans == "cis")
    trans_nodes <- sum(temp$cis_or_trans == "trans")
    
    # Assign community types based on the conditions
    communities$community_type[communities$group == group_id] <- dplyr::case_when(
      cis_nodes == 0 ~ "all_trans",
      cis_nodes == 1 & trans_nodes >= 1 ~ "one_cis_multi_trans",
      cis_nodes > 1 & trans_nodes > 1 ~ "multi_cis_multi_trans",
      cis_nodes > 1 & trans_nodes == 0 ~ "multi_cis",
      TRUE ~ "unknown"  # Fallback in case of unexpected conditions
    )
  }
}

cat("Saving the plots \n")
foldname <- paste(chr_hotspot,"_", sel_hotspot,"_20251119", sep ="")
dir.create(foldname)

num_colors <- length((unique(V(g_est)$group)))
colrs <- scales::hue_pal()(num_colors)

weight.community=function(row,weigth.within,weight.between){
  if(as.numeric(V(g_est)$group[V(g_est)$name == row[1]])==as.numeric(V(g_est)$group[V(g_est)$name == row[2]])){
    weight=weigth.within
  }else{
    weight=weight.between
  }
  return(weight)
}

E(g_est)$weight=apply(get.edgelist(g_est),1,weight.community,2,1)
V(g_est)$cis_trans <- NA
for(i in  1: length(V(g_est)$name)){
  temp_name <- V(g_est)$name[i]
  temp <- communities %>% filter(nodes == temp_name)
  if(nrow(temp)>0){
    V(g_est)$cis_trans[i] <- temp$cis_or_trans[1]}
}

mapping <- read_excel("/group/diangelantonio/users/alessia_mapelli/pQTL/INTERVAL/Hotspots and lonespots/New_coloc_results/supplementary_table_2.xlsx", sheet = 2)

mapping <- mapping[, c("SeqId","chromosome","TSS", "Target_Name", "UniProt_ID") ]
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

library(stringr)
LB$START <- as.numeric(word(LB$locus_START_END_37, 2, sep = "_"))
LB$START <- as.numeric(LB$START)
LB$CHR <- as.numeric(LB$CHR)
assoc2<-left_join(LB,minimap,by=c("SeqID"="study_id"))

data_map_cum <- assoc2|>
  group_by(chromosome) |>
  dplyr::summarise(max_bp_map = max(cis_start)) |>
  mutate(bp_add_map = lag(cumsum(as.numeric(max_bp_map)), default = 0)) |>
  select(chromosome, bp_add_map)

gwas_data <- assoc2 |>
  inner_join(data_map_cum, by = "chromosome") |>
  mutate(bp_cum_map = cis_start+ bp_add_map)
gwas_data$MLOG10P<-as.numeric(gwas_data$MLOG10P)

df3=data.frame(chr=1:22,center=NA,stringsAsFactors=F)

for(chr in 1:22){
  
  start=min(gwas_data$bp_cum_map[gwas_data$chromosome==chr],na.rm=T)
  
  end=max(gwas_data$bp_cum_map[gwas_data$chromosome==chr],na.rm=T)
  
  center=mean(c(start,end))
  
  df3[df3$chr==chr,"center"]=center
  
}

communities <- left_join(communities, minimap, by=c("SeqID"="study_id"))

communities <- communities %>%
  inner_join(data_map_cum, by = "chromosome") %>%
  mutate(bp_cum_map = cis_start+ bp_add_map)

library(stringr)
communities_adam_plot <- communities %>%
  select(SNPID, bp_cum_map, group, cis_or_trans, HARMONIZED_GENE_NAME) %>%
  mutate(POS = as.numeric(stringr::word(SNPID, 2, sep = ":")))

num_colors <- length(unique(V(g_est)$group))
palette <- scales::hue_pal()(num_colors)

tg <- as_tbl_graph(g_est)
layout <- create_layout(tg, layout = 'graphopt') 


communities_adam_plot <- communities_adam_plot %>%
  mutate(
    group = as.character(group),
    group = ifelse(is.na(group), "Not connected", group),
    group = factor(group)
  ) %>%
  rename(Communities=group)

library(ggplot2)
library(ggrepel)
library(scales)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


ggplot(
  communities_adam_plot,
  aes(x = POS, y = bp_cum_map, color = Communities, shape = cis_or_trans)
) +
  geom_point(size = 3) +
  geom_text_repel(
    data = subset(communities_adam_plot, cis_or_trans == "cis"),
    aes(label = HARMONIZED_GENE_NAME),
    size = 4,
    max.overlaps = 70
  ) +
  scale_shape_manual(
    values = c("cis" = 1, "trans" = 15),
    name = "Signal Type"   # <- Legend title for shapes
  ) +
  scale_x_continuous(labels = label_number()) +
  scale_y_continuous(breaks=df3$center,labels=c(1:22))+
  labs(
    x = "pQTL positions",
    y = "Coding gene position"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
    # legend.position = "none",
    # panel.grid = element_blank()
  )


plot_name <- paste(foldname, "/Adams_plot_",type_of_clustering,"_new_1.png", sep="")
ggsave(plot_name, width = 25, height = 20, units = "cm")

tg <- as_tbl_graph(g_est)
layout <- create_layout(tg, layout = 'graphopt') 
p <- ggraph(layout) +
  geom_edge_link(edge_width = 0.5, alpha = 0.2, color = "gray") +
  geom_node_point(aes(
    color = as.factor(V(g_est)$group),
    # shape = as.factor(V(g_est)$cis_trans),
    # fill = ifelse(V(g_est)$cis_trans == "cis", as.factor(V(g_est)$group), "white")
  ),shape = 15, size = 4) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = c(palette, white = "white")) +
  theme_void() +                                 # Remove background, grid, axes
  theme(legend.position = "none")                # Remove all legendsp
p
ggsave(p, filename=paste(foldname, "/HighQuality_Network_",type_of_clustering, "_new1.png", sep=""), width=8, height=8, dpi=300)

library(ggplot2)
library(ggrepel)
library(scales)

mb_format <- function(x) sprintf("%.0f Mb", x / 1e6)

labels_df <- communities_adam_plot %>%
  filter(
    cis_or_trans == "cis",
    !is.na(start),
    !is.na(bp_cum_map),
    !is.na(HARMONIZED_GENE_NAME)
  ) %>%
  arrange(bp_cum_map) %>%                     # "first signal" = smallest bp_cum_map
  distinct(HARMONIZED_GENE_NAME, .keep_all = TRUE)


ggplot(
  communities_adam_plot,
  aes(x = POS, y = bp_cum_map, color = Communities, shape = cis_or_trans)
) +
  geom_point(size = 3) +
  geom_text_repel(
    data = labels_df,
    aes(label = HARMONIZED_GENE_NAME),
    size = 4,
    max.overlaps = 70,
    nudge_y = 0.07 * diff(range(communities_adam_plot$bp_cum_map, na.rm = TRUE)),
    direction = "y",
    vjust = 0
  ) +
  scale_shape_manual(
    values = c("cis" = 16, "trans" = 15),
    name = "Signal Type"
  ) +
  scale_x_continuous(labels = mb_format) +
  scale_y_continuous(breaks=df3$center,labels=c(1:22))+
  labs(
    x = "\npQTL positions (Mb)",
    y = "\nCoding gene position (Mb)\n"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.position = "bottom", 
    legend.key.size = unit(1, "cm"),          # makes keys (symbols) larger
    legend.text = element_text(size = 14, face = "bold"),    # makes legend label text bigger
    legend.title = element_text(size = 14),    # makes legend title bigger
    legend.spacing.x = unit(2, "cm")
  )

plot_name <- paste(foldname, "/Adams_plot_",type_of_clustering,"_new_6.png", sep="")
ggsave(plot_name, width = 30, height = 20, units = "cm")


community_names <- communities_adam_plot %>%
  count(Communities, SNPID, name = "n") %>%        # count how often each SNPID appears per group
  group_by(Communities) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>% # keep SNPID with highest count in each group
  ungroup()

communities_res <- communities_adam_plot %>%
  mutate(
    Communities_recode = case_when(
      Communities== "5" ~ "12:7242740:C:G",  
      Communities== "6" ~ "12:7170336:A:G",  
      Communities== "8" ~ "12:7242204:C:T",  
      Communities == "Not connected"  ~ "Not colocalized signals",
      TRUE ~ "Others"                                   
    ),
    Communities_recode = factor(
      Communities_recode,
      levels = c("12:7242740:C:G","12:7170336:A:G","12:7242204:C:T","Others","Not colocalized signals")
    )
  )


community_colors <- c(
  "12:7242740:C:G" = gg_color_hue(3)[1],           # blue
  "12:7170336:A:G" = gg_color_hue(3)[2],           # green
  "12:7242204:C:T" = gg_color_hue(3)[3],           # red
  "Others" = "darkgrey",
  "Not colocalized signals" = "lightgrey"
)

labels_df <- communities_res %>%
  filter(
    cis_or_trans == "cis",
    !is.na(start),
    !is.na(bp_cum_map),
    !is.na(HARMONIZED_GENE_NAME)
  ) %>%
  arrange(bp_cum_map) %>%                     # "first signal" = smallest bp_cum_map
  distinct(HARMONIZED_GENE_NAME, .keep_all = TRUE)


library(ggplot2)
library(ggrepel)
library(scales)

mb_format <- function(x) sprintf("%.0f Mb", x / 1e6)

ggplot(
  communities_res,
  aes(x = POS, y = bp_cum_map, color = Communities_recode, shape = cis_or_trans)
) +
  geom_point(size = 3) +
  
  geom_text_repel(
    data = labels_df,
    aes(label = HARMONIZED_GENE_NAME),
    size = 4,
    max.overlaps = 70,
    nudge_y = 0.07 * diff(range(communities_res$bp_cum_map, na.rm = TRUE)),
    direction = "y",
    vjust = 0
  )+
  scale_color_manual(values = community_colors, name = "Communities") +
  
  scale_shape_manual(values = c("cis" = 16, "trans" = 15), name = "Signal Type") +
  
  scale_x_continuous(labels = mb_format) +
  scale_y_continuous(breaks=df3$center,labels=c(1:22))+
  
  labs(
    x = "\npQTL positions (Mb)",
    y = "\nCoding gene position (Mb)\n"
  ) +
  
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14),
    legend.spacing.x = unit(2, "cm")
  )

ggsave(
  paste0(foldname, "/Adams_plot_", type_of_clustering, "_new_7.png"),
  width = 30, height = 20, units = "cm"
)


communities_res <- communities_res %>%
  filter( Communities_recode %in% c("12:7242740:C:G","12:7170336:A:G","12:7242204:C:T"))

labels_df <- communities_res %>%
  filter(
    cis_or_trans == "cis",
    !is.na(start),
    !is.na(bp_cum_map),
    !is.na(HARMONIZED_GENE_NAME)
  ) %>%
  arrange(bp_cum_map) %>%                     # "first signal" = smallest bp_cum_map
  distinct(HARMONIZED_GENE_NAME, .keep_all = TRUE)

labels_df_trans_HG <- communities %>%
  filter(Target_Name %in% c("C4a", "C1q"))

labels_df_trans <- communities_res %>%
  filter(HARMONIZED_GENE_NAME %in% labels_df_trans_HG$HARMONIZED_GENE_NAME)

labels_df <- rbind(labels_df,labels_df_trans )
library(ggplot2)
library(ggrepel)
library(scales)

mb_format <- function(x) sprintf("%.2f Mb", x / 1e6)

ggplot(
  communities_res,
  aes(x = POS, y = bp_cum_map, color = Communities_recode, shape = cis_or_trans)
) +
  geom_point(size = 3) +
  
  geom_text_repel(
    data = labels_df,
    aes(label = HARMONIZED_GENE_NAME),
    size = 4,
    max.overlaps = 70,
    nudge_y = 0.07 * diff(range(communities_res$bp_cum_map, na.rm = TRUE)),
    direction = "y",
    vjust = 0
  )+
  scale_color_manual(values = community_colors, name = "Communities") +
  
  scale_shape_manual(values = c("cis" = 16, "trans" = 15), name = "Signal Type") +
  
  scale_x_continuous(labels = mb_format) +
  scale_y_continuous(breaks=df3$center,labels=c(1:22))+
  
  labs(
    x = "\npQTL positions (Mb)",
    y = "\nCoding gene position (Mb)\n"
  ) +
  
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.key.size = unit(1, "cm"),
    legend.box = "vertical",
    legend.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14),
    legend.spacing.x = unit(2, "cm")
  )

ggsave(
  paste0(foldname, "/Adams_plot_", type_of_clustering, "_new_8.png"),
  width = 30, height = 20, units = "cm"
)

# Plot just the subgraph on groups 5, 6, 8
keep_groups <- c(5, 6, 8)
nodes_keep  <- V(g_est)[group %in% keep_groups]

cols <- gg_color_hue(3)

community_colors <- c(
  "5" = cols[1],
  "6" = cols[2],
  "8" = cols[3]
)

g_sub <- induced_subgraph(g_est, vids = nodes_keep)
tg_sub  <- as_tbl_graph(g_sub)
layout_sub <- create_layout(tg_sub, layout = "graphopt",charge = 0.0000001, spring.length = 0.01,spring.constant = 0.1)
p <- ggraph(layout_sub) +
  geom_edge_link(edge_width = 0.5, alpha = 0.2, color = "gray") +
  geom_node_point(
    aes(color = as.factor(group)),
    shape = 15,
    size = 4
  ) +
  scale_color_manual(values = community_colors, drop = FALSE)+
  theme_void() +
  theme(legend.position = "none")

p
ggsave(p, filename=paste(foldname, "/HighQuality_Network_",type_of_clustering, "_red_v1.png", sep=""), width=8, height=8, dpi=300)


community_colors <- setNames(
  rep("darkgray", 18),
  as.character(1:18)
)
community_colors[c("5", "6", "8")] <- cols[1:3]


tg <- as_tbl_graph(g_est)
layout <- create_layout(tg, layout = 'graphopt',charge = 0.00001, spring.length = 0.01,spring.constant = 0.25) 
p <- ggraph(layout) +
  geom_edge_link(edge_width = 0.5, alpha = 0.2, color = "gray") +
  geom_node_point(aes(
    color = as.factor(V(g_est)$group),
    # shape = as.factor(V(g_est)$cis_trans),
    # fill = ifelse(V(g_est)$cis_trans == "cis", as.factor(V(g_est)$group), "white")
  ),shape = 15, size = 4) +
  scale_color_manual(values = community_colors, drop = FALSE)+
  theme_void() +                                 # Remove background, grid, axes
  theme(legend.position = "none")                # Remove all legendsp
p
ggsave(p, filename=paste(foldname, "/HighQuality_Network_",type_of_clustering, "_red_v2.png", sep=""), width=8, height=8, dpi=300)


