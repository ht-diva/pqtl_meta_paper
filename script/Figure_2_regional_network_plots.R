setwd("/group/diangelantonio/users/alessia_mapelli/pQTL/INTERVAL/Hotspots and lonespots/New_coloc_results/Coloc_networks_top_cond")
rm(list=ls(all=TRUE))

if (!requireNamespace("ggraph", quietly = TRUE)) {
  install.packages("ggraph")
}

library(readr)
library(readxl)
library(dplyr)
library(igraph)
library(ggplot2)
library(scales)
library(ggrepel)
library(Matrix)
library(ggraph)
library(tidygraph)
library(data.table)

colocalization_results_filtered_path ="/group/diangelantonio/users/alessia_mapelli/pQTL/INTERVAL/Hotspots and lonespots/New_coloc_results/Colocalization_results_filtered.csv"
hotspots_df_path = "/group/diangelantonio/users/alessia_mapelli/pQTL/INTERVAL/Hotspots and lonespots/New_coloc_results/Hotspots_df.csv"
type_of_clustering = "connected" 

colocalization_results_filtered <- fread(colocalization_results_filtered_path)
hotspot_df <- read_csv(hotspots_df_path)[,-1]

LB <- read_excel("/group/diangelantonio/users/alessia_mapelli/pQTL/INTERVAL/Hotspots and lonespots/New_coloc_results/supplementary_table_3.xlsx", sheet = 2)
colnames(LB) <- LB[3,]
LB <- LB[-c(1,2,3),]

#################################################
## Figure 2a
#################################################
LB_in_hotsposts <- LB %>%
  filter(hotspot== TRUE)

LB_in_hotsposts$CHR <- as.character(as.numeric(LB_in_hotsposts$CHR))

sel_hotspot <- "[25000000, 30000000]"
chr_hotspot <- "2"

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
comp <- components(g_est)
cat("Community computation done \n")

#Annotate results
V(g_est)$group <- comp$membership
communities <-data.frame(nodes = V(g_est)$name,
                         group = comp$membership)

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
foldname <- paste(chr_hotspot,"_", sel_hotspot,"/plots", sep ="")
dir.create(foldname)

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

communities_adam_plot <- communities %>%
  select(SNPID, bp_cum_map, group, cis_or_trans, HARMONIZED_GENE_NAME) %>%
  mutate(POS = as.numeric(stringr::word(SNPID, 2, sep = ":"))) %>%
  mutate(
    group = ifelse(is.na(group), "Not colocalized signals", as.character(group)),
    group = factor(group, levels = c("1", "Not colocalized signals"))
  )
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

mb_format <- function(x) sprintf("%.2f", x / 1e6)

label_data <- subset(
  communities_adam_plot,
  cis_or_trans == "cis" & !is.na(POS) & !is.na(bp_cum_map)
)

community_names <- communities_adam_plot %>%
  count(group, SNPID, name = "n") %>%        # count how often each SNPID appears per group
  group_by(group) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>% # keep SNPID with highest count in each group
  ungroup()

p <- ggplot(
  communities_adam_plot,
  aes(x = POS, y = bp_cum_map, color = group, shape = cis_or_trans)
) +
  geom_point(size = 4.5) +   # increased point size
  
  scale_color_manual(
    values = c("1" = gg_color_hue(1), 
               "Not colocalized signals" = "gray"),
    name = "Communities",
    labels = c("1" = "2:27730940:C:T", 
               "Not colocalized signals" = "Not colocalized signals")
  ) +
  
  # change square (15) to triangle (17)
  scale_shape_manual(
    values = c("cis" = 16,   # circle
               "trans" = 17), # triangle
    name = "Signal Type"
  ) +
  
  scale_x_continuous(
    labels = mb_format,
    limits = c(24000000, 32000000)
  ) +
  
  scale_y_continuous(
    breaks = df3$center,
    labels = c(1:22)
  ) +
  
  labs(
    x = "\npQTL positions (Mb)",
    y = "\nCoding gene position (Chr)\n"
  ) +
  
  theme_bw() +
  
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 13),
    
    # move legend inside bottom-left corner
    #legend.position = c(0.02, 0.02),
    #legend.justification = c(0, 0),
    legend.position = c(0.98, 0.5),
    legend.justification = c(1, 0.5),
    
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key.size = grid::unit(1, "cm"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.spacing.x = grid::unit(1, "cm")
  )
p
plot_name <- paste0(foldname, "/Regional_plot_", chr_hotspot,"_", sel_hotspot,".svg")
ggsave(plot_name, plot = p, width = 30, height = 23, units = "cm")


num_colors <- length(unique(V(g_est)$group))
palette <- scales::hue_pal()(num_colors)

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
ggsave(p, 
       filename=paste0(foldname, "/Network_plot_", chr_hotspot,"_", sel_hotspot,".svg"),
       width=8,
       height=8,
       dpi=300)

p_net <- ggraph(layout) +
  geom_edge_link(edge_width = 0.5, alpha = 0.2, color = "gray") +
  
  geom_node_point(
    aes(
      color = factor(group)
    ),
    fill="white",
    shape = 21,
    size = 3.5,
    stroke = 1.1
  ) +
  
  scale_color_manual(
    name = "Communities",
    values = palette,
    labels = c(
      "1" = "2:27730940:C:T"
    )
  ) +
  
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )

p_net
ggsave(
  filename = paste0(foldname, "/Network_plot_", chr_hotspot, "_", sel_hotspot, "_nodes_names.svg"),
  plot = p_net,
  width = 10,
  height = 10,
  dpi = 300
)

lead_snp_map <- c(
  "1" = "2:27730940:C:T"
)

node_meta <- communities %>%
  transmute(
    node_name = nodes,
    group_meta = as.character(group),
    cis_or_trans = cis_or_trans,
    HARMONIZED_GENE_NAME = HARMONIZED_GENE_NAME
  ) %>%
  distinct(node_name, .keep_all = TRUE) %>%
  mutate(
    snp_id = sub("^[^:]+:", "", node_name),
    highlight_fill = group_meta %in% names(lead_snp_map) & snp_id == lead_snp_map[group_meta],
    label_node = group_meta %in% c("1") & cis_or_trans == "cis" & !is.na(HARMONIZED_GENE_NAME)
  )

tg <- as_tbl_graph(g_est) %>%
  activate(nodes) %>%
  left_join(node_meta, by = c("name" = "node_name")) %>%
  mutate(
    group = ifelse(is.na(group), "Other", group),
    filled_group = ifelse(highlight_fill, group, NA_character_)
  )
community_cols <- c(
  "1" = palette
)
layout <- create_layout(tg, layout = "graphopt")

p_net <- ggraph(layout) +
  geom_edge_link(edge_width = 0.5, alpha = 0.2, color = "gray") +
  
  geom_node_point(
    aes(
      color = factor(group),
      fill = filled_group
    ),
    shape = 21,
    size = 3.5,
    stroke = 1.1
  ) +
  
  geom_node_text(
    data = layout %>% filter(label_node),
    aes(label = HARMONIZED_GENE_NAME),
    repel = TRUE,
    size = 3.5,
    point.padding = unit(0.2, "lines"),
    box.padding = unit(0.35, "lines"),
    max.overlaps = Inf
  ) +
  
  scale_color_manual(
    name = "Communities",
    values = community_cols,
    labels = c(
      "1" = "2:27730940:C:T",
      "Other" = "Other nodes"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "1" = gg_color_hue(1)
    ),
    na.value = "white",
    guide = "none"
  ) +
  
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )

p_net
ggsave(
  filename = paste0(foldname, "/Network_plot_", chr_hotspot, "_", sel_hotspot, "_nodes_names_fill_by_top_cond_snp.svg"),
  plot = p_net,
  width = 10,
  height = 10,
  dpi = 300
)

#################################################
## Figure 2c
#################################################
LB_in_hotsposts <- LB %>%
  filter(hotspot== TRUE)

LB_in_hotsposts$CHR <- as.character(as.numeric(LB_in_hotsposts$CHR))

sel_hotspot <- "[185000000, 190000000]"
chr_hotspot <- "3"

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
comp <- components(g_est)
cat("Community computation done \n")

#Annotate results
V(g_est)$group <- comp$membership
communities <-data.frame(nodes = V(g_est)$name,
                         group = comp$membership)

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
foldname <- paste(chr_hotspot,"_", sel_hotspot,"/plots", sep ="")
dir.create(foldname)

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
  select(nodes, SNPID, POS_37, bp_cum_map, group, cis_or_trans, HARMONIZED_GENE_NAME) %>%
  rename(POS = POS_37)
communities_adam_plot$POS <- as.numeric(communities_adam_plot$POS)
communities_adam_plot <- communities_adam_plot[communities_adam_plot$POS < 187000000, ]
communities_adam_plot <- communities_adam_plot[communities_adam_plot$POS > 186000000, ]
str(communities_adam_plot)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


communities_adam_plot <- communities_adam_plot %>%
  mutate(
    group = as.character(group),
    group = ifelse(is.na(group), "Not connected", group),
    group = factor(group)
  ) %>%
  rename(Communities=group)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

community_names <- communities_adam_plot %>%
  count(Communities, SNPID, name = "n") %>%        # count how often each SNPID appears per group
  group_by(Communities) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>% # keep SNPID with highest count in each group
  ungroup()

communities_res <- communities_adam_plot %>%
  mutate(
    Communities_recode = case_when(
      Communities== "6" ~ "3:186459227:A:G",  
      Communities== "1" ~ "3:186393786:A:G",  
      Communities== "12" ~ "3:186445052:G:T",  
      Communities == "Not connected"  ~ "Not colocalized signals",
      TRUE ~ "Others"                                   
    ),
    Communities_recode = factor(
      Communities_recode,
      levels = c("3:186459227:A:G","3:186393786:A:G","3:186445052:G:T","Others","Not colocalized signals")
    )
  )


community_colors <- c(
  "3:186459227:A:G" = gg_color_hue(3)[1],           # blue
  "3:186393786:A:G" = gg_color_hue(3)[2],           # green
  "3:186445052:G:T" = gg_color_hue(3)[3],           # red
  "Others" = "darkgrey",
  "Not colocalized signals" = "lightgrey"
)


mb_format <- function(x) sprintf("%.2f", x / 1e6)

labels_df <- communities_res %>%
  filter(Communities %in% c("6", "1", "12")) %>%
  filter(
    cis_or_trans == "cis",
    !is.na(start),
    !is.na(bp_cum_map),
    !is.na(HARMONIZED_GENE_NAME)
  ) %>%
  arrange(bp_cum_map) %>%                     # "first signal" = smallest bp_cum_map
  distinct(Communities,HARMONIZED_GENE_NAME, .keep_all = TRUE)

p <- ggplot(
  communities_res,
  aes(x = POS, y = bp_cum_map,
      color = Communities_recode,
      shape = cis_or_trans)
) +
  
  geom_point(size = 4.5) +   # increased point size
  
  geom_text_repel(
    data = labels_df,
    aes(label = HARMONIZED_GENE_NAME),
    size = 4,
    max.overlaps = 70,
    nudge_x = 0.07 * diff(range(communities_res$POS, na.rm = TRUE)),
    direction = "x",
    vjust = 0
  ) +
  
  scale_color_manual(
    values = community_colors,
    name = "Communities"
  ) +
  
  # square -> triangle
  scale_shape_manual(
    values = c("cis" = 16,      # circle
               "trans" = 17),   # triangle
    name = "Signal Type"
  ) +
  
  scale_x_continuous(labels = mb_format) +
  
  scale_y_continuous(
    breaks = df3$center,
    labels = c(1:22)
  ) +
  
  labs(
    x = "\npQTL positions (Mb)",
    y = "\nCoding gene position (Chr)\n"
  ) +
  
  theme_bw() +
  
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 13),
    
    # legend on right
    legend.position = c(0.98, 0.5),
    legend.justification = c(1, 0.5),
    legend.background = element_rect(fill = "white", color = "black"),
    
    legend.key.size = grid::unit(1, "cm"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    
    legend.box = "vertical"
  )
p
ggsave(
  paste0(foldname, "/Regional_plot_", chr_hotspot,"_", sel_hotspot,".svg"),
  plot = p,
  width = 30,
  height = 23,
  units = "cm"
)


keep_groups <- c(1, 6, 12)
nodes_keep  <- V(g_est)[group %in% keep_groups]

cols <- gg_color_hue(3)

community_colors <- setNames(
  rep("darkgray", 18),
  as.character(1:18)
)
community_colors[c("6", "1", "12")] <- cols[1:3]


tg <- as_tbl_graph(g_est)
layout <- create_layout(tg, layout = 'graphopt',charge = 0.0001) 
p <- ggraph(layout) +
  geom_edge_link(edge_width = 0.5, alpha = 0.2, color = "gray") +
  geom_node_point(aes(
    color = as.factor(V(g_est)$group),
    # shape = as.factor(V(g_est)$cis_trans),
    # fill = ifelse(V(g_est)$cis_trans == "cis", as.factor(V(g_est)$group), "white")
  ),shape = 15, size = 2) +
  scale_color_manual(values = community_colors, drop = FALSE)+
  theme_void() +                                 # Remove background, grid, axes
  theme(legend.position = "none")                # Remove all legendsp
p
ggsave(p, 
       filename=paste0(foldname, "/Network_plot_", chr_hotspot,"_", sel_hotspot,".svg"),
       width=8,
       height=8,
       dpi=300)



# communities of interest
highlight_map <- c(
  "6"  = "3:186459227:A:G",
  "1"  = "3:186393786:A:G",
  "12" = "3:186445052:G:T"
)

# Build node annotation table from the regional-plot data
node_meta <- communities_res %>%
  transmute(
    node_name = nodes,
    group_meta = as.character(Communities),
    cis_or_trans = cis_or_trans,
    HARMONIZED_GENE_NAME = HARMONIZED_GENE_NAME
  ) %>%
  distinct(node_name, .keep_all = TRUE) %>%
  mutate(
    snp_id = sub("^[^:]+:", "", node_name),
    highlight_fill = case_when(
      group_meta %in% names(highlight_map) & snp_id == highlight_map[group_meta] ~ TRUE,
      TRUE ~ FALSE
    ),
    label_node = group_meta %in% c("1", "6", "12") & cis_or_trans == "cis" & !is.na(HARMONIZED_GENE_NAME)
  )

tg <- as_tbl_graph(g_est) %>%
  activate(nodes) %>%
  left_join(node_meta, by = c("name" = "node_name")) %>%
  mutate(
    group = as.character(group),
    group = ifelse(is.na(group), "Other", group),
    
    # outline colors by community
    plot_color = case_when(
      group == "6"  ~ gg_color_hue(3)[1],
      group == "1"  ~ gg_color_hue(3)[2],
      group == "12" ~ gg_color_hue(3)[3],
      TRUE ~ "darkgray"
    ),
    
    # fill only the representative SNP of each highlighted community
    plot_fill = case_when(
      highlight_fill & group == "6"  ~ gg_color_hue(3)[1],
      highlight_fill & group == "1"  ~ gg_color_hue(3)[2],
      highlight_fill & group == "12" ~ gg_color_hue(3)[3],
      TRUE ~ "white"
    )
  )

layout <- create_layout(tg, layout = "graphopt", charge = 0.0001)
p_net <- ggraph(layout) +
  geom_edge_link(edge_width = 0.5, alpha = 0.2, color = "gray") +
  
  geom_node_point(
    aes(
      color = factor(group),
      fill = ifelse(highlight_fill, as.character(group), NA)
    ),
    shape = 21,
    size = 3.5,
    stroke = 1.1
  ) +
  
  geom_node_text(
    data = layout %>% filter(label_node),
    aes(label = HARMONIZED_GENE_NAME),
    repel = TRUE,
    size = 3.5,
    point.padding = unit(0.2, "lines"),
    box.padding = unit(0.35, "lines"),
    max.overlaps = Inf
  ) +
  
  scale_color_manual(
    name = "Communities",
    values = c(
      "6" = gg_color_hue(3)[1],
      "1" = gg_color_hue(3)[2],
      "12" = gg_color_hue(3)[3],
      "Other" = "darkgray"
    ),
    labels = c(
      "6" = "3:186459227:A:G",
      "1" = "3:186393786:A:G",
      "12" = "3:186445052:G:T",
      "Other" = "Other nodes"
    )
  ) +
  
  # 🔹 Fill (no legend)
  scale_fill_manual(
    values = c(
      "6" = gg_color_hue(3)[1],
      "1" = gg_color_hue(3)[2],
      "12" = gg_color_hue(3)[3]
    ),
    na.value = "white",
    guide = "none"  
  ) +
  
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )

p_net
ggsave(
  filename = paste0(foldname, "/Network_plot_", chr_hotspot, "_", sel_hotspot, "_nodes_names_fill_by_top_cond_snp.svg"),
  plot = p_net,
  width = 10,
  height = 10,
  dpi = 300
)

p_net <- ggraph(layout) +
  geom_edge_link(edge_width = 0.5, alpha = 0.2, color = "gray") +
  
  geom_node_point(
    aes(
      color = factor(group)
    ),
    fill="white",
    shape = 21,
    size = 3.5,
    stroke = 1.1
  ) +
  
  geom_node_text(
    data = layout %>% filter(label_node),
    aes(label = HARMONIZED_GENE_NAME),
    repel = TRUE,
    size = 3.5,
    point.padding = unit(0.2, "lines"),
    box.padding = unit(0.35, "lines"),
    max.overlaps = Inf
  ) +
  
  scale_color_manual(
    name = "Communities",
    values = c(
      "6" = gg_color_hue(3)[1],
      "1" = gg_color_hue(3)[2],
      "12" = gg_color_hue(3)[3],
      "Other" = "darkgray"
    ),
    labels = c(
      "6" = "3:186459227:A:G",
      "1" = "3:186393786:A:G",
      "12" = "3:186445052:G:T",
      "Other" = "Other nodes"
    )
  ) +
  
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )

p_net
ggsave(
  filename = paste0(foldname, "/Network_plot_", chr_hotspot, "_", sel_hotspot, "_nodes_names.svg"),
  plot = p_net,
  width = 10,
  height = 10,
  dpi = 300
)

