

library(data.table)
library(tidyverse)

# inputs
path_ld <- "/scratch/dariush.ghasemi/projects/ld_pipe/results/communities/ld/"
path_comun <- "/scratch/dariush.ghasemi/projects/ld_pipe/config/Hotspot_Community_Summary_Final_1.csv"

# outputs (save to Dariush home directory)
out_png   <- "18-Sep-25_ld_plot_communities_all_conditional_snps.png"
out_pdf_a <- "18-Sep-25_ld_plot_communities_all_cojo_snps_group_a.pdf"
out_pdf_b <- "18-Sep-25_ld_plot_communities_all_cojo_snps_group_b.pdf"


#---------------------------------------------#
#-----         LD plot PNG (test)      ------- 
#---------------------------------------------#

hotspot_files
test_hotspot <- "chr12_100_105"

ld_hotspot1 <- fread(glue(path_ld, test_hotspot, ".ld")) #"chr2_0_5.ld"

# read communities file
comunities <- fread(path_comun)


# ------   data wrangling   -----

# All SNPs
snps <- union(ld_hotspot1$SNP_A, ld_hotspot1$SNP_B)

# LD diagonal
ld_diagonal <- tibble(SNP_A = snps) %>%
  dplyr::mutate(
    CHR_A = as.integer(str_match(SNP_A, "(\\d+)")[,1]),
    BP_A  = as.integer(str_match(SNP_A, ":(\\d+):")[,2]),
    CHR_B = CHR_A,
    BP_B  = BP_A,
    SNP_B = SNP_A,
    R2 = 1
  ) %>%
  select(CHR_A, BP_A, SNP_A, ends_with("_B"), R2)

# missing LD segment
ld_missing <- ld_hotspot1 %>% 
  dplyr::select(
    CHR_A=CHR_B, BP_A=BP_B, SNP_A=SNP_B, CHR_B=CHR_A, BP_B=BP_A, SNP_B=SNP_A, R2
  )


# Extract SNPs and communities for the input hotspot
comunities_tidy <- comunities[hotspot2 == test_hotspot, ] %>%
  mutate(
    cojo_snp = str_remove_all(all_cond_SNP, " "),
    community = str_remove_all(community_id, "Community_") %>% as.integer()
    #C_A = as.integer(str_match(community_A, "(\\d+)")[,1]),
    #C_B = as.integer(str_match(community_B, "(\\d+)")[,1]),
  ) %>% 
  select(community, cojo_snp) %>%
  separate_rows(cojo_snp, sep = ",|;") %>% # separate communities with multiple conditional SNPs
  distinct()


# create a full LD matrix in long format
ld_long <- ld_hotspot1 %>%
  rbind(ld_missing) %>%
  rbind(ld_diagonal) %>%
  select(- CHR_A, - CHR_B) %>%
  mutate(rsquare = round(R2, 2) * 100) %>% # rescale LD
  left_join(comunities_tidy, join_by(SNP_A == cojo_snp)) %>%
  left_join(comunities_tidy, join_by(SNP_B == cojo_snp), suffix = c("_A", "_B"))


# Build SNP labels with community
snp_labels <- ld_long %>%
  distinct(SNP_A, BP_A, community_A) %>%
  mutate(label = paste0(SNP_A, " (", community_A, ")")) %>%
  arrange(community_A, BP_A)

# Extract order of labels
snp_order <- snp_labels$label



# -------  Heatmap -------

# draw LD plot with genom_tile
ld_long %>%
  mutate(
    snpcom_a = str_c(SNP_A, " (", community_A, ")"),
    snpcom_b = str_c(SNP_B, " (", community_B, ")"),
    label_a  = factor(snpcom_a, levels = snp_order),
    label_b  = factor(snpcom_b, levels = snp_order)
  ) %>%
  ggplot(aes(x = label_a, y = label_b, fill = rsquare)) +
  geom_tile(color = "black") +
  geom_text(aes(label = rsquare), size=3) +
  scale_fill_gradient(low = "white", high = "red") +
  #scale_fill_viridis_c(option = "plasma", limits = c(0,1)) +
  scale_x_discrete(position = "top") +
  labs(
    fill = expression(paste("LD (", r^2, ")")),
    title = paste0("Hotspot: ", test_hotspot)
    ) +
  # coord_fixed() +
  coord_equal() +
  #theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text  = element_text(size = 10, face  = 2),
    axis.text.x.top = element_text(angle = 45, hjust = 0),
    axis.ticks.length = unit(0.2, "cm"),
    #legend.title = ggtext::element_markdown(size = 10, face = 2, hjust = 1),
    legend.key.height = unit(2, "cm"),
    title = element_text(size = 12, face  = 2),
  )


# save in PNG
ggsave(filename = out_png,
       #height = 22, width = 23, dpi = 200
       height = 12, width = 12.5, dpi = 150
       )


#---------------------------------------------#
#-----      LD with Corrplot func.     ------- 
#---------------------------------------------#

# Make square matrix (wide form if you want)
ld_matrix <- ld_long %>%
  select(label_a, label_b, rsquare) %>%
  pivot_wider(names_from = label_b, values_from = rsquare, ) %>%
  column_to_rownames(var = "label_a")

library(corrplot)

corrplot.mixed(
  corr = ld_matrix %>% as.matrix(),
  is.corr = FALSE,
  upper.col = COL1('Reds', 5),
  lower.col = "black",
  tl.col = "black", # labels color
  tl.pos = "lt", # labels position
  #addCoef.col ='black',
  number.cex = 0.8,
  #diag = 'u',
  #type = 'full',
  lower = 'number',
  upper = 'square'
  )


#---------------------------------------------#
#-----          LD plot in PDF         ------- 
#---------------------------------------------#

# function to draw LD matrix plot 
create_ldplot <- function(ld_file){
  
  hotspot_name  <- gsub(".ld", "", ld_file)
  hotspot_label <- hotspot_dic[hotspot2 == hotspot_name, "hotspot"]
  
  full_path <- glue(path_ld, ld_file)
  
  # read LD from Plink 1.9
  ld_df <- fread(full_path)
  
  # All SNPs
  snps <- union(ld_df$SNP_A, ld_df$SNP_B)
  n_snps <- length(snps)
  
  # LD diagonal
  ld_diagonal <- tibble(SNP_A = snps) %>%
    dplyr::mutate(
      CHR_A = as.integer(str_match(SNP_A, "(\\d+)")[,1]),
      BP_A  = as.integer(str_match(SNP_A, ":(\\d+):")[,2]),
      CHR_B = CHR_A,
      BP_B  = BP_A,
      SNP_B = SNP_A,
      R2 = 1
    ) %>%
    select(CHR_A, BP_A, SNP_A, ends_with("_B"), R2)
  
  # missing LD segment
  ld_missing <- ld_df %>%
    dplyr::select(
      CHR_A = CHR_B,
      BP_A  = BP_B,
      SNP_A = SNP_B,
      CHR_B = CHR_A,
      BP_B  = BP_A,
      SNP_B = SNP_A,
      R2
    )
  
  # Extract SNPs and communities for the input hotspot
  comunities_tidy <- comunities[hotspot2 == hotspot_name, ] %>%
    mutate(
      cojo_snp  = str_remove_all(all_cond_SNP, " "),
      community = str_remove_all(community_id, "Community_") %>% as.integer()
    ) %>% 
    select(community, cojo_snp) %>%
    # separate communities with multiple conditional SNPs
    separate_rows(cojo_snp, sep = ",|;") %>%
    distinct()  # remove duplicated SNPs in a community

  
  # create a full LD matrix
  ld_long <- ld_df %>%
    rbind(ld_missing) %>%
    rbind(ld_diagonal) %>%
    select(- CHR_A, - CHR_B) %>%
    mutate(rsquare = round(R2, 2) * 100) %>% # rescale LD
    # add community number of each SNP to LD file
    left_join(comunities_tidy, join_by(SNP_A == cojo_snp)) %>%
    left_join(comunities_tidy, join_by(SNP_B == cojo_snp), suffix = c("_A", "_B"))
  
  
  # Build SNP labels with community
  snp_labels <- ld_long %>%
    distinct(SNP_A, BP_A, community_A) %>%
    mutate(label = paste0(SNP_A, " (", community_A, ")")) %>%
    arrange(community_A, BP_A)
  
  # number of distinct communities
  n_commun <- n_distinct(snp_labels$community_A)
  
  # Extract order of labels
  snp_order <- snp_labels$label
  
  # Heatmap
  ld_long %>%
    # create and sort pairs of SNP + community for x- and y-axis
    mutate(
      snpcom_a = str_c(SNP_A, " (", community_A, ")"),
      snpcom_b = str_c(SNP_B, " (", community_B, ")"),
      label_a  = factor(snpcom_a, levels = snp_order),
      label_b  = factor(snpcom_b, levels = snp_order)
    ) %>%
    ggplot(aes(x = label_a, y = label_b, fill = rsquare)) +
    geom_tile(color = "black") +
    geom_text(aes(label = rsquare), size=3) +
    scale_fill_gradient(low = "white", high = "red") +
    #scale_fill_viridis_c(option = "plasma", limits = c(0,1)) +
    scale_x_discrete(position = "top") +
    labs(
      fill = expression(paste("LD (", r^2, ")")),
      title = paste0("Hotspot: ", hotspot_label), # plot title
      subtitle = paste0(n_snps, " SNPs for ", n_commun, " community(ies)." )
    ) +
    coord_equal() +
    # coord_fixed() +
    # theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text  = element_text(size = 10, face  = 2),
      axis.text.x.top = element_text(angle = 45, hjust = 0),
      axis.ticks.length = unit(0.2, "cm"),
      legend.key.height = unit(2, "cm"),
      title = element_text(size = 12, face  = 2)
    )
  
}


#----------------------#
# create LD plots for 19 hotspots
pdf(
  out_pdf_b,
  #width = 13, height = 12 # group (a)
  width = 23, height = 22 # group (b)
  )

map(
  #hotspot_sorted[c(1,3,4,6,8:11,13:21)], # group (a)
  hotspot_sorted[c(2,5,7,12,22)],        # group (b)
  create_ldplot
  )

dev.off()


#----------------------#
# # Build complete SNP × SNP grid
# ld_long <- expand_grid(SNP_A = snps, SNP_B = snps) %>%
#   # Join your R2 values (both directions)
#   left_join(ld_hotspot1 %>% select(SNP_A, SNP_B, R2), 
#             by = c("SNP_A", "SNP_B")) %>%
#   left_join(ld_hotspot1 %>% 
#               transmute(SNP_A = SNP_B, SNP_B = SNP_A, R2_sym = R2),
#             by = c("SNP_A", "SNP_B")) %>%
#   mutate(R2 = coalesce(R2, R2_sym, NA)) %>%  # fill with sym value or 0 if missing
#   select(SNP_A, SNP_B, R2) %>%
#   # Force diagonal to 1
#   mutate(R2 = if_else(SNP_A == SNP_B, 1, R2)) %>%
#   distinct(SNP_A, SNP_B, .keep_all = TRUE)

#----------------------#
# check no. of COJO SNPs in each hotspot
map_dbl(
  #hotspot_sorted[c(1,3,4,6,8:11,13:21)],
  hotspot_sorted[c(2,5,7,12,22)],
  .f = function(ld_file){
    ld_df <- fread(glue(path_ld, ld_file))
    snps <- union(ld_df$SNP_A, ld_df$SNP_B)
    return(length(snps))
    }
  )


