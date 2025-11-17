
library(data.table)
library(tidyverse)

#----------------------#
# inputs
path_ld <- "/scratch/dariush.ghasemi/projects/ld_pipe/results/communities/ld/"
path_comun <- "/scratch/dariush.ghasemi/projects/ld_pipe/config/Hotspot_Community_Summary_Final_1.csv"

# outputs
plt_pdf <- "12-Sep-25_ld_plots_for_hotspots.pdf"

#----------------------#
# list of LD files
hotspot_files <- list.files(path = path_ld, pattern = ".ld$")

# sort hotspots for PDF
hotspot_sorted <- tibble(file = hotspot_files) %>%
  mutate(
    chr   = as.integer(str_match(file, "^chr(\\d+)_")[,2]),
    start = as.integer(str_match(file, "^chr\\d+_(\\d+)_")[,2]),
    end   = as.integer(str_match(file, "^chr\\d+_\\d+_(\\d+)\\.ld$")[,2])
  ) %>%
  arrange(chr, start, end) %>%
  pull(file)

#----------------------#
# read community file for building dictionary
comunities <- fread(path_comun)

# build dic for hotspot names
hotspot_dic <- comunities %>%
  distinct(hotspot, .keep_all = T) %>%
  mutate(chrom = str_remove_all(hotspot2, "_.*$")) %>%
  arrange(nchar(chrom)) %>%
  dplyr::select(hotspot, hotspot2)


#----------------------#
# discrete colors for LD plot
ld_colors <- c(
  "1"   = "red",
  "0.8" = "orange",
  "0.6" = "darkgreen",
  "0.4" = "skyblue",
  "0.2" = "blue"
  )


# function to draw LD matrix plot 
draw_ldplot <- function(ld_file){
  
  file_name <- gsub(".ld", "", ld_file)
  hotspot_name <- hotspot_dic[hotspot2 == file_name, "hotspot"]
  
  full_path = glue(path_ld, ld_file)
  
  ld_df <- fread(full_path)
  
  # All SNPs
  snps <- union(ld_df$SNP_A, ld_df$SNP_B)
  
  # Extract positions
  snp_order <- tibble(SNP = snps) %>%
    mutate(pos = as.integer(str_match(SNP, ":(\\d+):")[,2])) %>%
    arrange(pos) %>%
    pull(SNP)
  
  # Build complete SNP × SNP grid
  ld_long <- expand_grid(SNP_A = snps, SNP_B = snps) %>%
    # Join your R2 values (both directions)
    left_join(ld_df %>% select(SNP_A, SNP_B, R2), 
              by = c("SNP_A", "SNP_B")) %>%
    left_join(ld_df %>% 
                transmute(SNP_A = SNP_B, SNP_B = SNP_A, R2_sym = R2),
              by = c("SNP_A", "SNP_B")) %>%
    mutate(R2 = coalesce(R2, R2_sym, NA)) %>%  # fill with sym value or 0 if missing
    select(SNP_A, SNP_B, R2) %>%
    # Force diagonal to 1
    mutate(R2 = if_else(SNP_A == SNP_B, 1, R2)) %>%
    distinct(SNP_A, SNP_B, .keep_all = TRUE) %>%
    mutate(rsquare = cut(R2, 
                         breaks=c(-Inf, .2, .4, .6, .8, Inf), 
                         labels=c("0.2", "0.4", "0.6", "0.8", "1")),
           r2 = round(R2, 2))
  
  # create plot
  ld_long %>%
    ggplot(aes(x = factor(SNP_A, levels = snp_order),
               y = factor(SNP_B, levels = snp_order),
               fill = rsquare)) +
    geom_tile(color = "black") +
    geom_text(aes(label = r2)) +
    #scale_fill_gradient(low = "white", high = "red") +
    scale_fill_manual(
      values = ld_colors, na.value = "white",
      breaks = c("1", "0.8", "0.6", "0.4", "0.2"),
      labels = c("1", "0.8", "0.6", "0.4", "0.2"),
      ) +
    scale_x_discrete(position = "top") +
    coord_fixed() +
    labs(
      fill = expression(paste("LD (", r^2, ")")),
      title = paste0("Hotspot: ", hotspot_name)
      ) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text.x.top = element_text(angle = 45, hjust = 0)
      )
}

#----------------------#
# create LD plots for 19 hotspots
pdf(plt_pdf, width = 9.5, height = 9.5)
map(hotspot_sorted, draw_ldplot)
dev.off()

