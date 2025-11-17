
library(tidyverse)
library(circlize)
library(ComplexHeatmap) # for adding legend

#BiocManager::install("ComplexHeatmap") # install and load the package

#-----------------#
# inputs
path_mapping  <- "/exchange/healthds/pQTL/Reference_datasets_for_QC_proteomics/Cis_trans_mapping/somascan_tss_ncbi_grch37_version_20241212.txt"
path_comunity <- "/scratch/dariush.ghasemi/projects/Hotspot_communities_connected.csv"
path_base <- "/group/diangelantonio/users/alessia_mapelli/pQTL/INTERVAL/Hotspots and lonespots/New_coloc_results/Coloc_networks_top_cond/"


#-----------------#
# prepare mapping file for phenotype/seqid position
map <- read.delim(path_mapping, sep=",")

map_4join <- map %>%
  dplyr::mutate(
    target = paste0("seq.", gsub("-", ".", SeqId))
  ) %>%
  dplyr::select(
    phenotype_id = target,
    phenotype_chr = chromosome,
    phenotype_pos = TSS,
    UniProt_ID
  )


#-----------------#
# list comminities
hotspots <- list.files(path_base, pattern = "^\\d+_\\[\\d+, \\d+\\]")

# Extract chromosome and start range using regex
chr   <- as.integer(sub("_.*", "", hotspots))
start <- as.integer(sub(".*\\[(\\d+),.*", "\\1", hotspots))

# Order based on chromosome and start
hotspots_sorted <- hotspots[order(chr, start)]

# test
test_commun <- data.table::fread(
  paste0(
    path_base, 
    hotspots_sorted[4], 
    "/Hotspot_communities_connected.csv")
  ) 


#-----------------#
# adding position of seqid TSS
append_mapping <- function(lb){
  
  lb %>%
    dplyr::select(
      phenotype_id,
      UniProt_ID,
      new_somamer, # legend key, groupping variable
      snp_chr = chr,
      snp_pos = POS,
      group
      ) %>%
    dplyr::mutate(group = str_replace_na(group, "Unlabled")) %>%
    #add_count(group) %>%
    #dplyr::filter(n > 4) %>%
    left_join(
      map_4join,
      join_by(phenotype_id, UniProt_ID), 
      relationship = "many-to-many"
      ) %>%
    dplyr::mutate(
      across(c(snp_chr, phenotype_chr), ~ str_c("chr", .)),
      new_old = ifelse(
        new_somamer == TRUE, 
        "Aptamers new in 7k assay version", 
        "Aptamers already assessed in previous assay versions")
      )
}

#df %>% dplyr::select(chr:SNPID, phenotype_id, group) %>% tail(10)

#-----------------#
# Main plot

# Set group-based colors (assign colors to groups)
sort_communities <- function(df){
  
  # First, Separate "Unlabled" (or "NA") cleanly
  group_levels <- unique(na.omit(df$new_old))  # -----------> change this and line 147
  group_levels_no_unlabeled <- setdiff(group_levels, "Unlabled")
  
  # Identify which are numeric-like
  suppressWarnings({
    numeric_flags <- !is.na(as.numeric(group_levels_no_unlabeled))
  })
  
  # Extract and sort numeric and non-numeric separately
  numeric_parts <- group_levels_no_unlabeled[numeric_flags]
  non_numeric_parts <- group_levels_no_unlabeled[!numeric_flags]
  
  # Sort numeric parts numerically, but preserve original string values
  numeric_parts_sorted <- numeric_parts[order(as.numeric(numeric_parts))]
  
  # Combine final result
  group_levels_sorted <- c(numeric_parts_sorted, sort(non_numeric_parts)) #, "Unlabled"
  
  return(group_levels_sorted)
}

#-----------------#
# Determine chromosome ranges (optional: can use real sizes or a fixed window)
setup_ranges <- function(df){
  
  df %>%
    dplyr::select(chr = snp_chr, pos = snp_pos) %>%
    bind_rows(
      df %>% dplyr::select(chr = phenotype_chr, pos = phenotype_pos)
    ) %>%
    group_by(chr) %>%
    summarise(min_pos = max(min(pos) - 1e6, 0), max_pos = max(pos) + 1e6) %>%
    ungroup() %>%
    dplyr::filter(!is.na(chr))
}

#-----------------#
draw_circos <- function(lb){
  
  # add protein's TSS position to loci
  comun_mapped <- append_mapping(lb)
  
  # sort communities
  group_levels_sorted <- sort_communities(comun_mapped)
  
  # set up colors
  colourCount = length(group_levels_sorted)
  
  #getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
  #group_palette <- setNames(RColorBrewer::brewer.pal(colourCount, "Set2"), group_levels_sorted)
  
  if (colourCount < 3) {
    base_palette <- hcl.colors(colourCount, palette = "Dynamic")
    base_palette <- c("#69b3a2", "#404080")
  } else {
    base_palette <- RColorBrewer::brewer.pal(min(colourCount, 8), "Set2")
    if (colourCount > 8) {
      base_palette <- colorRampPalette(base_palette)(colourCount)
    }
  }
  group_palette <- setNames(base_palette, group_levels_sorted)
  
  # grey for NA group
  group_palette[names(group_palette) == "Unlabled"] <- "#999999"
  
  # define community's color
  comun_mapped$color <- group_palette[as.character(comun_mapped$new_old)]
  
  # Gather all chromosomes needed
  all_chrs <- union(comun_mapped$snp_chr, comun_mapped$phenotype_chr)
  
  # define genomic ranges for each chromosome
  chr_ranges <- setup_ranges(comun_mapped)
  
  
  # Clear existing plot
  circos.clear()

  # draw the circuit
  circos.initializeWithIdeogram()

  all_sectors <- get.all.sector.index()
  
  # Annotate with UniProt ID at protein location (chr2)
  # label_df <- test_commun %>%
  #   append_mapping() %>%
  #   dplyr::filter(!is.na(new_somamer)) %>%
  #   dplyr::select(
  #     chr = phenotype_chr,
  #     pos = phenotype_pos,
  #     new_old
  #   ) %>%
  #   dplyr::filter(chr %in% get.all.sector.index())
  # 
  
  # circos.trackPlotRegion(
  #   factors = label_df$chr,
  #   y = label_df$pos,
  #   track.height = 0.05,
  #   bg.border = NA,
  #   panel.fun = function(x, y) {
  #     chr_current <- CELL_META$sector.index
  #     labels_chr <- label_df %>% dplyr::filter(chr == chr_current)
  #     
  #     for (i in seq_len(nrow(labels_chr))) {
  #       
  #       xpos <- labels_chr$pos[i]
  #       start_y <- CELL_META$ylim[2]
  #       end_y <- start_y - mm_y(6.5)
  #       mid_y <- (start_y - end_y) / 2
  #       
  #       # Connector from tile to mid
  #       circos.lines(
  #         x = c(xpos, xpos),
  #         y = c(start_y, mid_y),
  #         col = "gray40",
  #         lwd = 0.5
  #       )
  #       
  #       # Text label at midpoint
  #       circos.text(
  #         x = xpos,
  #         y = CELL_META$ylim[2] - mm_y(2.5),
  #         labels = labels_chr$new_old[i],
  #         facing = "clockwise",
  #         niceFacing = TRUE,
  #         adj = c(0, 0.5),
  #         cex = 0.4
  #       )
  #     }
  #   }
  # )
  
  # Create a legend object
  lgd <- ComplexHeatmap::Legend(
    labels = names(group_palette),
    title = "", # "Regional association"
    legend_gp = gpar(fill = group_palette),
    by_row = FALSE  # optional but silences warning
  )

  # Draw it outside the Circos plot
  ComplexHeatmap::draw(
    lgd,
    #x = unit(1, "npc") - unit(0, "mm"),
    x = unit(8, "cm"),
    y = unit(1, "mm"),
    just = c("bottom") #"right", 
  )
  
  # Add links: from SNPs to protein positions
  n_snps <- nrow(comun_mapped)

  for (i in seq_len(n_snps)) {
    
    # parse columns
    chr1 <- comun_mapped$snp_chr[i]
    chr2 <- comun_mapped$phenotype_chr[i]
    pos1 <- comun_mapped$snp_pos[i]
    pos2 <- comun_mapped$phenotype_pos[i]
    color  <- comun_mapped$color[i]
    label <- comun_mapped$new_somamer[i]
    
    # for debug
    if (!(chr2 %in% all_sectors)) {
      message("Skipping invalid chr2: ", chr2)
      next
    }
    
    # Sanity checks
    if (is.na(pos1) || is.na(pos2) || is.na(chr1) || is.na(chr2)) next
    if (!(chr1 %in% all_sectors) || !(chr2 %in% all_sectors)) next
    
    # Get sector data to validate positions
    chr2_limits <- get.cell.meta.data("xlim", sector.index = chr2)
    
    if (pos2 < chr2_limits[1] || pos2 > chr2_limits[2]) next
    
    # Links
    circos.link(
      sector.index1 = chr1,
      point1 = c(pos1, pos1),
      sector.index2 = chr2,
      point2 = c(pos2, pos2),
      col = color,
      lwd = 1
    )
    
  }
  
}

#-----------------#

# test
draw_circos(test_commun)


# Initialize circos plot
# Set parameters with minimal padding to avoid narrow-sector issue
# circos.par(
#   start.degree = 90, 
#   gap.degree = 2,
#   cell.padding = c(0.02, 0, 0.02, 0) # left, bottom, right, top
#   )
# 
# circos.initialize(
#   factors = chr_ranges$chr,
#   xlim = chr_ranges[, c("min_pos", "max_pos")]
#   )
# 
# # Add dummy track with chromosome labels
# circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
#   circos.text(CELL_META$xcenter, CELL_META$ylim[1] + mm_y(5),
#               CELL_META$sector.index, cex = 0.6, facing = "bending.inside", niceFacing = TRUE)
# }, bg.border = NA)

#--------------------------------#

lapply(
  4,
  function(j){
    hotspot <- hotspots_sorted[j]
    hotspot_clean <- str_remove_all(hotspot, "\\[|\\]| ") %>% str_replace(",", "_")
    plotname <- paste0("06-Aug-25_plt_circos_hotspot_chr", hotspot_clean, ".png")
  
    communities <- data.table::fread(paste0(path_base, hotspot, "/Hotspot_communities_connected.csv")) 
  
    png(plotname, width = 15.5, height = 15.5, units = "cm", res = 500)

    draw_circos(communities)

    dev.off()
  #return(communities %>% dim())
})

#png("plt_circos_hotspot_chr17_5_10.png", width = 15.5, height = 15.5, units = "cm", res = 500)
#pdf("plt_circos_hotspot_chr17_5_10.pdf", width = 7.5, height = 7.5)


#dev.off()

