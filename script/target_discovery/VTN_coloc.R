
# A slight complication to the story, which we might be able to ignore,
# is that there are 2 seqIDs for VTN, which both have rs704 as the top variant. 
# But one seqID (seq.8280.238, which only has a weak cis-pQTL with MLOG10P=11)
# is the cis-pQTL in the hugely pleiotropic group 6, whereas the other 
# seqID (seq.13125.45, which has a strong cis-pQTL with MLOG10P>2000) 
# is the cis-pQTL in the other smaller groups 1,2,3 & 4). Looks like 
# we didn’t test coloc for seq.13125.45 with group 6 as rs704 doesn’t 
# end up in the one variant credible set for some reason, despite being 
# the index variant -> Dariush can you check? maybe we can force to run coloc?
  
path_rds_base <- "/scratch/dariush.ghasemi/projects/pqtl_pipeline_finemap/results/meta_correct_sdy/cojo/"
path_rds_loc1 <- "seq.13125.45/finemaping/{snpid}_locus_chr1_196494547_196984679_finemap.rds"
path_rds_loc2 <- "seq.8280.238/finemaping/{snpid}_locus_chr1_193837851_199237004_finemap.rds"

# check RDS contents
glue(
  path_rds_base,
  "seq.8280.238/finemaping/17:26694861:A:G_locus_chr17_26592946_26739127_finemap.rds"
  #"seq.13125.45/finemaping/17:26605158:G:T_locus_chr17_25283026_29975846_finemap.rds"
  ) %>% 
  readRDS() %>% head()



prepare4coloc <- function(data){
  temp  <- data %>% dplyr::rename(beta=bC, pvalues=pC)
  odata <- as.list(na.omit(temp))
  odata$type <- unique(odata$type)
  odata$sdY <- unique(odata$sdY)
  
  return(odata)
}



prepare_combinations <- function(trait_A, locus_A, trait_B, locus_B){
  
  path_A <- str_c(path_rds_base, trait_A, "/finemaping")
  path_B <- str_c(path_rds_base, trait_B, "/finemaping")
  
  name_A <- paste0('*._locus_chr', locus_A, "*._finemap.rds")
  name_B <- paste0('*._locus_chr', locus_B, "*._finemap.rds")
  
  files_A <- list.files(path = path_A, pattern = name_A, full.names = T)
  files_B <- list.files(path = path_B, pattern = name_B, full.names = T)
  
  comb_AB <- expand.grid(files_A, files_B, stringsAsFactors = FALSE)
  
  return(comb_AB)
}


iterate_coloc <- function(RDS_A, RDS_B){
  
  snp_A <- RDS_A %>% str_remove(".*finemaping/") %>% str_remove("_locus_chr.*")
  snp_B <- RDS_B %>% str_remove(".*finemaping/") %>% str_remove("_locus_chr.*")
  
  seq_A <- RDS_A %>% str_remove("/finemaping.*") %>% basename()
  seq_B <- RDS_B %>% str_remove("/finemaping.*") %>% basename()
  
  
  # run coloc standard
  res <- coloc::coloc.abf(
    dataset1 = readRDS(RDS_A) |> prepare4coloc(),
    dataset2 = readRDS(RDS_B) |> prepare4coloc()
  )
  
  res_h4 <- res$summary %>% t() %>% as.data.frame() %>% select(nsnps, PP.H4.abf)
  res_final <- data.frame(
    "seqid1" = seq_A,
    "seqid2" = seq_B,
    "SNP1" = snp_A, 
    "SNP2" = snp_B
  ) %>% 
    cbind(res_h4)
  
  return(res_final)
}


comb_all <- prepare_combinations(
  "seq.8280.238", "17_26592946_26739127",
  "seq.13125.45", "17_25283026_29975846"
)


res_combin <- map2_df(comb_all$Var1, comb_all$Var2, iterate_coloc)

# show results to the group
res_combin %>%
  #filter(PP.H4.abf > 0.7) %>%
  DT::datatable(
    #caption = 'Table 1: This colocalization results for two pQTLs at VTN locus.'
    caption = htmltools::tags$caption(
      style = 'caption-side: top; text-align: left;',
      'Table 1: ', htmltools::em('This colocalization results for two pQTLs at VTN locus.')
    )
  ) %>% 
  DT::formatSignif(columns = c('PP.H4.abf'), digits = 3) #formatRound('PP.H4.abf', 3)



