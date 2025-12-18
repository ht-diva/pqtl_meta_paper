

path_hole_3m <- "/scratch/dariush.ghasemi/projects/pqtl_pipeline_finemap/conf/path_meta_all.txt"
path_bybatch <- "/scratch/dariush.ghasemi/projects/path_meta_bybatch.txt"

path_both <- "/scratch/dariush.ghasemi/projects/pqtl_pipeline_finemap/conf/path_merged.txt"
path_out_combined <- "/scratch/dariush.ghasemi/projects/pqtl_pipeline_finemap/results/cor_beta/combined.tsv"


meta_hole_3m <- data.table::fread(path_hole_3m, header = F, col.names = "path_3m") #%>% head(10)
meta_bybatch <- data.table::fread(path_bybatch, header = F, col.names = "path_bb") #%>% head(10)

# extract the protein sequence id
meta_hole_3m <- meta_hole_3m %>%
  dplyr::mutate(
    seqid = gsub(".gwaslab.tsv.bgz", "", basename(path_3m))
  )


meta_bybatch <- meta_bybatch %>%
  dplyr::mutate(seqid = gsub(".gwaslab.tsv.bgz", "", basename(path_bb)))


# merge two path files
meta_both <- meta_hole_3m %>% 
  left_join(meta_bybatch, by = join_by(seqid)) %>%
  relocate(seqid)


# store path file for defining wildcards in SMK
data.table::fwrite(meta_both, path_both, quote = F)




compute_cor <- function(path_a, path_b){
  
  pwas_a <- data.table::fread(path_a) %>% dplyr::select(SNPID, BETA)
  pwas_b <- data.table::fread(path_b) %>% dplyr::select(SNPID, BETA)
  
  seqid <- path_a %>% str_extract("seq.\\d+.\\d+")
  n_a <- nrow(pwas_a)
  n_b <- nrow(pwas_b)
  
  pwas_ab <- inner_join(
    pwas_a,
    pwas_b,
    join_by(SNPID),
    suffix = c("_a", "_b")
  )
  
  n_overlap <- nrow(pwas_ab)
  
  r <- cor(
    pwas_ab$BETA_a,
    pwas_ab$BETA_b,
    method = "pearson"
    )
  
  res <- data.frame(
    "seqid" = seqid,
    "n_a" = n_a,
    "n_b" = n_b,
    "n_overlap" = n_overlap,
    "r" = r, 
    "path_a" = path_a,
    "path_b" = path_b
    )
  
  return(res)
}



map2_dfr(
  meta_both$path_3m[1],
  meta_both$path_bb[1],
  compute_cor
)


pwas_a <- data.table::fread(meta_both$path_3m[1])
pwas_b <- data.table::fread(meta_both$path_bb[1])


cor(pwas_ab$BETA_a, pwas_ab$BETA_b, method = "pearson")



data.table::fread(path_out_combined) %>% View()
  #as_tibble() %>% dplyr(nsnps_b)
  select(r_betas) %>% summary()
  


data.table::fread(path_out_combined) %>%
  ggplot(aes(x = r_betas)) +
  #geom_violin()
  geom_histogram(bins = 100, fill = 'steelblue', color = "steelblue3")+
  labs(x = "Correlaion between BETA coefficients") +
  theme_bw()


ggsave("~/output_plot/29-Apr-25_histogram_of_pearson_correlation_between_betas.png",
       plot = last_plot(), height = 5.5, width = 7.5, dpi = 150)
