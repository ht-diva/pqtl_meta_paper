#===========================#
# Validate correct genotype data
# harmonization impact on COJO
#===========================#

path_meta <- "/scratch/dariush.ghasemi/projects/pqtl_conditional/results/meta_correct_sdy/"
plt_scatter <- "22-Nov-25_cojo_corrected_vs._previous_harmonized_genotype.png"

# read collected loci
loci_meta <- data.table::fread(glue(path_meta, "break/collected_loci_excluding_mhc.csv"))
cojo_meta <- data.table::fread(glue(path_meta, "cojo/collected_credible_sets.csv"))

# read cojo results using corrected harmonized genotype in Interval
cojo_harm <- data.table::fread(glue(path_meta, "../validate_geno/combined_conditional_snps.csv"))

# save loci in CHR 21 and 22 for testing COJO sensitivity 
# against re-harmonized Interval geneotype
seqids_chr2122 <- loci_meta %>%
  group_by(seqid) %>%
  mutate(n_loci = n()) %>%
  ungroup() %>%
  filter(chr %in% c(21:22)) %>%
  #head() %>%
  pull(seqid) #%>% paste(collapse = "|")

# grep -E "seqid.1|seqid.2|...|seqid.52" conf/path_meta_all.txt > conf/meta_loci_chr2122.txt


# sanity checks
all.equal(cojo_meta, cojo_harm) 
sapply(names(cojo_harm), function(col) sum(cojo_meta[[col]] != cojo_harm[[col]])) 

# Finds exact row differences
cojo_meta[cojo_harm, on = names(cojo_meta)]      # rows in meta not in harm
cojo_harm[cojo_meta, on = names(cojo_meta)]      # rows in harm not in meta

# avoids false mismatches due to ordering.
setkeyv(cojo_meta, names(cojo_meta))
setkeyv(cojo_harm, names(cojo_harm))

identical(cojo_meta, cojo_harm)



# merge two versions
cojo_join <- cojo_meta %>%
  #filter(seqid %in% seqids_chr2122) %>%
  select(seqid, locus, SNP, bC, mlog10pC) %>%
  full_join(
    cojo_harm %>% select(seqid, locus, SNP, bC, mlog10pC),
    join_by(seqid, locus, SNP), suffix = c("_previous", "_corrected")
  )

# scatter plot
cojo_join %>%
  pivot_longer(
    cols = c(bC_previous, bC_corrected, mlog10pC_previous, mlog10pC_corrected),
    names_to = c("index", ".value"),
    names_pattern = "(bC|mlog10pC)_(previous|corrected)"
  ) %>%
  ggplot(aes(x = previous, y = corrected)) +
  geom_abline(slope = 1, linetype = 2, color = "grey40")+
  geom_hline(yintercept = 0, color = "grey50") +
  geom_vline(xintercept = 0, color = "grey50")+
  geom_point(color = "steelblue", alpha = .6) +
  facet_wrap(~index, scales = "free")+
  labs(title = "12,095 conditional SNPs with/out corrected harmonized genotype") +
  theme_light() +
  theme(
    strip.background = element_blank(),
    strip.placement = "inside",
    strip.text.x = element_text(size = 10, color = "blue3", face = 2)
  )

ggsave(filename = plt_scatter, last_plot(), 
       width = 9.5, height = 5.5, dpi = 300, units = "in")

