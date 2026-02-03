

path_pqtl_interval <- "/scratch/dariush.ghasemi/projects/pqtl_susie/results/interval/tmp/seq.24215.8_1_169258178_171220925_sumstat.csv"

headers <- c("CHR", "POS", "SNPID", "EA", "NEA", "EAF", "N", "BETA", "SE", "MLOG10P", "CHISQ")

# Read sumstat from INTERVAL GWAS
pqtl_interval <- fread(path_pqtl_interval, header = FALSE, col.names = headers, sep = "\t", data.table = FALSE)


#----------------------------------------#
#-----     Inspect GWAS Results     -----
#----------------------------------------#

# Check how many variants exist in eQTLs and 
# how many of them exist in Interval study genotype
length(pqtl_interval$POS)
length(prrx1_eqtl_hg37$pos_b37)

# No. common SNPs between eQTLs and pQTLs regions
intersect(prrx1_eqtl_hg37$pos_b37, pqtl_interval$POS) %>% length()

# Before harmonizing alleles, check alleles concordance 
pqtl_interval %>% filter(POS == "169836660") %>% select(SNPID, POS, EA, NEA)
prrx1_eqtl_hg37 %>% filter(pos_b37 == "169836660") %>% select(variant_id, pos_b37)


#----------------------------------------#
#-----    Sort GWAS Alphabetically  -----
#----------------------------------------#


# Harmonized sumstat with alphabetical order
prrx1_eqtl_harm <- prrx1_eqtl_hg37 %>%
  # define EA/NEA alleles based on variant_id from GTEx
  dplyr::mutate(
    CHR = str_split_fixed(variant_id, "_", 5)[,1] %>% str_remove("chr"),
    EA  = str_split_fixed(variant_id, "_", 5)[,4], # in GTEX, SNP ids are in chr_pos_ref_alt_build.
    NEA = str_split_fixed(variant_id, "_", 5)[,3],
    EAF = af,
    MLOG10P = -log10(pval_nominal)
  ) %>%
  # keep most significant one from multi-allelic variants
  group_by(pos_b37) %>% 
  slice_max(MLOG10P, n = 1) %>% 
  ungroup() %>% #distinct(pos_b37)
  # harmonization step
  mutate(
    swap = EA > NEA, # check if EA is not alphabetically preceding NEA 
    
    EA_new  = if_else(swap, NEA, EA),
    NEA_new = if_else(swap, EA, NEA),
    
    EAF_new    = if_else(swap, 1 - EAF, EAF),
    slope_new  = if_else(swap, -slope, slope),
    SNPID = str_c(CHR, pos_b37, EA_new, NEA_new, sep = ":")
  ) %>%
  # rename variables to match with susie pipe
  select(
    CHR,
    POS = pos_b37,
    SNPID,
    EA = EA_new,
    NEA = NEA_new,
    EAF = EAF_new,
    N = ma_samples, 
    BETA = slope_new,
    SE=slope_se,
    P=pval_nominal,
    MLOG10P
  )

# 8,292 SNPs left
# check lift-over
prrx1_eqtl_harm %>% filter(POS == "169836660")


#----------------------------------------#
#-----   Store SNPs to Subset PGEN  -----
#----------------------------------------#

# Having lifter over GTEx genomic positions to hg37, and sorting allels alphabetically,
# I will use SNPs list below to:
#   1. subset INTERVAL genotype in PGEN format
#   2. compute LD correlation matrix (Plink2: r-unphased) for susie fine-mapping

# Store list of SNP ids without multi-allelics
prrx1_eqtl_harm %>% select(SNPID) %>%
  write.table(
    "/scratch/dariush.ghasemi/projects/pqtl_susie/plink_ld/prrx1_eqtls_hg37_snps.list",
    quote = F, row.names = F, col.names = F)


#----------------------------------------#
#-----      GTEx eQTLs for SuSiE    -----
#----------------------------------------#

# Prepare sumstat for SuSiE fine-mapping
sumstat_gtex <- prrx1_eqtl %>%
  # GTEx variants are in CHR_POS_REF_ALT_b38 format
  dplyr::mutate(
    CHR = str_split_fixed(variant_id, "_", 5)[,1],
    POS = str_split_fixed(variant_id, "_", 5)[,2] %>% as.numeric(), # genomic position
    EA  = str_split_fixed(variant_id, "_", 5)[,4], # in GTEx, a2 or ALT is effect allele
    NEA = str_split_fixed(variant_id, "_", 5)[,3], # in GTEx, a1 or REF is non-effect allele
    CHR = str_remove(CHR, "chr"), # remove suffix for chromosome number
    MLOG10P = -log10(pval_nominal)
  ) %>%
  # add rsid
  inner_join(prrx1_snps, join_by(POS == pos)) %>% # filter(rsid == "rs7530405")
  # select most significant among multi-allelic SNPs 
  group_by(CHR, POS) %>%
  dplyr::slice_max(MLOG10P, n = 1) %>%
  ungroup() %>% #filter(rsid == "rs680084")
  # harmonization step
  mutate(
    swap = EA != ALT, # check if EA in GWAS does not match with ALT in GTEx genotype
    
    EA_new  = if_else(swap, NEA, EA),
    NEA_new = if_else(swap, EA, NEA),
    
    af_new    = if_else(swap, 1 - af, af),
    slope_new = if_else(swap, -slope, slope),
    SNPID = str_c(CHR, POS, EA_new, NEA_new, sep = ":")
  ) %>% #count(swap) #  TRUE: 856, FALSE: 5353
  # rename variables to match with susie pipe
  select(
    CHR,
    POS,
    SNPID, # = rsid,
    EA  = EA_new,
    NEA = NEA_new,
    EAF = af_new,
    N = ma_samples, 
    BETA = slope_new,
    SE = slope_se,
    P = pval_nominal,
    MLOG10P
    ) %>% # remove SNPs not matching with LD due to allele difference
  filter(!POS %in% c(170103075, 170128947))

