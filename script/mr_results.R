rm(list=ls())
library(data.table)
library(dplyr)
library(tidyr)
library(Rmpfr)
library(stringr)
library(readxl)
library(ggplot2)

file1 <- read.delim("/group/diangelantonio/users/giulia.pontali/file_to_share/unconditional_IVs_from_MVP_to_be_merge_w_new_results_pvalNA.txt")

no_phecode <- "/exchange/healthds/pQTL/results/META_CHRIS_INTERVAL/MVP_MR/new_results_20250515/INTERVAL_CHRIS_MR_COJO_NONPHECODE_07MAY2025.txt.gz"
phecode <- "/exchange/healthds/pQTL/results/META_CHRIS_INTERVAL/MVP_MR/new_results_20250515/INTERVAL_CHRIS_MR_COJO_PHECODE_07MAY2025.txt.gz"
data_no_phecode <- fread(no_phecode, sep = "\t")
data_phecode <- fread(phecode, sep = "\t")

data <- rbind(data_no_phecode, data_phecode)
data <- as.data.frame(data)
new_ivs <- read.delim("/group/diangelantonio/users/giulia.pontali/file_to_share/unconditional_IVs_to_MVP_to_rerun.txt")
new_ivs$CHRPOS_ID <- paste0(new_ivs$CHR,"_",new_ivs$POS_38,"_",new_ivs$SeqID)
data <- merge(data, new_ivs[, c("CHRPOS_ID", "locus_START_END_37", "SNP", "POS_37")], by = "CHRPOS_ID", all.x = TRUE)
data$locus_START_END_37 <- gsub("-", "_", data$locus_START_END_37)
data$locus_START_END_37 <- paste0("chr", data$CHR, "_", data$locus_START_END_37)
class(data)


colnames(file1)[colnames(file1) == "UniProt_ID"] <- "UNIPROT"
file1 <- file1[,c(-2,-3)]
colnames(file1)[colnames(file1) == "symbol"] <- "HARMONIZED_GENE_NAME"
colnames(file1)
file1 <- file1[,c(-26)]
colnames(data)
data <- data[,c(-3,-6,-28)]

file1 <- file1[, sort(colnames(file1))]
data <- data[, sort(colnames(data))]

# Combine
final <- rbind(file1, data)

dim(final)
head(final)

final <- final[!grepl("Height", final$PHENOTYPE), ]
final$seqID <- sub(".*(seq\\.\\d+\\.\\d+).*", "\\1", final$CHRPOS_ID)
unique_triplets <- final %>% distinct(seqID, UNIPROT, PHENOTYPE) %>% nrow()
final$sign <- final$PVAL < (0.05/3497303)

table2 <- read_xlsx("/home/giulia.pontali/supplementary_table_2 (7).xlsx", sheet=2)
ann_uniprot <- table2[,c(5,37)]
ann_uniprot <- unique(ann_uniprot, by = "UniProt_ID")

ann_seqid <- table2[,c(1,33)]
ann_seqid <- unique(ann_seqid, by = "SeqID")

# Left-join UniProt annotation onto final
final2 <- merge(
  final,
  ann_uniprot,
  by.x = "UNIPROT",
  by.y = "UniProt_ID",
  all.x = TRUE,
  sort = FALSE
)

# Left-join SeqID annotation onto final2
final3 <- merge(
  final2,
  ann_seqid,
  by.x = "seqID",
  by.y = "SeqID",
  all.x = TRUE,
  sort = FALSE
)

final <- final3
table4 <- read_xlsx("/home/giulia.pontali/supplementary_table_4 (14).xlsx", sheet=2, skip=2)
table4 <- table4[,c("SeqID", "SNPID","tissues_with_encoding_gene_coloc", "encoding_gene_coloc_in_region")]

merged <- merge(
  x = final,
  y = table4,
  by.x = c("seqID", "SNP"),
  by.y = c("SeqID", "SNPID"),
  all.x = TRUE
)

final <- merged

setDT(final)
collapsed_final<- final[, lapply(.SD, function(x) paste(unique(x), collapse = "|")), 
                                    by = .(seqID, locus_START_END_37, BETA, PROTEIN_NAME, PHENOTYPE, CHR, PVAL, SE)]

write.table(collapsed_final, "/scratch/giulia.pontali/meta_file/collapsed_final.txt", sep="\t", quote=F, row.names = F)
