# FN3K locus plots with local INTERVAL LD and LocusZoom-like gene annotation
# Version: fixed rs28485881 LD reference, hg37 INTERVAL pgen, GRCh37-to-GRCh38 mapping,
#          MR-variant preflight check, GENCODE gene track, LocusZoom-style options
#
# This version uses ONE fixed LD reference SNP across all panels:
#   rs28485881 = chr17:80690014:A:G (GRCh37) = chr17:82732138:G:A (GRCh38)
#
# LD is computed from a local INTERVAL PLINK2 pgen/pvar/psam subset in GRCh37.
# Plot coordinates can be GRCh37 or GRCh38. No UCSC liftOver/chain file is used.
# In this FN3K interval, the supplied POS/POS_hg38 map is used to estimate and
# validate a constant regional coordinate offset, which is then applied to all
# INTERVAL and GWAS variants. Variant matching prioritises allele-canonicalised
# CHR:POS:A1:A2 keys and falls back to unique physical positions.
#
# Important safeguards:
#   * The script stops if any of the three MR variants are absent from the INTERVAL panel.
#   * The script stops on duplicated INTERVAL variant IDs.
#   * Allele order is canonicalised, so A:G and G:A are treated as the same variant.
#   * LD output is joined by unique INTERVAL variant ID and then mapped to GRCh38;
#     it is never joined by a duplicated/missing ID or by position alone.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(coloc)
  library(ggplot2)
  library(patchwork)
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
  library(grid)
})

# ---------- User settings ----------
setwd("/Users/c.giambartolomei/Downloads/FN3K/")

a1c_file <- "A1C_Mean_INT_FN3K_locus.txt"
t2d_file <- "phe_250_2_FN3K_locus.txt"
magic_file <- "/Users/c.giambartolomei/Downloads/FN3K/MAGIC/Chen1000G/MAGIC1000G_FN3K_region/MAGIC1000G_2hGlu_EUR_FN3K_chr17_80462551_81047829_b37.tsv"
magic_trait_label <- "2-hour glucose"

cond_dir <- "/Users/c.giambartolomei/Downloads/FN3K/FN3K_conditional_gtex/"

cond_files <- c(
  file.path(cond_dir, "FN3K_seq.21147.9_17_80791139_C_T_Whole_Blood_for_coloc.txt"),
  file.path(cond_dir, "FN3K_seq.21147.9_17_80752787_G_T_Whole_Blood_for_coloc.txt"),
  file.path(cond_dir, "FN3K_seq.21147.9_17_80690014_A_G_Whole_Blood_for_coloc.txt")
)


# Prefix of the LOCAL COPY of the INTERVAL FN3K-region PLINK2 files.
# Expected files:
#   /Users/c.giambartolomei/Downloads/FN3K/FN3K_INTERVAL_LD_subset/HDS_BELIEVE_FN3K_hg37.pgen
#   /Users/c.giambartolomei/Downloads/FN3K/FN3K_INTERVAL_LD_subset/HDS_BELIEVE_FN3K_hg37.pvar
#   /Users/c.giambartolomei/Downloads/FN3K/FN3K_INTERVAL_LD_subset/HDS_BELIEVE_FN3K_hg37.psam
local_ld_pfile_prefix <- "/Users/c.giambartolomei/Downloads/FN3K/FN3K_INTERVAL_LD_subset/HDS_BELIEVE_FN3K_hg37"

# FN3K pQTL regional file containing GRCh37 POS and GRCh38 POS_hg38.
# This is used only as a build-conversion/variant-identity map for LD colouring.
build_map_file <- "/Users/c.giambartolomei/Downloads/FN3K/seq.21147.9_FN3K_chr17_80462551_81047829_hg38.tsv"

# Local GENCODE GRCh38 gene annotation file for the LocusZoom-like gene track.
# Download once, for example:
#   cd /Users/c.giambartolomei/Downloads/FN3K/
#   curl -L -o gencode.v50.annotation.gtf.gz \
#     https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz
gene_gtf_file_grch38 <- "/Users/c.giambartolomei/Downloads/FN3K/gencode.v50.annotation.gtf.gz"
# For GRCh37 plotting, use the GENCODE v19 hg19/GRCh37 annotation:
# curl -L -o /Users/c.giambartolomei/Downloads/FN3K/gencode.v19.annotation.gtf.gz \
#   https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
gene_gtf_file_grch37 <- "/Users/c.giambartolomei/Downloads/FN3K/gencode.v19.annotation.gtf.gz"

# If TRUE, download the required build-specific GENCODE GTF automatically when missing.
auto_download_gene_gtf <- TRUE
gene_gtf_url_grch37 <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
gene_gtf_url_grch38 <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz"

# No chain files are required. The script estimates the local GRCh37-to-GRCh38
# offset from the supplied pQTL build map and verifies that it is constant across
# the locus before applying it to variants absent from that map.

# Plot coordinate system. Set to 37 for GRCh37/hg19 or 38 for GRCh38.
plot_build <- 37L
if (!plot_build %in% c(37L, 38L)) stop("plot_build must be 37 or 38")
gene_gtf_file <- if (plot_build == 37L) gene_gtf_file_grch37 else gene_gtf_file_grch38
plot_build_label <- if (plot_build == 37L) "GRCh37" else "GRCh38"

# ---------- Plot options ----------
# "all" = produce one page per conditional pQTL file.
# "ld_reference_only" = produce only the conditional pQTL page corresponding to
# the fixed LD reference variant rs28485881 / 17:80690014:A:G b37.
# conditional_plot_mode <- "all"
conditional_plot_mode ="ld_reference_only"

# Draw vertical lines for the other two MR variants. Set to FALSE for a cleaner
# LocusZoom-like figure focused only on the LD reference MR variant.
show_additional_mr_variants <- FALSE

# Show gene symbols in the gene track. LocusZoom-style panels are often cleaner
# with this set to FALSE, especially in dense regions.
show_gene_labels <- TRUE

# Restrict the gene track to protein-coding genes only. Set TRUE for a simpler
# figure; FALSE shows the full GENCODE annotation in the region.
gene_track_protein_coding_only <- FALSE

# Keep the output visually close to LocusZoom: compact legends, minimal grid,
# transcript/exon gene track, and one fixed LD reference across all panels.

# Full FN3K interval represented in the INTERVAL subset (GRCh37 coordinates).
# These bounds match the supplied pQTL file and include all three MR instruments.
ld_panel_check_region_start37 <- 80462551L
ld_panel_check_region_end37   <- 81047829L

plink2_bin <- "plink2"

# Separate cache from the previous 1000G analysis so stale 1000G LD can never be reused.
cache_dir <- file.path(getwd(), "FN3K_INTERVAL_LD_cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# Increase this if the INTERVAL files, build map, or LD method change.
cache_version <- paste0("INTERVAL_unphased_regional_offset_plot_b", plot_build, "_v6")

# MR instruments. These are the LD lead variants for each conditional file.
instrument_table <- data.frame(
  instrument = c(
    "LD reference: rs28485881",
    "Other MR variant: 17:80752787:G:T",
    "Other MR variant: 17:80791139:C:T"
  ),
  snp37 = c("17:80690014:A:G", "17:80752787:G:T", "17:80791139:C:T"),
  pos37 = c(80690014, 80752787, 80791139),
  snp38 = c("17:82732138:G:A", "17:82794911", "17:82833263"),
  pos = c(82732138, 82794911, 82833263),
  rsid = c("rs28485881", NA, NA),
  color = c("purple", "grey45", "grey45"),
  line_type = c("solid", "dashed", "dashed"),
  stringsAsFactors = FALSE
)
instrument_table$pos38 <- instrument_table$pos
instrument_table$pos <- if (plot_build == 37L) instrument_table$pos37 else instrument_table$pos38

# Fixed LD reference for every plot.
# This is the MR instrument with evidence: instrument 1.
ld_reference <- list(
  chr = 17,
  pos37 = 80690014,
  pos38 = 82732138,
  rsid = "rs28485881",
  snp37 = "17:80690014:A:G",
  label = "rs28485881 / chr17:82732138 / 17:80690014:A:G b37"
)
ld_reference$plot_pos <- if (plot_build == 37L) ld_reference$pos37 else ld_reference$pos38

get_plot_instruments <- function(show_additional = show_additional_mr_variants) {
  if (isTRUE(show_additional)) {
    return(instrument_table)
  }

  instrument_table %>%
    dplyr::filter(snp37 == ld_reference$snp37)
}

select_conditional_files <- function(files, mode = conditional_plot_mode) {
  mode <- match.arg(mode, c("all", "ld_reference_only"))

  if (mode == "all") {
    return(files)
  }

  pattern <- gsub(":", "_", ld_reference$snp37)
  keep <- grepl(pattern, basename(files), fixed = TRUE)

  if (sum(keep) != 1) {
    stop(
      "conditional_plot_mode = 'ld_reference_only' expected exactly one file matching ",
      pattern, ", but found ", sum(keep), ".\nFiles were:\n",
      paste(basename(files), collapse = "\n")
    )
  }

  files[keep]
}


format_conditioned_label <- function(file_label) {
  # Use a manuscript-friendly title instead of the full filename.
  # Examples: 17_80690014_A_G -> 17:80690014:A:G.
  m <- regexec("17_([0-9]+)_([ACGT])_([ACGT])", file_label)
  hit <- regmatches(file_label, m)[[1]]
  if (length(hit) == 4) {
    return(paste0("Conditioned on 17:", hit[2], ":", hit[3], ":", hit[4], " (GRCh37)"))
  }
  file_label
}

# ---------- General helpers ----------
get_col <- function(df, candidates, label = "") {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) {
    stop("None of these columns found for ", label, ": ",
         paste(candidates, collapse = ", "),
         "\nAvailable columns: ", paste(names(df), collapse = ", "))
  }
  df[[hit[1]]]
}

extract_pos38 <- function(df) {
  if ("pos38" %in% names(df)) return(as.numeric(df$pos38))
  if ("POS_hg38" %in% names(df)) return(as.numeric(df$POS_hg38))
  if ("variant_id_38" %in% names(df)) {
    return(as.numeric(sapply(strsplit(as.character(df$variant_id_38), "_", fixed = TRUE), `[`, 2)))
  }
  if ("bp" %in% names(df)) return(as.numeric(df$bp))
  if ("POS" %in% names(df)) return(as.numeric(df$POS))
  stop("No position column found.")
}

diagnose_duplicates <- function(df, snp_col = "snp", label = "dataset", remove = TRUE) {
  dup <- duplicated(df[[snp_col]]) | duplicated(df[[snp_col]], fromLast = TRUE)
  dup_df <- df[dup, ]

  message("\n", label)
  message("  total rows = ", nrow(df))
  message("  duplicated rows = ", nrow(dup_df))
  message("  duplicated SNP IDs = ", length(unique(dup_df[[snp_col]])))

  if (nrow(dup_df) > 0) {
    diff_report <- dup_df %>%
      dplyr::group_by(.data[[snp_col]]) %>%
      dplyr::summarise(
        differing_columns = paste(
          names(.)[sapply(., function(x) dplyr::n_distinct(x, na.rm = FALSE) > 1)],
          collapse = ", "
        ),
        .groups = "drop"
      ) %>%
      dplyr::filter(differing_columns != "")

    message("  duplicated SNP IDs with differing columns = ", nrow(diff_report))
    if (nrow(diff_report) > 0) print(diff_report)
  }

  if (remove) df[!dup, ] else df
}

get_pp4_label <- function(pairs, pair_name) {
  if (!pair_name %in% names(pairs)) return("PP4 = NA")
  pp4 <- pairs[[pair_name]]$result$summary["PP.H4.abf"]
  paste0("PP4 = ", signif(pp4, 3))
}

safe_mlog10 <- function(p) {
  p <- as.numeric(p)
  p[p <= 0] <- .Machine$double.xmin
  -log10(p)
}

resolve_lead_from_instrument <- function(index_snp, reference_df, instruments = instrument_table) {
  hit <- instruments %>% dplyr::filter(snp37 == index_snp)

  if (nrow(hit) == 1) {
    return(list(
      pos38 = as.numeric(hit$pos[1]),
      rsid = hit$rsid[1],
      label = ifelse(!is.na(hit$rsid[1]) && hit$rsid[1] != "",
                     paste0(hit$rsid[1], " / chr17:", hit$pos[1]),
                     paste0("chr17:", hit$pos[1])),
      instrument = hit$instrument[1]
    ))
  }

  warning("Could not map index_snp = ", index_snp, "; using top p-value SNP instead.")
  top_pos <- reference_df %>%
    dplyr::filter(!is.na(pos), !is.na(p)) %>%
    dplyr::arrange(p) %>%
    dplyr::slice(1) %>%
    dplyr::pull(pos) %>%
    as.numeric()

  list(
    pos38 = top_pos,
    rsid = NA_character_,
    label = paste0("chr17:", top_pos),
    instrument = NA_character_
  )
}

# ---------- Local INTERVAL PLINK2 LD helpers ----------
check_local_pgen <- function(pfile_prefix) {
  needed <- paste0(pfile_prefix, c(".pgen", ".pvar", ".psam"))
  missing <- needed[!file.exists(needed)]
  if (length(missing) > 0) {
    stop(
      "Missing local INTERVAL PLINK2 files:\n",
      paste(missing, collapse = "\n"),
      "\n\nCopy the GRCh37 FN3K subset pgen/pvar/psam to the stated prefix first."
    )
  }
}

read_pvar_table <- function(pfile_prefix) {
  pvar_file <- paste0(pfile_prefix, ".pvar")
  pvar <- data.table::fread(pvar_file, data.table = FALSE, skip = "#CHROM")
  names(pvar)[1] <- sub("^#", "", names(pvar)[1])

  required <- c("CHROM", "POS", "ID", "REF", "ALT")
  missing <- setdiff(required, names(pvar))
  if (length(missing) > 0) {
    stop("PVAR file is missing required columns: ", paste(missing, collapse = ", "))
  }

  pvar <- pvar %>%
    dplyr::mutate(
      CHROM = gsub("^chr", "", as.character(CHROM)),
      POS = as.numeric(POS),
      ID = as.character(ID),
      REF = as.character(REF),
      ALT = as.character(ALT)
    )

  if (any(is.na(pvar$ID) | pvar$ID == "" | pvar$ID == ".")) {
    stop("INTERVAL PVAR contains missing/'.' variant IDs; unique IDs are required for safe LD parsing.")
  }

  dup_ids <- unique(pvar$ID[duplicated(pvar$ID)])
  if (length(dup_ids) > 0) {
    stop(
      "INTERVAL PVAR contains ", length(dup_ids), " duplicated variant ID(s).\n",
      "Examples: ", paste(head(dup_ids, 10), collapse = ", "),
      "\nRecreate the subset with duplicate checking before running this script."
    )
  }

  pvar
}

read_psam_n <- function(pfile_prefix) {
  nrow(data.table::fread(paste0(pfile_prefix, ".psam"), data.table = FALSE))
}

canonical_variant_key <- function(chr, pos, allele1, allele2) {
  chr <- gsub("^chr", "", as.character(chr))
  pos <- as.numeric(pos)
  allele1 <- toupper(as.character(allele1))
  allele2 <- toupper(as.character(allele2))
  lo <- ifelse(allele1 <= allele2, allele1, allele2)
  hi <- ifelse(allele1 <= allele2, allele2, allele1)
  paste(chr, pos, lo, hi, sep = ":")
}

parse_colon_variant <- function(x, label = "variant") {
  x <- as.character(x)
  x <- sub("^chr", "", x)
  spl <- strsplit(x, ":", fixed = TRUE)
  bad <- lengths(spl) < 4
  if (any(bad)) {
    stop("Could not parse ", label, " as CHR:POS:A1:A2. Examples: ",
         paste(head(x[bad], 5), collapse = ", "))
  }
  data.frame(
    chr = vapply(spl, `[`, character(1), 1),
    pos = as.numeric(vapply(spl, `[`, character(1), 2)),
    a1 = vapply(spl, `[`, character(1), 3),
    a2 = vapply(spl, `[`, character(1), 4),
    stringsAsFactors = FALSE
  )
}


read_build_map <- function(file = build_map_file) {
  if (!file.exists(file)) stop("Build-map pQTL file not found: ", file)
  x <- data.table::fread(file, data.table = FALSE)

  required <- c("SNPID", "POS", "POS_hg38")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0) {
    stop("Build-map file is missing columns: ", paste(missing, collapse = ", "))
  }

  parsed <- parse_colon_variant(x$SNPID, "SNPID in build_map_file")
  map <- data.frame(
    snpid37 = as.character(x$SNPID),
    chr37 = gsub("^chr", "", parsed$chr),
    pos37 = as.numeric(x$POS),
    a1 = parsed$a1,
    a2 = parsed$a2,
    pos38 = as.numeric(x$POS_hg38),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(key37 = canonical_variant_key(chr37, pos37, a1, a2)) %>%
    dplyr::filter(!is.na(pos37), !is.na(pos38), !is.na(key37))

  # Exact duplicated SNPIDs are allowed only when they carry the same mapping.
  conflicts <- map %>%
    dplyr::group_by(snpid37) %>%
    dplyr::summarise(
      n_key = dplyr::n_distinct(key37),
      n_pos38 = dplyr::n_distinct(pos38),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_key > 1 | n_pos38 > 1)
  if (nrow(conflicts) > 0) {
    stop("Conflicting duplicated SNPIDs in build_map_file. First examples:\n",
         paste(utils::capture.output(print(head(conflicts, 10))), collapse = "\n"))
  }

  map <- map %>% dplyr::distinct(key37, .keep_all = TRUE)

  # A canonical biallelic key must map to one GRCh38 position.
  key_conflicts <- map %>%
    dplyr::group_by(key37) %>%
    dplyr::summarise(n_pos38 = dplyr::n_distinct(pos38), .groups = "drop") %>%
    dplyr::filter(n_pos38 > 1)
  if (nrow(key_conflicts) > 0) {
    stop("Some canonical GRCh37 variant keys map to multiple GRCh38 positions.")
  }

  message("Build map loaded: ", nrow(map), " unique allele-canonicalised variants.")
  map
}

estimate_regional_build_offset <- function(build_map) {
  offsets <- as.numeric(build_map$pos38) - as.numeric(build_map$pos37)
  offsets <- offsets[is.finite(offsets)]
  if (length(offsets) == 0) stop("No valid POS/POS_hg38 pairs were available to estimate the regional build offset.")

  offset_tab <- sort(table(offsets), decreasing = TRUE)
  modal_offset <- as.numeric(names(offset_tab)[1])
  agreement <- mean(offsets == modal_offset)

  message(
    "Regional GRCh37->GRCh38 offset: +", modal_offset, " bp; agreement in build map: ",
    sum(offsets == modal_offset), "/", length(offsets), " (", sprintf("%.2f", 100 * agreement), "%)."
  )

  # This shortcut is safe only when essentially the entire small locus has one offset.
  if (agreement < 0.995) {
    bad <- build_map[offsets != modal_offset, c("snpid37", "pos37", "pos38")]
    stop(
      "The FN3K build map does not have a sufficiently constant regional offset (agreement <99.5%).
",
      "Do not use the offset conversion. First examples of discordant mappings:
",
      paste(utils::capture.output(print(head(bad, 10))), collapse = "
")
    )
  }
  modal_offset
}

annotate_pvar_with_build38 <- function(pvar, build_map, regional_offset) {
  pvar2 <- pvar %>%
    dplyr::mutate(
      pos37 = as.numeric(POS),
      is_multiallelic = grepl(",", ALT, fixed = TRUE),
      key37 = ifelse(is_multiallelic, NA_character_, canonical_variant_key(CHROM, POS, REF, ALT))
    ) %>%
    dplyr::left_join(
      build_map %>% dplyr::select(key37, pos38_map = pos38, snpid37),
      by = "key37", relationship = "many-to-one"
    ) %>%
    dplyr::mutate(
      pos38_offset = pos37 + regional_offset,
      pos38 = dplyr::coalesce(as.numeric(pos38_map), as.numeric(pos38_offset)),
      map_source = dplyr::case_when(
        !is.na(pos38_map) ~ "allele_map",
        !is.na(pos38_offset) ~ "regional_offset",
        TRUE ~ "unmapped"
      )
    )

  validation <- pvar2 %>%
    dplyr::filter(!is.na(pos38_map)) %>%
    dplyr::summarise(n = dplyr::n(), mismatches = sum(pos38_map != pos38_offset, na.rm = TRUE))
  if (validation$n > 0 && validation$mismatches > 0) {
    stop("Regional offset disagrees with ", validation$mismatches,
         " allele-aware mappings in the INTERVAL PVAR.")
  }

  message(
    "INTERVAL PVAR variants: ", nrow(pvar2),
    "; multi-ALT retained: ", sum(pvar2$is_multiallelic, na.rm = TRUE),
    "; allele-map validated: ", sum(pvar2$map_source == "allele_map"),
    "; regional-offset mapped: ", sum(pvar2$map_source == "regional_offset"),
    "; unmapped: ", sum(pvar2$map_source == "unmapped")
  )
  pvar2
}

find_variant_id_by_snp37 <- function(pvar_annotated, snp37) {
  q <- parse_colon_variant(snp37, "MR instrument")
  qchr <- gsub("^chr", "", as.character(q$chr[1]))
  qpos <- as.numeric(q$pos[1])
  qkey <- canonical_variant_key(qchr, qpos, q$a1[1], q$a2[1])

  # Preferred match: chromosome, position and allele-order-independent allele pair.
  hit <- pvar_annotated %>% dplyr::filter(key37 == qkey)
  if (nrow(hit) == 1) return(hit$ID[1])
  if (nrow(hit) > 1) {
    stop("More than one INTERVAL variant has the same canonical key as ", snp37,
         ". Matching IDs: ", paste(hit$ID, collapse = ", "))
  }

  # Some INTERVAL files encode the tested/effect allele pair differently from the
  # pQTL SNPID, while the physical variant is still uniquely identified by position.
  # Fall back only when there is exactly one PVAR row at that chromosome/position.
  by_pos <- pvar_annotated %>%
    dplyr::filter(CHROM == qchr, POS == qpos)

  if (nrow(by_pos) == 1) {
    warning(
      "No exact allele-pair match for ", snp37,
      "; using the unique INTERVAL variant at chr", qchr, ":", qpos,
      " (ID=", by_pos$ID[1], ", REF=", by_pos$REF[1], ", ALT=", by_pos$ALT[1], ")."
    )
    return(by_pos$ID[1])
  }

  # Last check: the ID itself may contain chr:pos:ref:alt even if REF/ALT metadata
  # was represented differently. Canonicalise parseable IDs and compare.
  parseable <- grepl("^(chr)?[0-9XYM]+:[0-9]+:[^:]+:[^:]+$", pvar_annotated$ID)
  if (any(parseable)) {
    idp <- parse_colon_variant(pvar_annotated$ID[parseable], "INTERVAL PVAR ID")
    idkeys <- canonical_variant_key(idp$chr, idp$pos, idp$a1, idp$a2)
    idhit <- pvar_annotated[which(parseable)[idkeys == qkey], , drop = FALSE]
    if (nrow(idhit) == 1) return(idhit$ID[1])
  }

  stop(
    "Could not identify ", snp37, " in the INTERVAL PVAR. ",
    "Exact canonical allele match: 0; rows at chr", qchr, ":", qpos, ": ", nrow(by_pos), ".\n",
    if (nrow(by_pos) > 0) paste(utils::capture.output(print(by_pos[, c("CHROM","POS","ID","REF","ALT")])), collapse = "\n") else "",
    "\nPlease inspect with: grep -P '^17\\t", qpos, "\\t' <file>.pvar"
  )
}

get_hg37_bounds_from_hg38 <- function(build_map, region_start38, region_end38) {
  m <- build_map %>%
    dplyr::filter(pos38 >= region_start38, pos38 <= region_end38)
  if (nrow(m) == 0) {
    stop("No build-map variants overlap GRCh38 region ", region_start38, "-", region_end38)
  }
  c(start = floor(min(m$pos37, na.rm = TRUE)), end = ceiling(max(m$pos37, na.rm = TRUE)))
}

check_mr_variants_in_ld_panel <- function(pfile_prefix = local_ld_pfile_prefix,
                                          build_map = read_build_map(),
                                          regional_offset = estimate_regional_build_offset(build_map)) {
  message("Preflight: checking the three MR variants in the INTERVAL GRCh37 LD panel...")
  check_local_pgen(pfile_prefix)
  pvar <- read_pvar_table(pfile_prefix)
  pvar_annotated <- annotate_pvar_with_build38(pvar, build_map, regional_offset)

  checked <- instrument_table %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      pvar_id = find_variant_id_by_snp37(pvar_annotated, snp37),
      mapped_pos38 = {
        z <- pvar_annotated %>% dplyr::filter(ID == pvar_id)
        if (nrow(z) == 1) z$pos38[1] else NA_real_
      },
      expected_pos38 = as.numeric(pos38),
      map_matches_expected = !is.na(mapped_pos38) && mapped_pos38 == expected_pos38
    ) %>%
    dplyr::ungroup()

  if (any(!checked$map_matches_expected)) {
    stop(
      "Preflight failed: one or more MR variants do not map to the expected GRCh38 position.\n",
      paste(utils::capture.output(print(checked)), collapse = "\n")
    )
  }

  message("Preflight passed: all three MR variants are present and map correctly.")
  print(checked %>% dplyr::select(instrument, snp37, pvar_id, mapped_pos38))
  invisible(list(pvar = pvar_annotated, build_map = build_map, checked = checked))
}

parse_plink_ld_file <- function(ld_file, pvar_annotated, lead_id, plot_build = plot_build) {
  ld <- data.table::fread(ld_file, data.table = FALSE)
  names(ld)[1] <- sub("^#", "", names(ld)[1])

  required <- c("ID_A", "ID_B")
  missing <- setdiff(required, names(ld))
  if (length(missing) > 0) {
    stop("PLINK LD output is missing columns: ", paste(missing, collapse = ", "),
         ". Columns were: ", paste(names(ld), collapse = ", "))
  }
  r2_col <- intersect(c("UNPHASED_R2", "R2", "PHASED_R2"), names(ld))[1]
  if (is.na(r2_col)) stop("Could not identify the r2 column in PLINK output.")

  a <- ld %>%
    dplyr::filter(ID_A == lead_id) %>%
    dplyr::transmute(ID = as.character(ID_B), ld_r2 = as.numeric(.data[[r2_col]]))
  b <- ld %>%
    dplyr::filter(ID_B == lead_id) %>%
    dplyr::transmute(ID = as.character(ID_A), ld_r2 = as.numeric(.data[[r2_col]]))
  out <- dplyr::bind_rows(a, b) %>% dplyr::filter(!is.na(ld_r2))

  if (nrow(out) == 0) {
    warning("No PLINK LD rows involving lead ID ", lead_id, ".")
    return(data.frame(pos37=numeric(0), pos38=numeric(0), pos=numeric(0), ld_r2=numeric(0)))
  }

  out %>%
    dplyr::left_join(
      pvar_annotated %>% dplyr::select(ID, pos37 = POS, pos38, key37),
      by = "ID", relationship = "many-to-one"
    ) %>%
    dplyr::mutate(pos = if (plot_build == 37L) as.numeric(pos37) else as.numeric(pos38)) %>%
    dplyr::filter(!is.na(pos)) %>%
    dplyr::select(key37, pos37, pos38, pos, ld_r2)
}

get_ld_for_region_local_plink <- function(chr,
                                          lead_snp37,
                                          region_start,
                                          region_end,
                                          pfile_prefix = local_ld_pfile_prefix,
                                          build_map = read_build_map(),
                                          regional_offset = estimate_regional_build_offset(build_map),
                                          plot_build = plot_build,
                                          plink2 = plink2_bin,
                                          cache_dir = cache_dir) {
  chr <- gsub("^chr", "", as.character(chr))
  check_local_pgen(pfile_prefix)

  pvar <- read_pvar_table(pfile_prefix)
  pvar_annotated <- annotate_pvar_with_build38(pvar, build_map, regional_offset)
  lead_id <- find_variant_id_by_snp37(pvar_annotated, lead_snp37)
  lead_plot_pos <- if (plot_build == 37L) ld_reference$pos37 else ld_reference$pos38

  # The INTERVAL subset is already the desired GRCh37 locus. Use its actual bounds
  # so indels/multiallelic variants are retained and no second subsetting is needed.
  bounds37 <- range(pvar$POS[pvar$CHROM == chr], na.rm = TRUE)

  cache_file <- file.path(
    cache_dir,
    paste0("LD_UNPHASED_", cache_version, "_chr", chr,
           "_lead37_", gsub(":", "_", lead_snp37),
           "_plot_", region_start, "_", region_end, ".tsv")
  )
  if (file.exists(cache_file)) {
    message("Reading cached INTERVAL LD: ", cache_file)
    return(data.table::fread(cache_file, data.table = FALSE))
  }

  message("Computing INTERVAL unphased genotype LD for ", lead_snp37,
          " (PVAR ID ", lead_id, "); GRCh37 bounds chr", chr, ":",
          bounds37[1], "-", bounds37[2], "; N samples = ", read_psam_n(pfile_prefix))

  out_prefix <- file.path(cache_dir,
    paste0("plink2_INTERVAL_ld_chr", chr, "_lead", ld_reference$pos37, "_", Sys.getpid()))
  args <- c(
    "--pfile", pfile_prefix,
    "--chr", chr,
    "--from-bp", as.integer(bounds37[1]),
    "--to-bp", as.integer(bounds37[2]),
    "--ld-snp", lead_id,
    "--r2-unphased", "allow-ambiguous-allele",
    "--ld-window-kb", "1000",
    "--ld-window", "999999",
    "--ld-window-r2", "0",
    "--out", out_prefix
  )
  status <- system2(plink2, args = args)
  if (!identical(status, 0L)) stop("PLINK2 failed while calculating INTERVAL LD. Check: ", out_prefix, ".log")

  ld_file <- paste0(out_prefix, ".vcor")
  ld_file_zst <- paste0(ld_file, ".zst")
  if (!file.exists(ld_file) && file.exists(ld_file_zst)) {
    zstd_status <- system2("zstd", args = c("-d", "-f", ld_file_zst))
    if (!identical(zstd_status, 0L)) stop("zstd decompression failed: ", ld_file_zst)
  }
  if (!file.exists(ld_file)) stop("PLINK2 did not create ", ld_file)

  ld_out <- parse_plink_ld_file(ld_file, pvar_annotated, lead_id, plot_build) %>%
    dplyr::filter(pos >= region_start, pos <= region_end) %>%
    dplyr::bind_rows(data.frame(
      key37 = canonical_variant_key(17, ld_reference$pos37, "A", "G"),
      pos37 = ld_reference$pos37, pos38 = ld_reference$pos38,
      pos = lead_plot_pos, ld_r2 = 1
    )) %>%
    dplyr::group_by(key37, pos37, pos38, pos) %>%
    dplyr::summarise(
      ld_r2 = max(ld_r2, na.rm = TRUE), .groups = "drop"
    )

  message("Mapped INTERVAL LD variants in plotted ", plot_build_label, " region: ", nrow(ld_out))
  print(table(cut(ld_out$ld_r2, breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf))))
  data.table::fwrite(ld_out, cache_file, sep = "\t")
  ld_out
}

add_ld_to_plot_df <- function(df, ld_tbl, lead_pos) {
  df <- df %>% dplyr::mutate(pos = as.numeric(pos), pos37 = as.numeric(pos37))

  # First, use allele-order-independent keys where available. This distinguishes
  # different variants at the same position and tolerates A:G versus G:A.
  ld_by_key <- ld_tbl %>%
    dplyr::filter(!is.na(key37), key37 != "") %>%
    dplyr::select(key37, ld_r2_key = ld_r2) %>%
    dplyr::distinct(key37, .keep_all = TRUE)
  df <- df %>% dplyr::left_join(ld_by_key, by = "key37", relationship = "many-to-one")

  # Position fallback is used only when the INTERVAL LD table has exactly one
  # variant at that GRCh37 position, so no allele is guessed.
  ld_by_pos <- ld_tbl %>%
    dplyr::filter(!is.na(pos37)) %>%
    dplyr::group_by(pos37) %>%
    dplyr::summarise(
      n_variants_at_pos = dplyr::n_distinct(key37, na.rm = TRUE),
      ld_r2_pos = if (dplyr::n_distinct(key37, na.rm = TRUE) <= 1) max(ld_r2, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )
  df <- df %>%
    dplyr::left_join(ld_by_pos, by = "pos37", relationship = "many-to-one") %>%
    dplyr::mutate(
      ld_r2 = dplyr::coalesce(ld_r2_key, ld_r2_pos),
      ld_match_method = dplyr::case_when(
        !is.na(ld_r2_key) ~ "allele_key",
        is.na(ld_r2_key) & !is.na(ld_r2_pos) ~ "unique_position",
        TRUE ~ "unmatched"
      ),
      ld_r2 = ifelse(pos == lead_pos, 1, ld_r2),
      ld_bin = dplyr::case_when(
        is.na(ld_r2) ~ "NA",
        ld_r2 >= 0.8 ~ "0.8–1.0",
        ld_r2 >= 0.6 ~ "0.6–0.8",
        ld_r2 >= 0.4 ~ "0.4–0.6",
        ld_r2 >= 0.2 ~ "0.2–0.4",
        TRUE ~ "0.0–0.2"
      ),
      ld_bin = factor(ld_bin, levels = c("0.8–1.0", "0.6–0.8", "0.4–0.6", "0.2–0.4", "0.0–0.2", "NA"))
    ) %>%
    dplyr::select(-ld_r2_key, -ld_r2_pos, -n_variants_at_pos)
}

report_ld_coverage <- function(df, label) {
  n <- nrow(df)
  mapped <- sum(!is.na(df$ld_r2))
  message(label, " LD coverage: ", mapped, "/", n, " (", sprintf("%.1f", 100 * mapped / max(n, 1)), "%)")
  if ("ld_match_method" %in% names(df)) {
    message("  matching methods: ", paste(names(table(df$ld_match_method)), as.integer(table(df$ld_match_method)), sep = "=", collapse = ", "))
  }
}

# ---------- Gene annotation: LocusZoom-like exon/gene track from local GENCODE GTF ----------
get_gene_track_gtf <- function(chr, region_start, region_end,
                               gtf_file = gene_gtf_file,
                               cache_dir = cache_dir,
                               protein_coding_only = FALSE) {
  chr_num <- gsub("^chr", "", as.character(chr))
  chr_name <- paste0("chr", chr_num)

  cache_file <- file.path(
    cache_dir,
    paste0("genes_gencode_gtf_chr", chr_num, "_", region_start, "_", region_end,
           ifelse(protein_coding_only, "_protein_coding", ""), ".tsv")
  )

  if (file.exists(cache_file)) {
    genes_cached <- data.table::fread(cache_file, data.table = FALSE)
    if (nrow(genes_cached) > 0) return(genes_cached)
    message("Cached GENCODE gene file was empty; refreshing: ", cache_file)
    unlink(cache_file, force = TRUE)
  }

  if (!file.exists(gtf_file)) {
    gtf_url <- if (plot_build == 37L) gene_gtf_url_grch37 else gene_gtf_url_grch38

    if (isTRUE(auto_download_gene_gtf)) {
      dir.create(dirname(gtf_file), recursive = TRUE, showWarnings = FALSE)
      message("GENCODE GTF is missing; downloading ", plot_build_label, " annotation to: ", gtf_file)
      download_status <- tryCatch({
        utils::download.file(gtf_url, destfile = gtf_file, mode = "wb", quiet = FALSE)
        0L
      }, error = function(e) {
        warning("Automatic GENCODE download failed: ", conditionMessage(e))
        1L
      })

      if (!identical(download_status, 0L) || !file.exists(gtf_file) || file.info(gtf_file)$size == 0) {
        stop(
          "GENCODE GTF file is required but could not be downloaded.\n",
          "Expected file: ", gtf_file, "\n",
          "Download URL: ", gtf_url, "\n",
          "Manual command:\n",
          "  curl -L -o ", shQuote(gtf_file), " ", shQuote(gtf_url)
        )
      }
    } else {
      stop(
        "GENCODE GTF file not found: ", gtf_file, "\n",
        "For plot_build = ", plot_build, ", download it with:\n",
        "  curl -L -o ", shQuote(gtf_file), " ", shQuote(gtf_url), "\n",
        "or set auto_download_gene_gtf <- TRUE."
      )
    }
  }


  region <- GenomicRanges::GRanges(
    seqnames = chr_name,
    ranges = IRanges::IRanges(
      start = as.integer(region_start),
      end = as.integer(region_end)
    )
  )

  message("Importing GENCODE gene annotation from local GTF for ",
          chr_name, ":", region_start, "-", region_end)

  gtf <- rtracklayer::import(gtf_file, which = region)
  gtf_df <- as.data.frame(gtf)

  if (nrow(gtf_df) == 0) {
    warning("No GENCODE features found in region ", chr_name, ":", region_start, "-", region_end)
    return(data.frame())
  }

  # GENCODE normally uses gene_type/transcript_type. Some GTFs use gene_biotype.
  if (!"gene_type" %in% names(gtf_df) && "gene_biotype" %in% names(gtf_df)) {
    gtf_df$gene_type <- gtf_df$gene_biotype
  }
  if (!"gene_type" %in% names(gtf_df)) {
    gtf_df$gene_type <- NA_character_
  }
  if (!"gene_name" %in% names(gtf_df)) {
    gtf_df$gene_name <- NA_character_
  }
  if (!"gene_id" %in% names(gtf_df)) {
    gtf_df$gene_id <- NA_character_
  }
  if (!"transcript_id" %in% names(gtf_df)) {
    gtf_df$transcript_id <- NA_character_
  }

  tx <- gtf_df %>%
    dplyr::filter(type == "transcript", !is.na(transcript_id)) %>%
    dplyr::mutate(
      gene_label = dplyr::coalesce(gene_name, gene_id),
      gene_biotype = gene_type,
      strand = ifelse(as.character(strand) == "+", 1, -1),
      tx_start = as.numeric(start),
      tx_end = as.numeric(end),
      tx_len = tx_end - tx_start
    ) %>%
    dplyr::filter(!is.na(gene_label), !is.na(tx_start), !is.na(tx_end))

  exons <- gtf_df %>%
    dplyr::filter(type == "exon", !is.na(transcript_id)) %>%
    dplyr::mutate(
      gene_label = dplyr::coalesce(gene_name, gene_id),
      gene_biotype = gene_type,
      strand = ifelse(as.character(strand) == "+", 1, -1),
      exon_chrom_start = as.numeric(start),
      exon_chrom_end = as.numeric(end)
    ) %>%
    dplyr::filter(!is.na(gene_label), !is.na(exon_chrom_start), !is.na(exon_chrom_end))

  if (protein_coding_only) {
    tx <- tx %>% dplyr::filter(gene_biotype == "protein_coding")
    exons <- exons %>% dplyr::filter(gene_biotype == "protein_coding")
  }

  # If transcript rows are absent in an unusual GTF, construct transcript spans from exons.
  if (nrow(tx) == 0 && nrow(exons) > 0) {
    tx <- exons %>%
      dplyr::group_by(gene_label, gene_id, transcript_id, strand, gene_biotype) %>%
      dplyr::summarise(
        tx_start = min(exon_chrom_start, na.rm = TRUE),
        tx_end = max(exon_chrom_end, na.rm = TRUE),
        tx_len = tx_end - tx_start,
        .groups = "drop"
      )
  }

  if (nrow(exons) == 0 || nrow(tx) == 0) {
    warning("No exon/transcript records found in GENCODE GTF for this region.")
    return(data.frame())
  }

  # One representative transcript per gene: choose the longest transcript.
  tx_keep <- tx %>%
    dplyr::arrange(gene_label, dplyr::desc(tx_len)) %>%
    dplyr::group_by(gene_label) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(tx_start, tx_end)

  # Greedy lane assignment by gene span.
  lane_ends <- numeric(0)
  lanes <- integer(nrow(tx_keep))

  for (i in seq_len(nrow(tx_keep))) {
    placed <- FALSE

    if (length(lane_ends) > 0) {
      for (lane in seq_along(lane_ends)) {
        if (tx_keep$tx_start[i] > lane_ends[lane] + 10000) {
          lanes[i] <- lane
          lane_ends[lane] <- tx_keep$tx_end[i]
          placed <- TRUE
          break
        }
      }
    }

    if (!placed) {
      lane_ends <- c(lane_ends, tx_keep$tx_end[i])
      lanes[i] <- length(lane_ends)
    }
  }

  tx_keep$lane <- lanes

  exons2 <- exons %>%
    dplyr::semi_join(tx_keep %>% dplyr::select(transcript_id), by = "transcript_id") %>%
    dplyr::left_join(
      tx_keep %>%
        dplyr::select(gene_label, transcript_id, tx_start, tx_end, lane),
      by = c("gene_label", "transcript_id")
    ) %>%
    dplyr::distinct(
      gene_label, transcript_id, strand, gene_biotype,
      tx_start, tx_end, exon_chrom_start, exon_chrom_end, lane
    ) %>%
    dplyr::arrange(lane, tx_start, exon_chrom_start)

  data.table::fwrite(exons2, cache_file, sep = "\t")
  exons2
}

plot_gene_track <- function(genes, region_start, region_end,
                            show_labels = show_gene_labels) {
  if (nrow(genes) == 0) {
    return(
      ggplot() +
        theme_void() +
        coord_cartesian(xlim = c(region_start, region_end)) +
        labs(title = NULL) +
        annotate(
          "text",
          x = mean(c(region_start, region_end)),
          y = 1,
          label = "No genes found in local GENCODE GTF",
          size = 3
        )
    )
  }

  genes <- genes %>%
    dplyr::mutate(
      y = -lane,
      label_x = (tx_start + tx_end) / 2,
      gene_label_plot = ifelse(grepl("^ENSG", gene_label), "", gene_label),
      arrow_start = ifelse(strand == 1, tx_start, tx_end),
      arrow_end = ifelse(strand == 1, tx_end, tx_start)
    )

  gene_spans <- genes %>%
    dplyr::group_by(gene_label, gene_label_plot, tx_start, tx_end, strand, lane, y, label_x, arrow_start, arrow_end) %>%
    dplyr::summarise(.groups = "drop")

  p <- ggplot() +
    geom_segment(
      data = gene_spans,
      aes(x = tx_start, xend = tx_end, y = y, yend = y),
      linewidth = 0.28,
      color = "black"
    ) +
    geom_segment(
      data = gene_spans,
      aes(x = arrow_start, xend = arrow_end, y = y, yend = y),
      linewidth = 0.28,
      color = "black",
      arrow = arrow(length = unit(0.045, "inches"), type = "closed")
    ) +
    geom_rect(
      data = genes,
      aes(
        xmin = pmax(exon_chrom_start, region_start),
        xmax = pmin(exon_chrom_end, region_end),
        ymin = y - 0.10,
        ymax = y + 0.10
      ),
      fill = "black",
      color = "black"
    ) +
    coord_cartesian(xlim = c(region_start, region_end), clip = "off") +
    theme_classic(base_size = 10) +
    labs(x = paste0("Position (", plot_build_label, ")"), y = NULL, title = NULL) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_blank(),
      plot.margin = margin(2, 95, 5.5, 5.5)
    )

  if (isTRUE(show_labels)) {
    p <- p +
      geom_text(
        data = gene_spans,
        aes(x = label_x, y = y + 0.27, label = gene_label_plot),
        size = 2.4,
        check_overlap = FALSE
      )
  }

  p
}

# ---------- Plotting ----------
lead_marker_df <- function(df, lead_pos) {
  df <- df %>%
    dplyr::mutate(pos = as.numeric(pos), logp = safe_mlog10(p))

  ymax <- max(df$logp, na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1

  hit <- df %>% dplyr::filter(pos == lead_pos)

  if (nrow(hit) > 0) {
    return(data.frame(
      pos = lead_pos,
      logp = max(hit$logp, na.rm = TRUE),
      marker_note = "observed"
    ))
  }

  # Some conditional pQTL datasets may not contain the fixed LD-reference SNP.
  # Still draw the purple diamond at the fixed reference position so every panel
  # makes clear that LD coloring is relative to the same MR variant.
  data.frame(
    pos = lead_pos,
    logp = ymax * 0.985,
    marker_note = "reference SNP not present in this panel; shown at top"
  )
}

plot_locus_ld <- function(df, title, region_start, region_end,
                          lead_pos,
                          lead_label = ld_reference$label,
                          instruments = get_plot_instruments(),
                          pp4_label = NULL,
                          show_ld_legend = TRUE,
                          show_reference_label = FALSE) {
  df <- df %>%
    dplyr::mutate(
      pos = as.numeric(pos),
      logp = safe_mlog10(p)
    )

  marker <- lead_marker_df(df, lead_pos)

  p <- ggplot(df, aes(pos, logp)) +
    geom_point(
      aes(fill = ld_bin),
      shape = 21,
      size = 1.45,
      alpha = 0.90,
      color = "black",
      stroke = 0.10
    ) +
    geom_vline(
      data = instruments,
      aes(xintercept = pos, color = instrument, linetype = instrument),
      linewidth = 0.65,
      alpha = 0.95
    ) +
    geom_point(
      data = marker,
      aes(pos, logp),
      inherit.aes = FALSE,
      shape = 23,
      size = 3.8,
      fill = "purple",
      color = "black",
      stroke = 0.45
    ) +
    scale_fill_manual(
      values = c(
        "0.8–1.0" = "#d73027",
        "0.6–0.8" = "#fc8d59",
        "0.4–0.6" = "#74c476",
        "0.2–0.4" = "#6baed6",
        "0.0–0.2" = "#2171b5",
        "NA" = "grey70"
      ),
      drop = FALSE
    ) +
    scale_color_manual(
      values = setNames(instruments$color, instruments$instrument)
    ) +
    scale_linetype_manual(
      values = setNames(instruments$line_type, instruments$instrument)
    ) +
    guides(
      fill = guide_legend(
        title = bquote(LD~(r^2)~to~rs28485881),
        override.aes = list(shape = 21, size = 2.8, color = "black"),
        order = 1
      ),
      color = "none",
      linetype = "none"
    ) +
    theme_classic(base_size = 10) +
    coord_cartesian(xlim = c(region_start, region_end), clip = "off") +
    labs(
      title = title,
      x = NULL,
      y = expression(-log[10](italic(P))),
      fill = NULL,
      color = NULL
    ) +
    theme(
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.30),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 10.5, face = "bold"),
      plot.margin = margin(5.5, 145, 5.5, 5.5),
      # Place the LD legend in the right margin, directly below the
      # "Reference conditional pQTL" annotation, like a LocusZoom side legend.
      legend.position = if (show_ld_legend) c(1.035, 0.82) else "none",
      legend.justification = c("left", "top"),
      legend.direction = "vertical",
      legend.background = element_rect(fill = grDevices::adjustcolor("white", alpha.f = 0.00), color = NA),
      legend.box.background = element_blank(),
      legend.title = element_text(size = 8.2),
      legend.text = element_text(size = 7.4),
      legend.key.size = unit(0.38, "lines"),
      legend.key.height = unit(0.38, "lines"),
      legend.spacing.y = unit(0.02, "cm"),
      legend.margin = margin(1, 2, 1, 2)
    )

  if (!is.null(pp4_label)) {
    pp4_x <- region_end + 0.065 * (region_end - region_start)

    p <- p +
      annotate(
        "text",
        x = pp4_x,
        y = Inf,
        label = pp4_label,
        hjust = 0,
        vjust = 1.3,
        size = 3.1
      )
  }

  p
}

# ---------- Preflight LD-panel checks ----------
check_local_pgen(local_ld_pfile_prefix)
build_map <- read_build_map(build_map_file)
regional_build_offset <- estimate_regional_build_offset(build_map)
preflight_ld <- check_mr_variants_in_ld_panel(
  pfile_prefix = local_ld_pfile_prefix,
  build_map = build_map,
  regional_offset = regional_build_offset
)

map_pos38_to37_complete <- function(pos38) as.numeric(pos38) - regional_build_offset
map_pos37_to38_complete <- function(pos37) as.numeric(pos37) + regional_build_offset
map_plot_pos_from38 <- function(pos38) {
  if (plot_build == 38L) as.numeric(pos38) else map_pos38_to37_complete(pos38)
}

# ---------- Load and format HbA1c / T2D ----------
a1c <- data.table::fread(a1c_file)
t2d <- data.table::fread(t2d_file)
if (!file.exists(magic_file)) stop("MAGIC file not found: ", magic_file)
magic <- data.table::fread(magic_file)

a1c2 <- a1c %>%
  dplyr::mutate(
    snp = paste0(CHR, ":", POS),
    pos38 = as.numeric(POS),
    pos37 = map_pos38_to37_complete(POS),
    pos = if (plot_build == 37L) pos37 else pos38,
    key37 = canonical_variant_key(CHR, pos37, Allele1, Allele2),
    beta = Effect,
    varbeta = StdErr^2,
    p = `P-value`,
    maf = pmin(Freq1, 1 - Freq1),
    n_ukbb = ifelse(!is.na(ukbb_beta), dplyr::coalesce(ukbb_num_samples, ukbb_num_cases), 0),
    n_mvp = ifelse(!is.na(mvp_beta), mvp_num_samples, 0),
    N_row = n_ukbb + n_mvp
  ) %>%
  dplyr::filter(!is.na(pos), !is.na(beta), !is.na(varbeta), !is.na(p), !is.na(maf),
                maf > 0, maf < 0.5, N_row > 0) %>%
  diagnose_duplicates("snp", "HbA1c", remove = TRUE)

t2d2 <- t2d %>%
  dplyr::mutate(
    snp = paste0(CHR, ":", POS),
    pos38 = as.numeric(POS),
    pos37 = map_pos38_to37_complete(POS),
    pos = if (plot_build == 37L) pos37 else pos38,
    key37 = canonical_variant_key(CHR, pos37, Allele1, Allele2),
    beta = Effect,
    varbeta = StdErr^2,
    p = `P-value`,
    maf = pmin(Freq1, 1 - Freq1),
    n_ukbb = ifelse(!is.na(ukbb_beta), ukbb_num_samples, 0),
    c_ukbb = ifelse(!is.na(ukbb_beta), ukbb_num_cases, 0),
    n_mvp = ifelse(!is.na(mvp_beta), mvp_num_samples, 0),
    c_mvp = ifelse(!is.na(mvp_beta), mvp_num_cases, 0),
    n_finngen = ifelse(!is.na(finngen_beta), finngen_num_samples, 0),
    c_finngen = ifelse(!is.na(finngen_beta), finngen_num_cases, 0),
    N_row = n_ukbb + n_mvp + n_finngen,
    cases_row = c_ukbb + c_mvp + c_finngen
  ) %>%
  dplyr::filter(!is.na(pos), !is.na(beta), !is.na(varbeta), !is.na(p), !is.na(maf),
                maf > 0, maf < 0.5, N_row > 0) %>%
  diagnose_duplicates("snp", "T2D", remove = TRUE)

magic2 <- magic %>%
  dplyr::mutate(
    pos37 = as.numeric(base_pair_location),
    pos38 = map_pos37_to38_complete(pos37),
    pos = if (plot_build == 37L) pos37 else pos38,
    snp = paste0(chromosome, ":", pos38),
    key37 = canonical_variant_key(chromosome, pos37, effect_allele, other_allele),
    beta = as.numeric(beta),
    varbeta = as.numeric(standard_error)^2,
    p = as.numeric(p_value),
    maf = pmin(as.numeric(effect_allele_frequency), 1 - as.numeric(effect_allele_frequency)),
    N_row = as.numeric(sample_size)
  ) %>%
  dplyr::filter(!is.na(pos), !is.na(pos37), !is.na(pos38), !is.na(beta), !is.na(varbeta),
                !is.na(p), !is.na(maf), maf > 0, maf < 0.5, N_row > 0) %>%
  diagnose_duplicates("snp", magic_trait_label, remove = TRUE)

N_a1c <- median(a1c2$N_row, na.rm = TRUE)
N_t2d <- median(t2d2$N_row, na.rm = TRUE)
s_t2d <- sum(t2d2$cases_row, na.rm = TRUE) / sum(t2d2$N_row, na.rm = TRUE)
N_magic <- median(magic2$N_row, na.rm = TRUE)

# ---------- Conditional pQTL + whole blood eQTL ----------
format_conditional_file <- function(file) {
  x <- data.table::fread(file)
  x$pos38 <- extract_pos38(x)
  x$pos37 <- if ("bp" %in% names(x)) as.numeric(x$bp) else if ("POS" %in% names(x)) as.numeric(x$POS) else map_pos38_to37_complete(x$pos38)
  variant37_source <- if ("SNP" %in% names(x)) x$SNP else if ("SNPID" %in% names(x)) x$SNPID else if ("SNP2" %in% names(x)) x$SNP2 else NA_character_
  if (all(is.na(variant37_source))) {
    x$key37 <- NA_character_
  } else {
    vp <- parse_colon_variant(variant37_source, "conditional-file GRCh37 variant")
    x$key37 <- canonical_variant_key(vp$chr, x$pos37, vp$a1, vp$a2)
  }

  pqtl_c <- x %>%
    dplyr::mutate(
      snp = paste0(chr, ":", pos38),
      pos37 = as.numeric(pos37),
      pos38 = as.numeric(pos38),
      pos = if (plot_build == 37L) pos37 else pos38,
      key37 = key37,
      beta = bC,
      varbeta = bC_se^2,
      p = pC,
      maf_raw = get_col(x, c("fre", "freq", "MAF", "maf"), "conditional pQTL MAF"),
      maf = pmin(maf_raw, 1 - maf_raw),
      N_row = n
    ) %>%
    dplyr::filter(!is.na(beta), !is.na(varbeta), !is.na(p), !is.na(maf),
                  maf > 0, maf < 0.5, N_row > 0) %>%
    diagnose_duplicates("snp", paste0(basename(file), " conditional pQTL"), remove = TRUE)

  eqtl_wb <- x %>%
    dplyr::mutate(
      snp = paste0(chr, ":", pos38),
      pos37 = as.numeric(pos37),
      pos38 = as.numeric(pos38),
      pos = if (plot_build == 37L) pos37 else pos38,
      key37 = key37,
      beta = slope,
      varbeta = slope_se^2,
      p = pval_nominal,
      maf = pmin(maf, 1 - maf),
      N_row = N
    ) %>%
    dplyr::filter(!is.na(beta), !is.na(varbeta), !is.na(p), !is.na(maf),
                  maf > 0, maf < 0.5, N_row > 0) %>%
    diagnose_duplicates("snp", paste0(basename(file), " whole blood eQTL"), remove = TRUE)

  list(
    file = file,
    file_label = basename(file),
    index_snp = unique(x$cojo_snp)[1],
    pqtl = pqtl_c,
    eqtl = eqtl_wb
  )
}

run_coloc <- function(d1, d2, name1, name2, type1, type2,
                      N1, N2, s1 = NULL, s2 = NULL,
                      sdY1 = NULL, sdY2 = NULL) {

  m <- dplyr::inner_join(
    d1 %>% dplyr::select(snp, pos, beta1 = beta, varbeta1 = varbeta, maf1 = maf),
    d2 %>% dplyr::select(snp, beta2 = beta, varbeta2 = varbeta, maf2 = maf),
    by = "snp"
  ) %>%
    dplyr::mutate(maf = ifelse(!is.na(maf1), maf1, maf2)) %>%
    dplyr::filter(!is.na(beta1), !is.na(beta2), !is.na(varbeta1), !is.na(varbeta2))

  D1 <- list(snp = m$snp, position = m$pos, beta = m$beta1,
             varbeta = m$varbeta1, MAF = m$maf, N = N1, type = type1)

  D2 <- list(snp = m$snp, position = m$pos, beta = m$beta2,
             varbeta = m$varbeta2, MAF = m$maf, N = N2, type = type2)

  if (!is.null(s1)) D1$s <- s1
  if (!is.null(s2)) D2$s <- s2
  if (!is.null(sdY1)) D1$sdY <- sdY1
  if (!is.null(sdY2)) D2$sdY <- sdY2

  list(
    name = paste(name1, "vs", name2),
    merged = m,
    result = coloc::coloc.abf(D1, D2)
  )
}

# ---------- Main ----------
all_summary <- list()

plot_cond_files <- select_conditional_files(cond_files, conditional_plot_mode)
message("Conditional plotting mode: ", conditional_plot_mode, "; plotting ", length(plot_cond_files), " conditional file(s).")
message("Show additional MR variant lines: ", show_additional_mr_variants)
message("Show gene labels: ", show_gene_labels)

output_pdf <- paste0("FN3K_conditional_pQTL_whole_blood_HbA1c_T2D_MAGIC2hGlu_INTERVAL_LD_fixed_rs28485881_LocusZoomStyle_no_liftover_build", plot_build, ".pdf")
pdf(
  output_pdf,
  width = 10.5,
  height = ifelse(conditional_plot_mode == "ld_reference_only", 14.0, 14.5)
)

for (file in plot_cond_files) {

  dat <- format_conditional_file(file)

  pqtl_c <- dat$pqtl
  eqtl_wb <- dat$eqtl

  region_start <- max(min(pqtl_c$pos), min(eqtl_wb$pos), min(a1c2$pos), min(t2d2$pos), min(magic2$pos))
  region_end   <- min(max(pqtl_c$pos), max(eqtl_wb$pos), max(a1c2$pos), max(t2d2$pos), max(magic2$pos))

  pqtl_c_region <- pqtl_c %>% dplyr::filter(pos >= region_start, pos <= region_end)
  eqtl_wb_region <- eqtl_wb %>% dplyr::filter(pos >= region_start, pos <= region_end)
  a1c_region <- a1c2 %>% dplyr::filter(pos >= region_start, pos <= region_end)
  t2d_region <- t2d2 %>% dplyr::filter(pos >= region_start, pos <= region_end)
  magic_region <- magic2 %>% dplyr::filter(pos >= region_start, pos <= region_end)

  N_pqtl_c <- median(pqtl_c_region$N_row, na.rm = TRUE)
  N_eqtl_wb <- median(eqtl_wb_region$N_row, na.rm = TRUE)

  main_title <- format_conditioned_label(dat$file_label)

  pairs <- list(
    pqtlC_eqtlWB = run_coloc(
      pqtl_c_region, eqtl_wb_region,
      "conditional pQTL", "whole blood eQTL",
      "quant", "quant",
      N1 = N_pqtl_c, N2 = N_eqtl_wb
    ),
    pqtlC_a1c = run_coloc(
      pqtl_c_region, a1c_region,
      "conditional pQTL", "HbA1c",
      "quant", "quant",
      N1 = N_pqtl_c, N2 = N_a1c,
      sdY2 = 1
    ),
    pqtlC_t2d = run_coloc(
      pqtl_c_region, t2d_region,
      "conditional pQTL", "Type 2 diabetes",
      "quant", "cc",
      N1 = N_pqtl_c, N2 = N_t2d,
      s2 = s_t2d
    ),
    pqtlC_magic = run_coloc(
      pqtl_c_region, magic_region,
      "conditional pQTL", magic_trait_label,
      "quant", "quant",
      N1 = N_pqtl_c, N2 = N_magic,
      sdY2 = 1
    ),
    eqtlWB_a1c = run_coloc(
      eqtl_wb_region, a1c_region,
      "whole blood eQTL", "HbA1c",
      "quant", "quant",
      N1 = N_eqtl_wb, N2 = N_a1c,
      sdY2 = 1
    ),
    eqtlWB_t2d = run_coloc(
      eqtl_wb_region, t2d_region,
      "whole blood eQTL", "Type 2 diabetes",
      "quant", "cc",
      N1 = N_eqtl_wb, N2 = N_t2d,
      s2 = s_t2d
    ),
    a1c_t2d = run_coloc(
      a1c_region, t2d_region,
      "HbA1c", "Type 2 diabetes",
      "quant", "cc",
      N1 = N_a1c, N2 = N_t2d,
      sdY1 = 1, s2 = s_t2d
    )
  )

  # Use the SAME LD reference SNP for all pages/panels, matching instrument 1 MR.
  # INTERVAL LD is calculated in GRCh37; plot coordinates follow plot_build.
  chr_region <- ld_reference$chr
  lead_plot_pos <- ld_reference$plot_pos

  message("Fixed INTERVAL LD reference: ", ld_reference$label)

  ld_tbl <- get_ld_for_region_local_plink(
    chr = chr_region,
    lead_snp37 = ld_reference$snp37,
    region_start = region_start,
    region_end = region_end,
    pfile_prefix = local_ld_pfile_prefix,
    build_map = build_map,
    regional_offset = regional_build_offset,
    plot_build = plot_build,
    plink2 = plink2_bin,
    cache_dir = cache_dir
  )

  pqtl_c_region_ld <- add_ld_to_plot_df(pqtl_c_region, ld_tbl, lead_plot_pos)
  eqtl_wb_region_ld <- add_ld_to_plot_df(eqtl_wb_region, ld_tbl, lead_plot_pos)
  a1c_region_ld <- add_ld_to_plot_df(a1c_region, ld_tbl, lead_plot_pos)
  t2d_region_ld <- add_ld_to_plot_df(t2d_region, ld_tbl, lead_plot_pos)
  magic_region_ld <- add_ld_to_plot_df(magic_region, ld_tbl, lead_plot_pos)

  report_ld_coverage(pqtl_c_region_ld, "Conditional pQTL")
  report_ld_coverage(eqtl_wb_region_ld, "Whole-blood eQTL")
  report_ld_coverage(a1c_region_ld, "HbA1c")
  report_ld_coverage(t2d_region_ld, "T2D")
  report_ld_coverage(magic_region_ld, magic_trait_label)

  genes_region <- get_gene_track_gtf(
    chr = chr_region,
    region_start = region_start,
    region_end = region_end,
    gtf_file = gene_gtf_file,
    cache_dir = cache_dir,
    protein_coding_only = gene_track_protein_coding_only
  )

  locus_plot <-
    plot_locus_ld(
      pqtl_c_region_ld,
      paste0(
        main_title,
        "\nFN3K SomaScan conditional pQTL"
      ),
      region_start,
      region_end,
      lead_pos = lead_plot_pos,
      lead_label = ld_reference$label,
      pp4_label = "Reference conditional pQTL",
      show_ld_legend = TRUE,
      show_reference_label = FALSE
    ) /
    plot_locus_ld(
      eqtl_wb_region_ld,
      "FN3K GTEx whole blood eQTL",
      region_start,
      region_end,
      lead_pos = lead_plot_pos,
      lead_label = ld_reference$label,
      pp4_label = get_pp4_label(pairs, "pqtlC_eqtlWB"),
      show_ld_legend = FALSE
    ) /
    plot_locus_ld(
      a1c_region_ld,
      "HbA1c",
      region_start,
      region_end,
      lead_pos = lead_plot_pos,
      lead_label = ld_reference$label,
      pp4_label = get_pp4_label(pairs, "pqtlC_a1c"),
      show_ld_legend = FALSE
    ) /
    plot_locus_ld(
      t2d_region_ld,
      "Type 2 diabetes",
      region_start,
      region_end,
      lead_pos = lead_plot_pos,
      lead_label = ld_reference$label,
      pp4_label = get_pp4_label(pairs, "pqtlC_t2d"),
      show_ld_legend = FALSE
    ) /
    plot_locus_ld(
      magic_region_ld,
      magic_trait_label,
      region_start,
      region_end,
      lead_pos = lead_plot_pos,
      lead_label = ld_reference$label,
      pp4_label = get_pp4_label(pairs, "pqtlC_magic"),
      show_ld_legend = FALSE
    ) /
    plot_gene_track(
      genes_region,
      region_start,
      region_end,
      show_labels = show_gene_labels
    ) +
    patchwork::plot_layout(
      heights = c(1, 1, 1, 1, 1, ifelse(show_gene_labels, 0.62, 0.42)),
      guides = "keep"
    )

  print(locus_plot)

  summary_one <- data.table::rbindlist(lapply(names(pairs), function(nm) {
    data.frame(
      file = dat$file_label,
      index_snp = dat$index_snp,
      ld_lead = ld_reference$label,
      comparison = nm,
      trait1 = sub(" vs .*", "", pairs[[nm]]$name),
      trait2 = sub(".* vs ", "", pairs[[nm]]$name),
      nsnps = nrow(pairs[[nm]]$merged),
      t(as.data.frame(pairs[[nm]]$result$summary))
    )
  }), fill = TRUE)

  all_summary[[dat$file_label]] <- summary_one
}

dev.off()

summary_table <- data.table::rbindlist(all_summary, fill = TRUE)

data.table::fwrite(
  summary_table,
  "FN3K_conditional_whole_blood_pairwise_coloc_summary.tsv",
  sep = "\t"
)

print(summary_table)

message("Done.")
message("Outputs:")
message("  ", output_pdf)
message("  FN3K_conditional_whole_blood_pairwise_coloc_summary.tsv")
message("  Cached INTERVAL LD/gene files in: ", cache_dir)
