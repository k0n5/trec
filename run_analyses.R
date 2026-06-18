suppressPackageStartupMessages({
  library(TcellExTRECT)
  library(ImmuneLENS)
  library(data.table)
  library(GenomicAlignments)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
parse_opt <- function(args, flag, default) {
  idx <- which(args == flag)
  if (length(idx) == 0L) return(default)
  if (length(args) <= idx) stop("Missing value after ", flag)
  args[idx + 1L]
}

config <- list(
  bam_file        = args[1],
  regions_bed     = args[2],
  windows_tsv     = args[3],
  output_prefix   = args[4],
  sample_name     = args[5],
  run_extrect     = as.logical(parse_opt(args, "--run_extrect", FALSE)),
  run_immunelens  = as.logical(parse_opt(args, "--run_immunelens", FALSE)),
  run_bam_read_counts  = as.logical(parse_opt(args, "--run_bam_read_counts", FALSE)),
  run_bam_mapping_stats  = as.logical(parse_opt(args, "--run_bam_mapping_stats", FALSE)),
  run_find_trec_krec = as.logical(parse_opt(args, "--run_find_trec_krec", FALSE)),
  make_circle_funnel = as.logical(parse_opt(args, "--make_circle_funnel", FALSE)),
  run_find_dup_del = as.logical(parse_opt(args, "--run_find_dup_del", FALSE)),
  min_mapq        = as.integer(parse_opt(args, "--min_mapq", 40L)),
  min_distance    = as.integer(parse_opt(args, "--min_distance", 1000L)),
  max_nm          = as.integer(parse_opt(args, "--max_nm", 15)),
  max_unaligned_frac = as.numeric(parse_opt(args, "--max_unaligned_frac", 0.20))
)


if (config$make_circle_funnel & !config$run_find_trec_krec) stop("If make_circle_funnel is set, run_find_trec_krec must also be set.")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Recognises only canonical human contigs (1..22, X, Y, M, MT). Non-canonical
# names (decoys, alts, HLA contigs) are left untouched and warned about 
reads_to_chr_style <- function(reads) {
  stopifnot(
    "reads must be a data.frame/tibble" = is.data.frame(reads),
    "reads must contain `rname`"        = "rname" %in% names(reads)
  )
  
  canonical_no_prefix <- c(as.character(1:22), "X", "Y")
  
  # Vectorised renamer for a character vector of chromosome names.
  rename_chr_vec <- function(x) {
    if (length(x) == 0L) return(x)
    out <- x
    needs <- !is.na(x) & !startsWith(x, "chr")
    # MT / M -> chrM
    mt_idx <- needs & x %in% c("MT", "M")
    out[mt_idx] <- "chrM"
    # 1..22 / X / Y -> chr1..chr22 / chrX / chrY
    std_idx <- needs & x %in% canonical_no_prefix
    out[std_idx] <- paste0("chr", x[std_idx])
    # Track what was left alone (non-canonical, non-NA, non-chr)
    untouched <- needs & !mt_idx & !std_idx
    if (any(untouched)) {
      examples <- unique(x[untouched])
      warning(sprintf(
        "Left %d non-canonical chromosome name(s) untouched: %s",
        length(examples),
        paste(head(examples, 5), collapse = ", ")))
    }
    out
  }
  
  # --- rname ---------------------------------------------------------------
  reads$rname <- rename_chr_vec(as.character(reads$rname))
  
  # --- mrnm (mate ref name) ------------------------------------------------
  if ("mrnm" %in% names(reads)) {
    reads$mrnm <- rename_chr_vec(as.character(reads$mrnm))
  }
  
  # --- sa tag (SA:Z) -------------------------------------------------------
  # Format per SAM spec: "rname,pos,strand,CIGAR,mapQ,NM;" repeated.
  # We rewrite only the leading rname field of each semicolon-delimited record.
  if ("sa" %in% names(reads) && !is.null(reads$sa)) {
    sa_chr <- as.character(reads$sa)
    has_sa <- !is.na(sa_chr) & nzchar(sa_chr)
    
    if (any(has_sa)) {
      # Split each read's SA string on ';', drop empty trailing fragments,
      # rewrite the first comma-delimited field, then re-paste.
      sa_fixed <- vapply(sa_chr[has_sa], function(s) {
        parts <- strsplit(s, ";", fixed = TRUE)[[1]]
        parts <- parts[nzchar(parts)]
        rewritten <- vapply(parts, function(rec) {
          fields <- strsplit(rec, ",", fixed = TRUE)[[1]]
          if (length(fields) < 1L) return(rec)
          fields[1] <- rename_chr_vec(fields[1])
          paste(fields, collapse = ",")
        }, character(1), USE.NAMES = FALSE)
        # Preserve trailing ';' which SAM spec mandates
        paste0(paste(rewritten, collapse = ";"), ";")
      }, character(1), USE.NAMES = FALSE)
      sa_chr[has_sa] <- sa_fixed
    }
    reads$sa <- sa_chr
  }
  
  reads
}

cigar_read_length <- function(cigar) {
  stopifnot(is.character(cigar))
  suppressWarnings(cigarWidthAlongQuerySpace(cigar, before.hard.clipping = TRUE))
}

cigar_aligned_bases <- function(cigar) {
  stopifnot(is.character(cigar))
  suppressWarnings({
    ops_list    <- explodeCigarOps(cigar)
    widths_list <- explodeCigarOpLengths(cigar)
  })
  mapply(function(ops, widths) {
    sum(widths[ops %in% c("M", "I", "=", "X")])
  }, ops_list, widths_list, USE.NAMES = FALSE)
}

cigar_clips <- function(cigar) {
  stopifnot(is.character(cigar))
  suppressWarnings({
    ops_list    <- explodeCigarOps(cigar)
    widths_list <- explodeCigarOpLengths(cigar)
  })
  result <- mapply(function(ops, widths) {
    n <- length(ops)
    lead <- 0L
    for (i in seq_len(n)) {
      if (ops[i] %in% c("S", "H")) lead <- lead + widths[i] else break
    }
    trail <- 0L
    for (i in rev(seq_len(n))) {
      if (ops[i] %in% c("S", "H")) trail <- trail + widths[i] else break
    }
    c(lead, trail)
  }, ops_list, widths_list)
  tibble(lead_clip = result[1, ], trail_clip = result[2, ])
}

# Per-CIGAR: the [start, end] interval on the original (pre-hard-clip) read
# that this alignment covers.
cigar_query_interval <- function(cigar) {
  stopifnot(is.character(cigar))
  clips   <- cigar_clips(cigar)            # your existing helper: lead, trail
  aligned <- cigar_aligned_bases(cigar)    # M/I/=/X
  tibble(
    q_start = clips$lead_clip + 1L,
    q_end   = clips$lead_clip + aligned
  )
}

# Truly unaligned bases per read = read_length - union of query intervals
unaligned_frac_primary_supplementary_pair <- function(cigar_p, cigar_s) {
  stopifnot(length(cigar_p) == length(cigar_s))
  
  read_len   <- cigar_read_length(cigar_p)
  read_len_s <- cigar_read_length(cigar_s)
  if (!all(read_len == read_len_s)) {
    mismatched <- which(read_len != read_len_s)
    stop(sprintf(
      "Primary and supplementary CIGAR imply different original read lengths for %d records (e.g. row %d: %d vs %d).",
      length(mismatched), mismatched[1],
      read_len[mismatched[1]], read_len_s[mismatched[1]]
    ))
  }
  
  qi_p <- cigar_query_interval(cigar_p)
  qi_s <- cigar_query_interval(cigar_s)
  
  width_p <- qi_p$q_end - qi_p$q_start + 1L
  width_s <- qi_s$q_end - qi_s$q_start + 1L
  
  # Overlap of two intervals: max(0, min(ends) - max(starts) + 1)
  overlap <- pmax(0L,
                  pmin(qi_p$q_end,   qi_s$q_end) -
                    pmax(qi_p$q_start, qi_s$q_start) + 1L)
  
  covered <- width_p + width_s - overlap
  
  stopifnot(
    "covered bases exceed read length" = all(covered <= read_len),
    "covered bases negative"           = all(covered >= 0L)
  )
  
  unaligned_frac  = (read_len - covered) / read_len
  return(unaligned_frac)
}

annotate_region <- function(chr, pos, regions_gr, sep = "_") {
  
  stopifnot(
    "chr and pos must have equal length"  = length(chr) == length(pos),
    "regions_gr must be a GRanges"        = inherits(regions_gr, "GRanges"),
    "pos must be integer-like"            = is.numeric(pos) | is.na(pos),
    "chr must be character"               = is.character(chr) | is.na(chr)
  )
  
  # Sanity: chromosome naming must match between query and regions
  query_chrs <- unique(chr)
  region_chrs <- unique(as.character(seqnames(regions_gr)))
  if (length(intersect(query_chrs, region_chrs)) == 0L) {
    stop("No shared seqnames between query coordinates and regions. ",
         "Query uses e.g. '", head(query_chrs, 1),
         "', regions use e.g. '", head(region_chrs, 1),
         "'. Check chr-prefix.")
  }
  
  query_gr <- GRanges(
    seqnames = chr,
    ranges   = IRanges(start = as.integer(pos), width = 1L)
  )
  hits <- GenomicRanges::findOverlaps(query_gr, regions_gr,
                                      ignore.strand = TRUE)
  
  # For each query, collapse all overlapping region names into one string
  result <- rep("", length(chr))
  if (length(hits) > 0L) {
    qh <- queryHits(hits)
    sh <- subjectHits(hits)
    # Collapse per query index
    by_query <- split(regions_gr$name[sh], qh)
    # Sort within each so output is deterministic regardless of region order
    collapsed <- vapply(by_query,
                        function(x) paste(sort(unique(x)), collapse = sep),
                        character(1))
    result[as.integer(names(by_query))] <- collapsed
  }
  
  result
}

# ============================================================================
# READ BED
# ============================================================================
regions_gr <- rtracklayer::import(config$regions_bed, format = "BED")

# ============================================================================
# READ BAM
# ============================================================================
if (config$run_bam_read_counts | config$run_bam_mapping_stats | config$run_find_trec_krec | config$run_find_dup_del){
  stopifnot("BAM file does not exist"  = file.exists(config$bam_file))
  stopifnot("windows_tsv file does not exist"  = file.exists(config$windows_tsv))
  # Ensure output directory exists
  out_dir <- dirname(config$output_prefix)
  if (!dir.exists(out_dir) && out_dir != ".") {
    dir.create(out_dir, recursive = TRUE)
  }
  
  cat(sprintf("BAM:          %s\n", config$bam_file))
  cat(sprintf("WINDOWS:      %s\n", config$windows_tsv))
  cat(sprintf("Output:       %s.*\n", config$output_prefix))
  cat(sprintf("Min mapQ:     %d\n", config$min_mapq))
  cat(sprintf("Min distance: %d bp\n", config$min_distance))
  cat(sprintf("Max NM:       %d\n", config$max_nm))
  cat(sprintf("Max unaligned: %.2f\n", config$max_unaligned_frac))
  cat("\n")
  
  cat("Reading BAM ...\n")
  
  bam_param <- ScanBamParam(
    what = c("qname", "flag", "rname", "pos", "mapq", "cigar",
             "mrnm", "mpos", "isize"),
    tag  = c("NM", "SA"),
    flag = scanBamFlag(
      isUnmappedQuery      = FALSE,
      isSecondaryAlignment = FALSE,
      isDuplicate = FALSE
    )
  )
  
  raw <- scanBam(BamFile(config$bam_file), param = bam_param)[[1]]
  
  reads <- tibble(
    qname  = raw$qname,
    flag   = raw$flag,
    rname  = as.character(raw$rname),
    pos    = raw$pos,
    mapq   = raw$mapq,
    cigar  = raw$cigar,
    mrnm   = as.character(raw$mrnm),
    mpos   = raw$mpos,
    isize  = raw$isize,
    nm     = raw$tag$NM,
    sa     = raw$tag$SA
  ) %>%
    mutate(
      is_paired        = bamFlagTest(flag, "isPaired"),
      is_r1            = bamFlagTest(flag, "isFirstMateRead"),
      is_r2            = bamFlagTest(flag, "isSecondMateRead"),
      is_supplementary = bamFlagTest(flag, "isSupplementaryAlignment"),
      strand           = ifelse(bamFlagTest(flag, "isMinusStrand"), "-", "+"),
      mate_strand      = ifelse(bamFlagTest(flag, "isMateMinusStrand"), "-", "+")
    )
  
  if (nrow(reads)>0 & !grepl("chr", reads$rname[1])){
    reads = reads_to_chr_style(reads)
  }
  
  rm(raw); gc()
  cat(sprintf("  Loaded %d alignment records\n", nrow(reads)))
  
  # --- Remove reads with multiple supplementary alignments on either mate ---
  multi_suppl_qnames <- reads %>%
    filter(is_supplementary) %>%
    group_by(qname, is_r1) %>%
    summarise(n_suppl = dplyr::n(), .groups = "drop") %>%
    filter(n_suppl > 1L) %>%
    pull(qname) %>%
    unique()
  
  if (length(multi_suppl_qnames) > 0L) {
    n_before <- nrow(reads)
    reads <- reads %>% filter(!qname %in% multi_suppl_qnames)
    cat(sprintf("  Removed %d records (%d qnames) with multiple supplementary alignments\n",
                n_before - nrow(reads), length(multi_suppl_qnames)))
  }
}

# ------------------------------------------------------------------------------
# run_bam_read_counts
# ------------------------------------------------------------------------------
if (config$run_bam_read_counts){
  
  # --- Build read GRanges ---------------------------------------------------
  suppressWarnings({
    reads_gr <- GRanges(
      seqnames = reads$rname,
      ranges   = IRanges(
        start = reads$pos,
        width = GenomicAlignments::cigarWidthAlongReferenceSpace(reads$cigar)
      ),
      strand   = "*",
      qname    = reads$qname,
      mapq     = reads$mapq,
      nm       = reads$nm,
      is_supplementary = reads$is_supplementary
    )
  })
  
  # --- Sanity check: chromosome naming convention ---------------------------
  common_chr <- intersect(as.character(seqnames(reads_gr)),
                          as.character(seqnames(regions_gr)))
  if (length(common_chr) == 0L) {
    stop("No shared seqnames between BAM reads and regions BED. ",
         "BAM uses e.g. '", head(unique(as.character(seqnames(reads_gr))), 1),
         "', BED uses e.g. '", head(unique(as.character(seqnames(regions_gr))), 1),
         "'. Fix chromosome naming (chr-prefix mismatch?) before continuing.")
  }
  
  # --- Count overlaps -------------------------------------------------------
  hits <- GenomicRanges::findOverlaps(reads_gr, regions_gr,
                                      ignore.strand = TRUE)
  
  if (length(hits) == 0L) {
    warning("No reads overlap any region. Check coordinates and chromosome naming.")
  }
  
  reads_in_regions <- tibble(
    region        = regions_gr$name[subjectHits(hits)],
    qname         = reads_gr$qname[queryHits(hits)],
    is_supplementary = reads_gr$is_supplementary[queryHits(hits)]
  )
  
  # --- Summarise: alignments, unique reads, post-filter unique reads --------
  summary_per_region <- reads_in_regions %>%
    group_by(region) %>%
    summarise(
      n_unique_reads   = dplyr::n_distinct(qname),
      n_unique_primary = dplyr::n_distinct(qname[!is_supplementary]),
      .groups = "drop"
    )
  
  # Add regions with zero overlapping reads back in
  summary_per_region <- 
    tibble(region = regions_gr$name) %>%
    left_join(summary_per_region, by = "region") %>%
    mutate(across(starts_with("n_"), ~ tidyr::replace_na(., 0L)))
  
  # Write output 
  out_file <- paste0(config$output_prefix, ".region_read_counts.tsv")
  data.table::fwrite(summary_per_region, out_file, sep = "\t")
  
}

# ------------------------------------------------------------------------------
# run_bam_mapping_stats
# ------------------------------------------------------------------------------
if (config$run_bam_mapping_stats){
  out_file <- paste0(config$output_prefix, ".properpair_mapping_stats.tsv")
  reads %>%
    dplyr::filter(bamFlagTest(flag, "isProperPair")) %>%
    dplyr::filter(!is_supplementary) %>%
    dplyr::mutate(unaligned_frac = 1 - ((cigar_aligned_bases(cigar) / cigar_read_length(cigar)))) %>%
    dplyr::summarise(
      mean_isize = mean(abs(isize), na.rm = TRUE),
      sd_isize = sd(abs(isize), na.rm = TRUE),
      median_isize = median(abs(isize), na.rm = TRUE),
      q95_isize = quantile(abs(isize), 0.95, na.rm = TRUE),
      q99_isize = quantile(abs(isize), 0.99, na.rm = TRUE),
      mean_mapq = mean(mapq, na.rm = TRUE),
      sd_mapq = sd(mapq, na.rm = TRUE),
      median_mapq = median(mapq, na.rm = TRUE),
      q5_mapq = quantile(mapq, 0.05, na.rm = TRUE),
      q1_mapq = quantile(mapq, 0.01, na.rm = TRUE),
      mean_nm = mean(nm, na.rm = TRUE),
      sd_nm = sd(nm, na.rm = TRUE),
      median_nm = median(nm, na.rm = TRUE),
      q95_nm = quantile(nm, 0.95, na.rm = TRUE),
      q99_nm = quantile(nm, 0.99, na.rm = TRUE),
      mean_unaligned_frac = mean(unaligned_frac, na.rm = TRUE),
      sd_unaligned_frac = sd(unaligned_frac, na.rm = TRUE),
      median_unaligned_frac = median(unaligned_frac, na.rm = TRUE),
      q95_unaligned_frac = quantile(unaligned_frac, 0.95, na.rm = TRUE),
      q99_unaligned_frac = quantile(unaligned_frac, 0.99, na.rm = TRUE)
    ) %>%
    data.table::fwrite(out_file, sep = "\t")
}

# ------------------------------------------------------------------------------
# find_trec_krec
# ------------------------------------------------------------------------------
if (config$run_find_trec_krec){
  
  # ============================================================================
  # REARRANGEMENT CHECK
  # ============================================================================
  are_combinable_vec <- function(type1, group1, type2, group2) {
    same_group <- !is.na(group1) & !is.na(group2) & (group1 == group2)
    types_sorted <- ifelse(type1 <= type2,
                           paste0(type1, "-", type2),
                           paste0(type2, "-", type1))
    valid_combo <- types_sorted %in% c("D-J", "D-V", "J-V")
    same_group & valid_combo
  }
  
  # ============================================================================
  # FIND CIRCLE-SUPPORTING DISCORDANT PAIRS
  # ============================================================================
  find_circle_discordant_pairs <- function(reads, config) {
    cat("Scanning for circle-supporting discordant pairs ...\n")
    
    if (nrow(reads)==0){
      cat(sprintf("  find_circle_discordant_pairs: reads is empty, exiting and return empty table."))
      return(tibble())
    }
    
    primary <- reads %>% dplyr::filter(is_paired, !is_supplementary)
    r1 <- primary %>% dplyr::filter(is_r1)
    r1$unaligned_frac = 1-(cigar_aligned_bases(r1$cigar) / cigar_read_length(r1$cigar))
    r2 <- primary %>% dplyr::filter(is_r2)
    r2$unaligned_frac = 1-(cigar_aligned_bases(r2$cigar) / cigar_read_length(r2$cigar))
    
    stopifnot("Duplicate R1 primary qnames" = !any(duplicated(r1$qname)))
    stopifnot("Duplicate R2 primary qnames" = !any(duplicated(r2$qname)))
    
    pairs <- inner_join(r1, r2, by = "qname", suffix = c(".r1", ".r2"))
    
    result <- pairs %>%
      filter(rname.r1 == rname.r2) %>% 
      mutate(
        r1r2_distance = abs(pos.r1 - pos.r2),
        unaligned_frac = (unaligned_frac.r1 + unaligned_frac.r2)/2
      ) %>% 
      mutate(
        is_circle = ifelse(
          pos.r1 < pos.r2,
          strand.r1 == "-" & strand.r2 == "+",
          strand.r1 == "+" & strand.r2 == "-"
        )
      ) %>%
      filter(is_circle)
    
    result$region.r1 = annotate_region(result$rname.r1, result$pos.r1, regions_gr)
    result$region.r2 = annotate_region(result$rname.r2, result$pos.r2, regions_gr)
    
    cat(sprintf("  Found %d circle-supporting discordant pairs\n", nrow(result)))
    
    result
  }
  
  # ============================================================================
  # FIND CIRCLE-SUPPORTING SPLIT READS
  # ============================================================================
  
  find_circle_split_reads <- function(reads, config) {
    cat("Scanning for circle-supporting split reads ...\n")
    
    if (nrow(reads)==0){
      cat(sprintf("  find_circle_split_reads: reads is empty, exiting and return empty table."))
      return(tibble())
    }
    
    primary <- reads %>%
      filter(is_paired, !is_supplementary) %>%
      dplyr::select(qname, rname, pos, mapq, strand, cigar,
                    mrnm, mpos, mate_strand, nm, is_r1, is_r2)
    
    suppl <- reads %>%
      filter(is_paired, is_supplementary) %>%
      dplyr::select(qname, rname, pos, mapq, strand, cigar, nm, is_r1, is_r2)
    
    suppl_counts <- suppl %>%
      group_by(qname, is_r1) %>%
      summarise(n_suppl = dplyr::n(), .groups = "drop") %>%
      filter(n_suppl == 1L)
    suppl <- semi_join(suppl, suppl_counts, by = c("qname", "is_r1"))
    
    paired <- inner_join(primary, suppl,
                         by = c("qname", "is_r1", "is_r2"),
                         suffix = c(".p", ".s"))
    
    cat(sprintf("  Candidate primary-supplementary pairs: %d\n", nrow(paired)))
    paired <- paired %>%
      filter(
        rname.p == rname.s, 
        strand.p == strand.s
      ) %>%
      mutate(
        ps_distance = abs(pos.p - pos.s)
      )
    cat(sprintf("  After geometry filters: %d\n", nrow(paired)))
    if (nrow(paired) == 0L) return(tibble())
    
    # check that each read only has one primary-supplementary pair
    pair_check <- paired %>%
      group_by(qname, is_r1) %>%
      summarise(n_rows = dplyr::n(), .groups = "drop")
    stopifnot(
      "Each (qname, is_r1) in `paired` must have exactly one primary-supplementary pair" =
        all(pair_check$n_rows == 1L)
    )
    
    # Unaligned fraction for this read
    paired$unaligned_frac.ps = unaligned_frac_primary_supplementary_pair(paired$cigar.p, paired$cigar.s)
    
    # First-on-read via both CIGARs
    p_clips <- cigar_clips(paired$cigar.p)
    s_clips <- cigar_clips(paired$cigar.s)
    
    paired <- paired %>%
      mutate(
        p_query_start = ifelse(strand.p == "+", p_clips$lead_clip, p_clips$trail_clip),
        s_query_start = ifelse(strand.s == "+", s_clips$lead_clip, s_clips$trail_clip),
        primary_is_first = p_query_start < s_query_start
      ) %>%
      filter(p_query_start != s_query_start)
    cat(sprintf("  After resolving first-on-read: %d\n", nrow(paired)))
    if (nrow(paired) == 0L) return(tibble())
    
    # Circle signature
    paired <- paired %>%
      mutate(
        first_pos  = ifelse(primary_is_first, pos.p, pos.s),
        second_pos = ifelse(primary_is_first, pos.s, pos.p),
        is_circle  = case_when(
          strand.p == "+" ~ first_pos > second_pos,
          strand.p == "-" ~ first_pos < second_pos
        )
      ) %>%
      filter(is_circle)
    cat(sprintf("  With circle/tandem-dup signature: %d\n", nrow(paired)))
    if (nrow(paired) == 0L) return(tibble())
    
    # Mate consistency check
    suppressWarnings({
      ref_end_p <- paired$pos.p + cigarWidthAlongReferenceSpace(paired$cigar.p) - 1L
      ref_end_s <- paired$pos.s + cigarWidthAlongReferenceSpace(paired$cigar.s) - 1L
    })
    circle_start <- pmin(paired$pos.p, paired$pos.s)
    circle_end   <- pmax(ref_end_p, ref_end_s)
    paired <- paired %>%
      mutate(
        circle_start = circle_start,
        circle_end   = circle_end,
        mate_consistent =
          (mrnm == rname.p) &
          (mate_strand != strand.p) &
          (mpos >= circle_start) &
          (mpos <= circle_end)
      ) %>%
      filter(mate_consistent | is.na(mrnm))
    cat(sprintf("  After mate-consistency check: %d\n", nrow(paired)))
    
    # Mate unaligned fraction
    # Mate primary (may be missing if mate is unmapped or otherwise absent)
    mate_primary <- reads %>%
      dplyr::filter(is_paired, !is_supplementary) %>%
      dplyr::select(qname, is_r1, cigar.mate = cigar, mapq.mate = mapq, nm.mate = nm)
    
    paired <- paired %>%
      mutate(mate_is_r1 = !is_r1) %>%
      left_join(
        mate_primary,
        by = c("qname" = "qname", "mate_is_r1" = "is_r1")
      ) %>%
      mutate(has_mate = !is.na(cigar.mate)) # mark singleton reads
    
    # Unaligned fraction on mate: NA when mate primary is absent
    paired$unaligned_frac.mate <- NA_real_
    if (any(paired$has_mate)) {
      cig <- paired$cigar.mate[paired$has_mate]
      paired$unaligned_frac.mate[paired$has_mate] <-
        1 - cigar_aligned_bases(cig) / cigar_read_length(cig)
    }
    
    paired$unaligned_frac = ifelse(is.na(paired$unaligned_frac.mate),
                                   paired$unaligned_frac.ps,
                                   (paired$unaligned_frac.ps + paired$unaligned_frac.mate) / 2)
    
    # Annotate regions
    paired$region.p = annotate_region(paired$rname.p, paired$pos.p, regions_gr)
    paired$region.s = annotate_region(paired$rname.s, paired$pos.s, regions_gr)
    paired$region.mate <- NA_character_
    if (any(paired$has_mate)) {
      mate_chr <- paired$mrnm[paired$has_mate]
      mate_pos <- paired$mpos[paired$has_mate]
      paired$region.mate[paired$has_mate] <-
        annotate_region(mate_chr, mate_pos, regions_gr)
    }
    
    paired
  }
  
  # ============================================================================
  # RUN CIRCLE DETECTION
  # ============================================================================
  discordant  <- find_circle_discordant_pairs(reads, config)
  split_reads <- find_circle_split_reads(reads, config)
  
  # ============================================================================
  # FILTER DISCORDANT AND SPLIT READS USING DEFAULTS 
  # ============================================================================
  cat(sprintf("  Discordant reads before default filtering: %d\n", nrow(discordant)))
  cat(sprintf("  Split reads before default filtering: %d\n", nrow(split_reads)))
  
  if (file.exists(paste0(config$output_prefix, ".discordant_reads.tsv"))) file.remove(paste0(config$output_prefix, ".discordant_reads.tsv"))
  if (nrow(discordant)>0){
    discordant %>%
      dplyr::select(
        -mrnm.r1, -mpos.r1, -isize.r1, -mate_strand.r1, -is_r1.r1, -is_r2.r1, -is_paired.r1, -is_supplementary.r1,
        -mrnm.r2, -mpos.r2, -isize.r2, -mate_strand.r2, -is_r1.r2, -is_r2.r2, -is_paired.r2, -is_supplementary.r2,
        -is_circle
      ) %>%
      data.table::fwrite(
        paste0(config$output_prefix, ".discordant_reads.tsv"),
        sep = "\t", quote = F, row.names = F, col.names = T
      )
    discordant = 
      discordant %>%
      filter(
        pmin(mapq.r1,mapq.r2) >= config$min_mapq,
        (nm.r1+nm.r2) <= config$max_nm,
        r1r2_distance >= config$min_distance,
        unaligned_frac <= config$max_unaligned_frac
      )
  } else {
    data.table::fwrite(
      tibble(),
      paste0(config$output_prefix, ".discordant_reads.tsv"),
      sep = "\t", quote = F, row.names = F, col.names = T
    )
  }
  
  if (file.exists(paste0(config$output_prefix, ".split_reads.tsv"))) file.remove(paste0(config$output_prefix, ".split_reads.tsv"))
  if (nrow(split_reads)>0){
    split_reads %>%
      dplyr::select(
        -is_circle, -mate_consistent
      ) %>%
      data.table::fwrite(
        paste0(config$output_prefix, ".split_reads.tsv"),
        sep = "\t", quote = F, row.names = F, col.names = T 
      )
    split_reads = 
      split_reads %>%
      filter(
        pmin(mapq.p, mapq.s) >= config$min_mapq,
        is.na(mapq.mate) | (mapq.mate >= config$min_mapq),
        (nm.p + nm.s + ifelse(is.na(nm.mate), 0L, nm.mate)) <= config$max_nm,
        ps_distance >= config$min_distance,
        unaligned_frac <= config$max_unaligned_frac
      )
  } else {
    data.table::fwrite(
      tibble(),
      paste0(config$output_prefix, ".split_reads.tsv"),
      sep = "\t", quote = F, row.names = F, col.names = T 
    )
  }
  
  cat(sprintf("  Discordant reads after default filtering: %d\n", nrow(discordant)))
  cat(sprintf("  Split reads after default filtering: %d\n", nrow(split_reads)))
  
  # ============================================================================
  # UNIFY BREAKPOINTS FOR ANNOTATION
  # ============================================================================
  # For each circle-supporting read, extract the two breakpoint positions.
  # Discordant pairs: breakpoints are at the alignment positions of each mate.
  # Split reads: breakpoints are at the two alignment positions (primary + suppl).
  # We use the full alignment span for overlap with RSS windows.
  
  cat("\nAnnotating circle-supporting reads ...\n")
  breakpoints <- bind_rows(
    # --- Discordant pairs ---
    {
      if (nrow(discordant) > 0L) {
        suppressWarnings({
          ref_end_r1 <- discordant$pos.r1 +
            cigarWidthAlongReferenceSpace(discordant$cigar.r1) - 1L
          ref_end_r2 <- discordant$pos.r2 +
            cigarWidthAlongReferenceSpace(discordant$cigar.r2) - 1L
        })
        discordant %>%
          transmute(
            qname,
            evidence   = "discordant",
            chr        = rname.r1,
            bp1_start  = pos.r1,
            bp1_end    = ref_end_r1,
            bp2_start  = pos.r2,
            bp2_end    = ref_end_r2
          )
      } else {
        tibble()
      }
    },
    # --- Split reads ---
    {
      if (nrow(split_reads) > 0L) {
        suppressWarnings({
          ref_end_p <- split_reads$pos.p +
            cigarWidthAlongReferenceSpace(split_reads$cigar.p) - 1L
          ref_end_s <- split_reads$pos.s +
            cigarWidthAlongReferenceSpace(split_reads$cigar.s) - 1L
        })
        split_reads %>%
          transmute(
            qname,
            evidence   = "split",
            chr        = rname.p,
            bp1_start  = pos.p,
            bp1_end    = ref_end_p,
            bp2_start  = pos.s,
            bp2_end    = ref_end_s
          )
      } else {
        tibble()
      }
    }
  )
  
  if (nrow(breakpoints) == 0L) {
    
    cat("No circle-supporting reads found. Writing empty outputs.\n")
    # Write empty summary
    summary_out <- tibble(
      sample = basename(config$bam_file),
      trec_a   = 0L,
      trec_b   = 0L,
      trec_d   = 0L,
      trec_g   = 0L,
      krec   = 0L,
      igh_circle   = 0L,
      igl_circle   = 0L,
      other  = 0L,
      total  = 0L,
    )
    if (file.exists(paste0(config$output_prefix, ".circle_summary.tsv"))) file.remove(paste0(config$output_prefix, ".circle_summary.tsv"))
    write.table(summary_out,
                file = paste0(config$output_prefix, ".circle_summary.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    if (file.exists(paste0(config$output_prefix, ".circle_reads.tsv"))) file.remove(paste0(config$output_prefix, ".circle_reads.tsv"))
    write.table(tibble(),
                file = paste0(config$output_prefix, ".circle_reads.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
  } else {
    
    # Prefer split-read evidence when a qname has both, make sure every read is only counted once
    breakpoints <- breakpoints %>%
      arrange(qname, desc(evidence == "split")) %>%
      distinct(qname, .keep_all = TRUE)
    
    cat(sprintf("  %d unique circle-supporting reads to annotate\n", nrow(breakpoints)))
    
    # ============================================================================
    # OVERLAP BREAKPOINTS WITH RSS WINDOWS
    # ============================================================================
    # For each read, check if bp1 and bp2 each overlap an RSS window, and whether
    # the two windows are from combinable segments (V-J, V-D, D-J) in the same
    # locus group.
    
    cat("Reading RSS Windows ...\n")
    windows_df = read.table(
      config$windows_tsv, 
      sep = "\t", header = FALSE, stringsAsFactors = FALSE,
      col.names = c("chr", "start", "end", "segment_name", "seg_type", "group", "side")
    ) 
    windows_gr <- makeGRangesFromDataFrame(windows_df, keep.extra.columns = TRUE)
    
    present_chrs <- unique(reads$rname)
    windows_gr <- windows_gr[as.character(seqnames(windows_gr)) %in% present_chrs]
    
    cat("Overlap circle reads with RSS Windows ...\n")
    bp1_gr <- GRanges(seqnames = breakpoints$chr,
                      ranges = IRanges(start = breakpoints$bp1_start,
                                       end   = breakpoints$bp1_end))
    bp2_gr <- GRanges(seqnames = breakpoints$chr,
                      ranges = IRanges(start = breakpoints$bp2_start,
                                       end   = breakpoints$bp2_end))
    
    hits1 <- findOverlaps(bp1_gr, windows_gr)
    hits2 <- findOverlaps(bp2_gr, windows_gr)
    
    hit_tbl1 <- tibble(
      read_idx     = queryHits(hits1),
      seg_type.bp1 = windows_gr$seg_type[subjectHits(hits1)],
      group.bp1    = windows_gr$group[subjectHits(hits1)],
      segment.bp1  = windows_gr$segment_name[subjectHits(hits1)]
    )
    
    hit_tbl2 <- tibble(
      read_idx     = queryHits(hits2),
      seg_type.bp2 = windows_gr$seg_type[subjectHits(hits2)],
      group.bp2    = windows_gr$group[subjectHits(hits2)],
      segment.bp2  = windows_gr$segment_name[subjectHits(hits2)]
    )
    
    # Cross-join within each read: check all (bp1_window, bp2_window) pairs
    combos <- inner_join(hit_tbl1, hit_tbl2, by = "read_idx", relationship = "many-to-many")
    combos <- combos %>%
      mutate(combinable = are_combinable_vec(seg_type.bp1, group.bp1,
                                             seg_type.bp2, group.bp2)) %>%
      filter(combinable)
    
    # For each read, pick one valid combination (prefer V-J over V-D / D-J)
    if (nrow(combos) > 0L) {
      combos <- combos %>%
        mutate(
          combo_type = ifelse(seg_type.bp1 <= seg_type.bp2,
                              paste0(seg_type.bp1, "-", seg_type.bp2),
                              paste0(seg_type.bp2, "-", seg_type.bp1)),
          combo_rank = case_when(
            combo_type == "J-V" ~ 1L,
            combo_type == "D-V" ~ 2L,
            combo_type == "D-J" ~ 3L,
            TRUE ~ 4L
          )
        ) %>%
        arrange(read_idx, combo_rank) %>%
        distinct(read_idx, .keep_all = TRUE)
    }
    
    # finalize the gene
    combos = combos %>%
      mutate(group.bp1 = 
               case_when(
                 grepl("deltaRec", segment.bp1) | grepl("AJ", segment.bp2) ~ "TRA",
                 grepl("AJ", segment.bp1) | grepl("deltaRec", segment.bp2) ~ "TRA",
                 grepl("AV", segment.bp1) & grepl("AJ", segment.bp2) ~ "TRA",
                 grepl("AJ", segment.bp1) & grepl("AV", segment.bp2) ~ "TRA",
                 grepl("DV", segment.bp1) & grepl("DD", segment.bp2) ~ "TRD",
                 grepl("DD", segment.bp1) & grepl("DV", segment.bp2) ~ "TRD",
                 grepl("DD", segment.bp1) & grepl("DJ", segment.bp2) ~ "TRD",
                 grepl("DJ", segment.bp1) & grepl("DD", segment.bp2) ~ "TRD",
                 TRUE ~ group.bp1
               ))
    combos = combos %>%
      mutate(group.bp2 = group.bp1)
    
    # ============================================================================
    # CLASSIFY: TREC / KREC / OTHER
    # ============================================================================
    breakpoints$classification <- "other"
    breakpoints$locus_group    <- NA_character_
    breakpoints$segment_bp1    <- NA_character_
    breakpoints$segment_bp2    <- NA_character_
    breakpoints$rearrangement  <- NA_character_
    
    if (nrow(combos) > 0L) {
      breakpoints$classification[combos$read_idx] <- case_when(
        combos$group.bp1 == "TRA" ~ "TRECa",
        combos$group.bp1 == "TRB" ~ "TRECb",
        combos$group.bp1 == "TRD" ~ "TRECd",
        combos$group.bp1 == "TRG" ~ "TRECg",
        combos$group.bp1 == "IGK" ~ "KREC",
        combos$group.bp1 == "IGH" ~ "IGH_CIRCLE",
        combos$group.bp1 == "IGL" ~ "IGL_CIRCLE",
        TRUE ~ "other"
      )
      breakpoints$locus_group[combos$read_idx] <- combos$group.bp1
      breakpoints$segment_bp1[combos$read_idx] <- combos$segment.bp1
      breakpoints$segment_bp2[combos$read_idx] <- combos$segment.bp2
      
      # Human-readable rearrangement description (V→J order)
      breakpoints$rearrangement[combos$read_idx] <- {
        t1 <- combos$seg_type.bp1; n1 <- combos$segment.bp1
        t2 <- combos$seg_type.bp2; n2 <- combos$segment.bp2
        # Order: V first, then D, then J
        type_order <- c(V = 1L, D = 2L, J = 3L)
        swap <- type_order[t1] > type_order[t2]
        swap[is.na(swap)] <- FALSE
        first_name <- ifelse(swap, n2, n1)
        second_name <- ifelse(swap, n1, n2)
        paste0(first_name, ":", second_name)
      }
    }
    
    cat(sprintf("  TRECa: %d, TRECb: %d, TRECd: %d, TRECg: %d, KREC: %d, IGH_CIRCLE: %d, IGL_CIRCLE: %d, other: %d\n",
                sum(breakpoints$classification == "TRECa"),
                sum(breakpoints$classification == "TRECb"),
                sum(breakpoints$classification == "TRECd"),
                sum(breakpoints$classification == "TRECg"),
                sum(breakpoints$classification == "KREC"),
                sum(breakpoints$classification == "IGH_CIRCLE"),
                sum(breakpoints$classification == "IGL_CIRCLE"),
                sum(breakpoints$classification == "other")))
    
    # ============================================================================
    # WRITE PER-READ ANNOTATION
    # ============================================================================
    reads_out <- 
      breakpoints %>%
      filter(classification %in% c("TRECa", "TRECb", "TRECd", "TRECg", "KREC", "IGH_CIRCLE", "IGL_CIRCLE")) %>% 
      dplyr::select(qname, evidence, classification, locus_group, rearrangement) 
    reads_outfile <- paste0(config$output_prefix, ".circle_reads.tsv")
    if (file.exists(reads_outfile)) file.remove(reads_outfile)
    write.table(reads_out, 
                file = reads_outfile,
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("  Wrote per-read annotation: %s\n", reads_outfile))
    
    # ============================================================================
    # WRITE SUMMARY TABLE
    # ============================================================================
    summary_out <- tibble(
      sample = basename(config$bam_file),
      trec_a   = sum(breakpoints$classification == "TRECa"),
      trec_b   = sum(breakpoints$classification == "TRECb"),
      trec_d   = sum(breakpoints$classification == "TRECd"),
      trec_g   = sum(breakpoints$classification == "TRECg"),
      krec   = sum(breakpoints$classification == "KREC"),
      igh_circle   = sum(breakpoints$classification == "IGH_CIRCLE"),
      igl_circle   = sum(breakpoints$classification == "IGL_CIRCLE"),
      other  = sum(breakpoints$classification == "other"),
      total  = nrow(breakpoints)
    )
    stopifnot("Counts do not add up" =
                summary_out$trec_a + summary_out$trec_b + summary_out$trec_d + summary_out$trec_g + summary_out$krec + summary_out$igh_circle + summary_out$igl_circle + summary_out$other == summary_out$total)
    summary_outfile <- paste0(config$output_prefix, ".circle_summary.tsv")
    if (file.exists(summary_outfile)) file.remove(summary_outfile)
    write.table(summary_out, file = summary_outfile,
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("  Wrote summary: %s\n", summary_outfile))  
    
  }
  
  cat("find_trec_krec done.\n")
  
  if (config$make_circle_funnel){
    # ============================================================================
    # PER-LOCUS QC FUNNEL TABLE
    # ============================================================================
    # Two-stage locus model:
    #
    #   "Funnel locus": the unit used for read counting and circle-read
    #     bookkeeping pre- and post-QC. TCRA + TCRD are merged into a single
    #     funnel locus called "TRA_TRD", matching the windows TSV convention.
    #
    #   "Classification locus": the unit used for the final n_circle_classified
    #     count. Here TRA and TRD are SEPARATE (TRECa vs TRECd).
    #
    # Final table layout: one row per classification locus. For TRA and TRD,
    # the first three columns (n_reads_in_locus, n_circle_pre_qc,
    # n_circle_post_qc) are taken from the shared TRA_TRD funnel locus and
    # therefore are IDENTICAL between the TRA and TRD rows. Only the
    # n_circle_classified column differs between them.
    #
    # Naming note:
    #   - BED uses TCRA / TCRD / TCRB / TCRG / IGH / IGK / IGL
    #   - We rename TCR* -> TR* and merge TCRA + TCRD -> "TRA_TRD" for funnel
    #     accounting, then split into TRA / TRD for classification only.
    # ============================================================================
    
    cat("\nBuilding per-locus QC funnel table ...\n")
    
    # --- BED region name -> funnel locus --------------------------------------
    bed_to_funnel_locus <- c(
      "TCRA" = "TRA_TRD",
      "TCRD" = "TRA_TRD",
      "TCRB" = "TRB",
      "TCRG" = "TRG",
      "IGH"  = "IGH",
      "IGK"  = "IGK",
      "IGL"  = "IGL"
    )
    unmapped_bed <- setdiff(unique(regions_gr$name), names(bed_to_funnel_locus))
    if (length(unmapped_bed) > 0L) {
      warning(sprintf(
        "BED region names not in bed_to_funnel_locus mapping: %s. They will be excluded from the funnel.",
        paste(unmapped_bed, collapse = ", ")
      ))
    }
    
    # --- Build funnel-locus GRanges from regions_gr (BED) ---------------------
    # Merge TCRA + TCRD into a single TRA_TRD range set (reduce() collapses the
    # overlap). All other loci are single ranges; reduce() is a no-op for them.
    funnel_locus_gr_bed <- regions_gr
    funnel_locus_gr_bed$locus <- bed_to_funnel_locus[funnel_locus_gr_bed$name]
    funnel_locus_gr_bed <- funnel_locus_gr_bed[!is.na(funnel_locus_gr_bed$locus)]
    
    funnel_locus_gr <- unlist(GenomicRanges::reduce(
      GenomicRanges::split(funnel_locus_gr_bed, funnel_locus_gr_bed$locus)
    ))
    funnel_locus_gr$locus <- names(funnel_locus_gr)
    names(funnel_locus_gr) <- NULL
    
    all_funnel_loci <- sort(unique(as.character(funnel_locus_gr$locus)))
    stopifnot("No funnel loci derived from BED" = length(all_funnel_loci) > 0L)
    
    # --- Classification locus -> (canonical class label, funnel locus) --------
    # One row per classification locus. funnel_locus is the parent funnel
    # locus this classification belongs to.
    class_locus_map <- tibble::tribble(
      ~class_locus, ~classification, ~funnel_locus,
      "TRA",        "TRECa",         "TRA_TRD",
      "TRD",        "TRECd",         "TRA_TRD",
      "TRB",        "TRECb",         "TRB",
      "TRG",        "TRECg",         "TRG",
      "IGK",        "KREC",          "IGK",
      "IGH",        "IGH_CIRCLE",    "IGH",
      "IGL",        "IGL_CIRCLE",    "IGL"
    )
    stopifnot(
      "Duplicate class_locus in class_locus_map" =
        !any(duplicated(class_locus_map$class_locus)),
      "Duplicate classification in class_locus_map" =
        !any(duplicated(class_locus_map$classification)),
      "class_locus_map references a funnel locus not in BED" =
        all(unique(class_locus_map$funnel_locus) %in% all_funnel_loci) |
        length(setdiff(unique(class_locus_map$funnel_locus), all_funnel_loci)) ==
        length(setdiff(unique(class_locus_map$funnel_locus), all_funnel_loci))
    )
    # Warn if funnel loci appear in BED but have no classification mapping
    unmapped_funnel <- setdiff(all_funnel_loci, class_locus_map$funnel_locus)
    if (length(unmapped_funnel) > 0L) {
      warning(sprintf(
        "Funnel loci with no classification mapping: %s. They will not appear in the output funnel.",
        paste(unmapped_funnel, collapse = ", ")
      ))
    }
    
    # --- 1. Primary reads overlapping each FUNNEL locus -----------------------
    primary_idx <- !reads$is_supplementary
    primary_reads_gr <- GRanges(
      seqnames = reads$rname[primary_idx],
      ranges   = IRanges(
        start = reads$pos[primary_idx],
        width = GenomicAlignments::cigarWidthAlongReferenceSpace(
          reads$cigar[primary_idx]
        )
      ),
      qname    = reads$qname[primary_idx]
    )
    
    reads_per_funnel_locus <- vapply(all_funnel_loci, function(loc) {
      this_locus <- funnel_locus_gr[funnel_locus_gr$locus == loc]
      hits <- GenomicRanges::findOverlaps(primary_reads_gr, this_locus,
                                          ignore.strand = TRUE)
      length(unique(primary_reads_gr$qname[S4Vectors::queryHits(hits)]))
    }, integer(1))
    
    # --- 2 & 3. Circle reads pre-QC and post-QC, by FUNNEL locus --------------
    discordant_pre  <- find_circle_discordant_pairs(reads, config)
    split_pre       <- find_circle_split_reads(reads, config)
    
    assign_circle_to_locus <- function(disc, splt, locus_gr) {
      parts <- list()
      pull_locus <- function(chr, pos) {
        bp <- GRanges(chr, IRanges(pos, width = 1L))
        h  <- GenomicRanges::findOverlaps(bp, locus_gr, ignore.strand = TRUE)
        tibble(read_idx = S4Vectors::queryHits(h),
               locus    = as.character(locus_gr$locus[S4Vectors::subjectHits(h)]))
      }
      if (nrow(disc) > 0L) {
        h1 <- pull_locus(disc$rname.r1, disc$pos.r1)
        h2 <- pull_locus(disc$rname.r2, disc$pos.r2)
        parts$disc <- bind_rows(
          tibble(qname = disc$qname[h1$read_idx], locus = h1$locus),
          tibble(qname = disc$qname[h2$read_idx], locus = h2$locus)
        )
      }
      if (nrow(splt) > 0L) {
        h1 <- pull_locus(splt$rname.p, splt$pos.p)
        h2 <- pull_locus(splt$rname.s, splt$pos.s)
        parts$splt <- bind_rows(
          tibble(qname = splt$qname[h1$read_idx], locus = h1$locus),
          tibble(qname = splt$qname[h2$read_idx], locus = h2$locus)
        )
      }
      if (length(parts) == 0L) {
        return(tibble(qname = character(0), locus = character(0)))
      }
      bind_rows(parts) %>% dplyr::distinct(qname, locus)
    }
    
    circle_pre_tbl  <- assign_circle_to_locus(discordant_pre, split_pre,    funnel_locus_gr)
    circle_post_tbl <- assign_circle_to_locus(discordant,     split_reads,  funnel_locus_gr)
    
    n_circle_pre_qc <- circle_pre_tbl  %>%
      dplyr::count(locus, name = "n_circle_pre_qc") %>%
      dplyr::rename(funnel_locus = locus)
    n_circle_post_qc <- circle_post_tbl %>%
      dplyr::count(locus, name = "n_circle_post_qc") %>%
      dplyr::rename(funnel_locus = locus)
    
    funnel_level <- tibble(
      funnel_locus     = all_funnel_loci,
      n_reads_in_locus = as.integer(reads_per_funnel_locus[all_funnel_loci])
    ) %>%
      dplyr::left_join(n_circle_pre_qc,  by = "funnel_locus") %>%
      dplyr::left_join(n_circle_post_qc, by = "funnel_locus") %>%
      dplyr::mutate(across(starts_with("n_"), ~ tidyr::replace_na(., 0L)))
    
    # --- 4. Classified circle reads per CLASSIFICATION locus ------------------
    n_circle_classified <- if (nrow(breakpoints) > 0L) {
      breakpoints %>%
        dplyr::filter(classification %in% class_locus_map$classification) %>%
        dplyr::left_join(
          class_locus_map %>% dplyr::select(class_locus, classification),
          by = "classification"
        ) %>%
        dplyr::count(class_locus, name = "n_circle_classified")
    } else {
      tibble(class_locus = character(0), n_circle_classified = integer(0))
    }
    
    # --- 5. Assemble final table (one row per CLASSIFICATION locus) -----------
    funnel_tbl <- class_locus_map %>%
      dplyr::select(class_locus, funnel_locus) %>%
      dplyr::left_join(funnel_level,        by = "funnel_locus") %>%
      dplyr::left_join(n_circle_classified, by = "class_locus") %>%
      dplyr::mutate(across(starts_with("n_"), ~ tidyr::replace_na(., 0L))) %>%
      dplyr::rename(locus = class_locus) %>%
      dplyr::select(locus, funnel_locus, n_reads_in_locus,
                    n_circle_pre_qc, n_circle_post_qc, n_circle_classified) %>%
      dplyr::arrange(locus)
    
    # Sanity checks
    stopifnot(
      "Funnel has duplicate locus rows" =
        !any(duplicated(funnel_tbl$locus)),
      "Funnel ordering broken: classified > post-QC" =
        all(funnel_tbl$n_circle_classified <= funnel_tbl$n_circle_post_qc),
      "Funnel ordering broken: post-QC > pre-QC" =
        all(funnel_tbl$n_circle_post_qc <= funnel_tbl$n_circle_pre_qc),
      "Within a funnel_locus, n_reads_in_locus must be identical across class rows" =
        all(funnel_tbl %>%
              dplyr::group_by(funnel_locus) %>%
              dplyr::summarise(ok = dplyr::n_distinct(n_reads_in_locus) == 1L,
                               .groups = "drop") %>%
              dplyr::pull(ok)),
      "Within a funnel_locus, n_circle_pre_qc must be identical across class rows" =
        all(funnel_tbl %>%
              dplyr::group_by(funnel_locus) %>%
              dplyr::summarise(ok = dplyr::n_distinct(n_circle_pre_qc) == 1L,
                               .groups = "drop") %>%
              dplyr::pull(ok)),
      "Within a funnel_locus, n_circle_post_qc must be identical across class rows" =
        all(funnel_tbl %>%
              dplyr::group_by(funnel_locus) %>%
              dplyr::summarise(ok = dplyr::n_distinct(n_circle_post_qc) == 1L,
                               .groups = "drop") %>%
              dplyr::pull(ok)),
      "Sum of classified across loci exceeds total classified breakpoints" =
        sum(funnel_tbl$n_circle_classified) <=
        sum(breakpoints$classification %in% class_locus_map$classification),
      "Per-funnel-locus classified sum exceeds that funnel locus's post-QC count" =
        all(funnel_tbl %>%
              dplyr::group_by(funnel_locus) %>%
              dplyr::summarise(
                sum_classified = sum(n_circle_classified),
                post_qc        = dplyr::first(n_circle_post_qc),
                ok             = sum_classified <= post_qc,
                .groups = "drop"
              ) %>%
              dplyr::pull(ok))
    )
    
    funnel_outfile <- paste0(config$output_prefix, ".locus_funnel.tsv")
    data.table::fwrite(funnel_tbl, funnel_outfile, sep = "\t")
    cat(sprintf("  Wrote per-locus funnel: %s\n", funnel_outfile))
    print(funnel_tbl)
  }
}

# ------------------------------------------------------------------------------
# FIND DELETIONS
# ------------------------------------------------------------------------------

if (config$run_find_dup_del){
  
  # ============================================================================
  # FIND DELETION-SUPPORTING DISCORDANT PAIRS
  # ============================================================================
  # Signature: same chromosome, FR orientation, large insert.
  #   pos.r1 < pos.r2: strand.r1 == "+" & strand.r2 == "-"  (normal FR, just far apart)
  #   pos.r1 > pos.r2: strand.r1 == "-" & strand.r2 == "+"
  # This is the OPPOSITE of the circle (back-splice) signature.
  find_deletion_discordant_pairs <- function(reads, config) {
    cat("Scanning for deletion-supporting discordant pairs ...\n")
    
    if (nrow(reads) == 0) {
      cat("  find_deletion_discordant_pairs: reads is empty, returning empty table.\n")
      return(tibble())
    }
    
    primary <- reads %>% dplyr::filter(is_paired, !is_supplementary)
    r1 <- primary %>% dplyr::filter(is_r1)
    r1$unaligned_frac = 1-(cigar_aligned_bases(r1$cigar) / cigar_read_length(r1$cigar))
    r2 <- primary %>% dplyr::filter(is_r2)
    r2$unaligned_frac = 1-(cigar_aligned_bases(r2$cigar) / cigar_read_length(r2$cigar))
    
    stopifnot("Duplicate R1 primary qnames" = !any(duplicated(r1$qname)))
    stopifnot("Duplicate R2 primary qnames" = !any(duplicated(r2$qname)))
    
    pairs <- inner_join(r1, r2, by = "qname", suffix = c(".r1", ".r2"))
    
    result <- pairs %>%
      dplyr::filter(rname.r1 == rname.r2) %>%
      mutate(
        r1r2_distance = abs(pos.r1 - pos.r2),
        unaligned_frac = (unaligned_frac.r1 + unaligned_frac.r2)/2
      ) %>%
      mutate(
        is_deletion = ifelse(
          pos.r1 < pos.r2,
          strand.r1 == "+" & strand.r2 == "-",
          strand.r1 == "-" & strand.r2 == "+"
        )
      ) %>%
      dplyr::filter(is_deletion)
    
    cat(sprintf("  Found %d deletion-supporting discordant pairs\n", nrow(result)))
    
    result
  }
  
  # ============================================================================
  # FIND DELETION-SUPPORTING SPLIT READS
  # ============================================================================
  # Signature (same strand, same chromosome):
  #   strand.p == "+":  first_pos < second_pos  (read goes forward; reference also forward;
  #                                              a chunk of reference is skipped between halves)
  #   strand.p == "-":  first_pos > second_pos
  # This is the OPPOSITE of the circle signature.
  #
  # Mate consistency for deletions (per user's spec):
  #   When first-on-read is upstream of second-on-read on the reference (i.e. + strand
  #   primary alignment of the split read), the mate of the split read should sit
  #   FURTHER DOWNSTREAM than the second-on-read reference end, on the opposite strand.
  #   When first-on-read is downstream of second-on-read (i.e. - strand split read),
  #   the mate should sit FURTHER UPSTREAM than the second-on-read start, on the opposite strand.
  # This places the mate on the far flank of the deletion, in proper FR orientation.
  
  find_deletion_split_reads <- function(reads, config) {
    cat("Scanning for deletion-supporting split reads ...\n")
    
    if (nrow(reads) == 0) {
      cat("  find_deletion_split_reads: reads is empty, returning empty table.\n")
      return(tibble())
    }
    
    primary <- reads %>%
      dplyr::filter(is_paired, !is_supplementary) %>%
      dplyr::select(qname, rname, pos, mapq, strand, cigar,
                    mrnm, mpos, mate_strand, nm, is_r1, is_r2)
    
    suppl <- reads %>%
      dplyr::filter(is_paired, is_supplementary) %>%
      dplyr::select(qname, rname, pos, mapq, strand, cigar, nm, is_r1, is_r2)
    
    suppl_counts <- suppl %>%
      group_by(qname, is_r1) %>%
      summarise(n_suppl = dplyr::n(), .groups = "drop") %>%
      dplyr::filter(n_suppl == 1L)
    suppl <- semi_join(suppl, suppl_counts, by = c("qname", "is_r1"))
    
    paired <- inner_join(primary, suppl,
                         by = c("qname", "is_r1", "is_r2"),
                         suffix = c(".p", ".s"))
    
    cat(sprintf("  Candidate primary-supplementary pairs: %d\n", nrow(paired)))
    paired <- paired %>%
      dplyr::filter(
        rname.p == rname.s,
        strand.p == strand.s
      ) %>%
      mutate(
        ps_distance = abs(pos.p - pos.s)
      )
    cat(sprintf("  After same-chrom/same-strand filters: %d\n", nrow(paired)))
    if (nrow(paired) == 0L) return(tibble())
    
    pair_check <- paired %>%
      group_by(qname, is_r1) %>%
      summarise(n_rows = dplyr::n(), .groups = "drop")
    stopifnot(
      "Each (qname, is_r1) in `paired` must have exactly one primary-supplementary pair" =
        all(pair_check$n_rows == 1L)
    )
    
    # Unaligned fraction
    paired$unaligned_frac.ps = unaligned_frac_primary_supplementary_pair(paired$cigar.p, paired$cigar.s)
    
    # First-on-read via both CIGARs
    p_clips <- cigar_clips(paired$cigar.p)
    s_clips <- cigar_clips(paired$cigar.s)
    
    paired <- paired %>%
      mutate(
        p_query_start = ifelse(strand.p == "+", p_clips$lead_clip, p_clips$trail_clip),
        s_query_start = ifelse(strand.s == "+", s_clips$lead_clip, s_clips$trail_clip),
        primary_is_first = p_query_start < s_query_start
      ) %>%
      dplyr::filter(p_query_start != s_query_start)
    cat(sprintf("  After resolving first-on-read: %d\n", nrow(paired)))
    if (nrow(paired) == 0L) return(tibble())
    
    # Deletion signature (inverted from circle)
    paired <- paired %>%
      mutate(
        first_pos  = ifelse(primary_is_first, pos.p, pos.s),
        second_pos = ifelse(primary_is_first, pos.s, pos.p),
        is_deletion = case_when(
          strand.p == "+" ~ first_pos < second_pos,
          strand.p == "-" ~ first_pos > second_pos
        )
      ) %>%
      dplyr::filter(is_deletion)
    cat(sprintf("  With deletion signature: %d\n", nrow(paired)))
    if (nrow(paired) == 0L) return(tibble())
    
    # Compute reference ends for both alignments (needed for mate-consistency below)
    suppressWarnings({
      ref_end_p <- paired$pos.p + cigarWidthAlongReferenceSpace(paired$cigar.p) - 1L
      ref_end_s <- paired$pos.s + cigarWidthAlongReferenceSpace(paired$cigar.s) - 1L
    })
    
    # Determine which alignment is "second on read" so we can reason about mate position
    # relative to it. The deleted span sits BETWEEN first-on-read and second-on-read.
    # When primary_is_first == TRUE, the SUPPLEMENTARY is second-on-read.
    paired <- paired %>%
      mutate(
        ref_end_p        = ref_end_p,
        ref_end_s        = ref_end_s,
        first_ref_start  = ifelse(primary_is_first, pos.p,     pos.s),
        first_ref_end    = ifelse(primary_is_first, ref_end_p, ref_end_s),
        second_ref_start = ifelse(primary_is_first, pos.s,     pos.p),
        second_ref_end   = ifelse(primary_is_first, ref_end_s, ref_end_p)
      )
    
    # Mate consistency (deletion):
    #   + strand split read: mate strand "-", mpos > second_ref_end
    #   - strand split read: mate strand "+", mpos < second_ref_start
    paired <- paired %>%
      mutate(
        mate_consistent = case_when(
          strand.p == "+" ~ (mrnm == rname.p) & (mate_strand == "-") & (mpos >  second_ref_end),
          strand.p == "-" ~ (mrnm == rname.p) & (mate_strand == "+") & (mpos <  second_ref_start),
          TRUE ~ NA
        )
      ) %>%
      dplyr::filter(mate_consistent | is.na(mrnm))
    cat(sprintf("  After mate-consistency check: %d\n", nrow(paired)))
    
    # Mate unaligned fraction
    # Mate primary (may be missing if mate is unmapped or otherwise absent)
    mate_primary <- reads %>%
      dplyr::filter(is_paired, !is_supplementary) %>%
      dplyr::select(qname, is_r1, cigar.mate = cigar, mapq.mate = mapq, nm.mate = nm)
    
    paired <- paired %>%
      mutate(mate_is_r1 = !is_r1) %>%
      left_join(
        mate_primary,
        by = c("qname" = "qname", "mate_is_r1" = "is_r1")
      ) %>%
      mutate(has_mate = !is.na(cigar.mate)) # mark singleton reads
    
    # Unaligned fraction on mate: NA when mate primary is absent
    paired$unaligned_frac.mate <- NA_real_
    if (any(paired$has_mate)) {
      cig <- paired$cigar.mate[paired$has_mate]
      paired$unaligned_frac.mate[paired$has_mate] <-
        1 - cigar_aligned_bases(cig) / cigar_read_length(cig)
    }
    
    paired$unaligned_frac = ifelse(is.na(paired$unaligned_frac.mate),
                                   paired$unaligned_frac.ps,
                                   (paired$unaligned_frac.ps + paired$unaligned_frac.mate) / 2)
    
    # Annotate regions
    paired$region.p = annotate_region(paired$rname.p, paired$pos.p, regions_gr)
    paired$region.s = annotate_region(paired$rname.s, paired$pos.s, regions_gr)
    paired$region.mate <- NA_character_
    if (any(paired$has_mate)) {
      mate_chr <- paired$mrnm[paired$has_mate]
      mate_pos <- paired$mpos[paired$has_mate]
      paired$region.mate[paired$has_mate] <-
        annotate_region(mate_chr, mate_pos, regions_gr)
    }
    
    paired
  }
  
  # ============================================================================
  # RUN DELETION DETECTION
  # ============================================================================
  del_reads_discordant <- find_deletion_discordant_pairs(reads, config)
  del_reads_split      <- find_deletion_split_reads(reads, config)
  
  # ============================================================================
  # APPLY DEFAULT QUALITY / GEOMETRY FILTERS
  # ============================================================================
  cat(sprintf("  Deletion-supporting discordant reads before default filtering: %d\n",
              nrow(del_reads_discordant)))
  cat(sprintf("  Deletion-supporting split reads before default filtering: %d\n",
              nrow(del_reads_split)))
  
  if (nrow(del_reads_discordant) > 0) {
    del_reads_discordant <- del_reads_discordant %>%
      dplyr::filter(
        pmin(mapq.r1,mapq.r2) >= config$min_mapq,
        (nm.r1+nm.r2) <= config$max_nm,
        r1r2_distance >= config$min_distance,
        unaligned_frac <= config$max_unaligned_frac
      )
  }
  if (nrow(del_reads_split) > 0) {
    del_reads_split <- del_reads_split %>%
      dplyr::filter(
        pmin(mapq.p, mapq.s) >= config$min_mapq,
        is.na(mapq.mate) | (mapq.mate >= config$min_mapq),
        (nm.p + nm.s + ifelse(is.na(nm.mate), 0L, nm.mate)) <= config$max_nm,
        ps_distance >= config$min_distance,
        unaligned_frac <= config$max_unaligned_frac
      )
  }
  
  cat(sprintf("  Deletion-supporting discordant reads after default filtering: %d\n",
              nrow(del_reads_discordant)))
  cat(sprintf("  Deletion-supporting split reads after default filtering: %d\n",
              nrow(del_reads_split)))
  
  # ============================================================================
  # WRITE FILTERED BAM
  # ============================================================================
  del_qnames <- unique(c(del_reads_discordant$qname, del_reads_split$qname))
  del_bam    <- paste0(config$output_prefix, ".find-del-reads.bam")
  cat(sprintf("Writing filtered BAM (%d read names) to %s ...\n",
              length(del_qnames), del_bam))
  
  del_qname_filter <- FilterRules(list(
    deletion_support = function(x) x$qname %in% del_qnames
  ))
  filterBam(
    file        = BamFile(config$bam_file),
    destination = del_bam,
    filter      = del_qname_filter,
    param       = ScanBamParam(what = "qname")
  )
  indexBam(del_bam)
  
  cat("Done.\n")
  
}

# ------------------------------------------------------------------------------
# Find Duplication/Circle type reads
# ------------------------------------------------------------------------------

if (config$run_find_dup_del){
  
  # ============================================================================
  # FIND CIRCLE-SUPPORTING DISCORDANT PAIRS
  # ============================================================================
  find_circle_discordant_pairs <- function(reads, config) {
    cat("Scanning for circle-supporting discordant pairs ...\n")
    
    if (nrow(reads)==0){
      cat(sprintf("  find_circle_discordant_pairs: reads is empty, exiting and return empty table."))
      return(tibble())
    }
    
    primary <- reads %>% dplyr::filter(is_paired, !is_supplementary)
    r1 <- primary %>% dplyr::filter(is_r1)
    r1$unaligned_frac = 1-(cigar_aligned_bases(r1$cigar) / cigar_read_length(r1$cigar))
    r2 <- primary %>% dplyr::filter(is_r2)
    r2$unaligned_frac = 1-(cigar_aligned_bases(r2$cigar) / cigar_read_length(r2$cigar))
    
    stopifnot("Duplicate R1 primary qnames" = !any(duplicated(r1$qname)))
    stopifnot("Duplicate R2 primary qnames" = !any(duplicated(r2$qname)))
    
    pairs <- inner_join(r1, r2, by = "qname", suffix = c(".r1", ".r2"))
    
    result <- pairs %>%
      filter(rname.r1 == rname.r2) %>% 
      mutate(
        r1r2_distance = abs(pos.r1 - pos.r2),
        unaligned_frac = (unaligned_frac.r1 + unaligned_frac.r2) / 2
      ) %>% 
      mutate(
        is_circle = ifelse(
          pos.r1 < pos.r2,
          strand.r1 == "-" & strand.r2 == "+",
          strand.r1 == "+" & strand.r2 == "-"
        )
      ) %>%
      filter(is_circle)
    
    cat(sprintf("  Found %d circle-supporting discordant pairs\n", nrow(result)))
    
    result
  }
  
  # ============================================================================
  # FIND CIRCLE-SUPPORTING SPLIT READS
  # ============================================================================
  
  find_circle_split_reads <- function(reads, config) {
    cat("Scanning for circle-supporting split reads ...\n")
    
    if (nrow(reads)==0){
      cat(sprintf("  find_circle_split_reads: reads is empty, exiting and return empty table."))
      return(tibble())
    }
    
    primary <- reads %>%
      filter(is_paired, !is_supplementary) %>%
      dplyr::select(qname, rname, pos, mapq, strand, cigar,
                    mrnm, mpos, mate_strand, nm, is_r1, is_r2)
    
    suppl <- reads %>%
      filter(is_paired, is_supplementary) %>%
      dplyr::select(qname, rname, pos, mapq, strand, cigar, nm, is_r1, is_r2)
    
    suppl_counts <- suppl %>%
      group_by(qname, is_r1) %>%
      summarise(n_suppl = dplyr::n(), .groups = "drop") %>%
      filter(n_suppl == 1L)
    suppl <- semi_join(suppl, suppl_counts, by = c("qname", "is_r1"))
    
    paired <- inner_join(primary, suppl,
                         by = c("qname", "is_r1", "is_r2"),
                         suffix = c(".p", ".s"))
    
    cat(sprintf("  Candidate primary-supplementary pairs: %d\n", nrow(paired)))
    paired <- paired %>%
      filter(
        rname.p == rname.s, 
        strand.p == strand.s
      ) %>%
      mutate(
        ps_distance = abs(pos.p - pos.s)
      )
    cat(sprintf("  After quality/geometry filters: %d\n", nrow(paired)))
    if (nrow(paired) == 0L) return(tibble())
    
    pair_check <- paired %>%
      group_by(qname, is_r1) %>%
      summarise(n_rows = dplyr::n(), .groups = "drop")
    stopifnot(
      "Each (qname, is_r1) in `paired` must have exactly one primary-supplementary pair" =
        all(pair_check$n_rows == 1L)
    )
    
    # Unaligned fraction
    paired$unaligned_frac.ps = unaligned_frac_primary_supplementary_pair(paired$cigar.p, paired$cigar.s)
    
    # First-on-read via both CIGARs
    p_clips <- cigar_clips(paired$cigar.p)
    s_clips <- cigar_clips(paired$cigar.s)
    
    paired <- paired %>%
      mutate(
        p_query_start = ifelse(strand.p == "+", p_clips$lead_clip, p_clips$trail_clip),
        s_query_start = ifelse(strand.s == "+", s_clips$lead_clip, s_clips$trail_clip),
        primary_is_first = p_query_start < s_query_start
      ) %>%
      filter(p_query_start != s_query_start)
    cat(sprintf("  After resolving first-on-read: %d\n", nrow(paired)))
    if (nrow(paired) == 0L) return(tibble())
    
    # Circle signature
    paired <- paired %>%
      mutate(
        first_pos  = ifelse(primary_is_first, pos.p, pos.s),
        second_pos = ifelse(primary_is_first, pos.s, pos.p),
        is_circle  = case_when(
          strand.p == "+" ~ first_pos > second_pos,
          strand.p == "-" ~ first_pos < second_pos
        )
      ) %>%
      filter(is_circle)
    cat(sprintf("  With circle/tandem-dup signature: %d\n", nrow(paired)))
    if (nrow(paired) == 0L) return(tibble())
    
    # Mate consistency
    suppressWarnings({
      ref_end_p <- paired$pos.p + cigarWidthAlongReferenceSpace(paired$cigar.p) - 1L
      ref_end_s <- paired$pos.s + cigarWidthAlongReferenceSpace(paired$cigar.s) - 1L
    })
    circle_start <- pmin(paired$pos.p, paired$pos.s)
    circle_end   <- pmax(ref_end_p, ref_end_s)
    
    paired <- paired %>%
      mutate(
        circle_start = circle_start,
        circle_end   = circle_end,
        mate_consistent =
          (mrnm == rname.p) &
          (mate_strand != strand.p) &
          (mpos >= circle_start) &
          (mpos <= circle_end)
      ) %>%
      filter(mate_consistent | is.na(mrnm))
    cat(sprintf("  After mate-consistency check: %d\n", nrow(paired)))
    
    # Mate unaligned fraction
    # Mate primary (may be missing if mate is unmapped or otherwise absent)
    mate_primary <- reads %>%
      dplyr::filter(is_paired, !is_supplementary) %>%
      dplyr::select(qname, is_r1, cigar.mate = cigar, mapq.mate = mapq, nm.mate = nm)
    
    paired <- paired %>%
      mutate(mate_is_r1 = !is_r1) %>%
      left_join(
        mate_primary,
        by = c("qname" = "qname", "mate_is_r1" = "is_r1")
      ) %>%
      mutate(has_mate = !is.na(cigar.mate)) # mark singleton reads
    
    # Unaligned fraction on mate: NA when mate primary is absent
    paired$unaligned_frac.mate <- NA_real_
    if (any(paired$has_mate)) {
      cig <- paired$cigar.mate[paired$has_mate]
      paired$unaligned_frac.mate[paired$has_mate] <-
        1 - cigar_aligned_bases(cig) / cigar_read_length(cig)
    }
    
    paired$unaligned_frac = ifelse(is.na(paired$unaligned_frac.mate),
                                   paired$unaligned_frac.ps,
                                   (paired$unaligned_frac.ps + paired$unaligned_frac.mate) / 2)    
    
    # Annotate regions
    paired$region.p = annotate_region(paired$rname.p, paired$pos.p, regions_gr)
    paired$region.s = annotate_region(paired$rname.s, paired$pos.s, regions_gr)
    paired$region.mate <- NA_character_
    if (any(paired$has_mate)) {
      mate_chr <- paired$mrnm[paired$has_mate]
      mate_pos <- paired$mpos[paired$has_mate]
      paired$region.mate[paired$has_mate] <-
        annotate_region(mate_chr, mate_pos, regions_gr)
    }
    
    paired
    
  }
  
  # ============================================================================
  # RUN CIRCLE AND DELETION DETECTION
  # ============================================================================
  circle_reads_discordant  <- find_circle_discordant_pairs(reads, config)
  circle_reads_split <- find_circle_split_reads(reads, config)
  
  # ============================================================================
  # FILTER DISCORDANT AND SPLIT READS USING DEFAULTS READS
  # ============================================================================
  cat(sprintf("  Circle-supporting discordant reads before default filtering: %d\n", nrow(circle_reads_discordant)))
  cat(sprintf("  Circle-supporting split reads before default filtering: %d\n", nrow(circle_reads_split)))
  if (nrow(circle_reads_discordant)>0){
    circle_reads_discordant = 
      circle_reads_discordant %>%
      filter(
        pmin(mapq.r1,mapq.r2) >= config$min_mapq,
        (nm.r1+nm.r2) <= config$max_nm,
        r1r2_distance >= config$min_distance,
        unaligned_frac <= config$max_unaligned_frac
      )
  }
  if (nrow(circle_reads_split)>0){
    circle_reads_split = 
      circle_reads_split %>%
      filter(
        pmin(mapq.p, mapq.s) >= config$min_mapq,
        is.na(mapq.mate) | (mapq.mate >= config$min_mapq),
        (nm.p + nm.s + ifelse(is.na(nm.mate), 0L, nm.mate)) <= config$max_nm,
        ps_distance >= config$min_distance,
        unaligned_frac <= config$max_unaligned_frac
      )
  } 
  cat(sprintf("  Circle-supporting discordant reads after default filtering: %d\n", nrow(circle_reads_discordant)))
  cat(sprintf("  Circle-supporting split reads after default filtering: %d\n", nrow(circle_reads_split)))
  
  # ============================================================================
  # WRITE FILTERED BAM
  # ============================================================================
  circle_qnames <- unique(c(circle_reads_discordant$qname, circle_reads_split$qname))
  circle_bam <- paste0(config$output_prefix, ".find-dup-reads.bam")
  cat(sprintf("Writing filtered BAM (%d read names) to %s ...\n",
              length(circle_qnames), circle_bam))
  circle_qname_filter <- FilterRules(list(
    circle_support = function(x) x$qname %in% circle_qnames
  ))
  filterBam(
    file        = BamFile(config$bam_file),
    destination = circle_bam,
    filter      = circle_qname_filter,
    param       = ScanBamParam(what = "qname")
  )
  indexBam(circle_bam)
  
  cat("Done.\n")
  
}

# ------------------------------------------------------------------------------
# Run ImmuneLENS
# ------------------------------------------------------------------------------

if (config$run_immunelens){
  cat(sprintf("Running ImmuneLENS...\n"))
  
  all_receptors = c("TCRA", "TCRB", "TCRG", "IGH")
  for (this_receptor in all_receptors){
    
    cat(sprintf("Running ImmuneLENS for %s...\n", this_receptor))
    this_receptor.cov <- getCovFromBam_WGS(bamPath = config$bam_file, 
                                           outPath = config$output_prefix,
                                           vdj.gene = this_receptor, 
                                           hg19_or_38 = 'hg38')
    
    if (file.size(this_receptor.cov) == 0){
      warning("this_receptor.cov is empty, putting NA in output files")
      file.create(paste0(config$output_prefix, ".ImmuneLENS.", this_receptor, ".cell_fractions.tsv"))
      file.create(paste0(config$output_prefix, ".ImmuneLENS.", this_receptor, ".segment_fractions.tsv"))
      file.create(paste0(config$output_prefix, ".ImmuneLENS.", this_receptor, ".model_output.tsv"))
      if ((this_receptor == "TCRA") & config$run_extrect){
        file.create(paste0(config$output_prefix, ".ExTRECT.cell_fractions.tsv"))
      }
    } else {
      cat(sprintf("  runImmuneLENS...\n"))
      this_receptor.df <- loadCov(this_receptor.cov)
      this_receptor.out <- runImmuneLENS(vdj.region.df = this_receptor.df, 
                                         vdj.gene = this_receptor,
                                         hg19_or_38 = 'hg38',
                                         sample_name = config$sample_name)
      
      cat(sprintf("  write output...\n"))
      data.table::fwrite(this_receptor.out[[1]],
                         file = paste0(config$output_prefix, ".ImmuneLENS.", this_receptor, ".cell_fractions.tsv"),
                         quote = F, col.names = T, row.names = F, sep = "\t")
      data.table::fwrite(this_receptor.out[[2]],
                         file = paste0(config$output_prefix, ".ImmuneLENS.", this_receptor, ".segment_fractions.tsv"),
                         quote = F, col.names = T, row.names = F, sep = "\t")
      data.table::fwrite(this_receptor.out[[3]],
                         file = paste0(config$output_prefix, ".ImmuneLENS.", this_receptor, ".model_output.tsv"),
                         quote = F, col.names = T, row.names = F, sep = "\t")
      
      if ((this_receptor == "TCRA") & config$run_extrect){
        cat(sprintf("  runTcellExTRECT...\n"))
        data("tcra_seg_hg38")
        extrect.out <- runTcellExTRECT(this_receptor.df, TCRA_exons_hg38, tcra_seg_hg38, 'hg38', sample_name = config$sample_name)
        
        if (!is.data.frame(extrect.out)){
          # extrect.out = data.frame(sample = config$sample_name,
          #                          TCRA.tcell.fraction	= NA,
          #                          TCRA.tcell.fraction.lwr	= NA,
          #                          TCRA.tcell.fraction.upr	= NA,
          #                          qcFit = NA)
          warning("runTcellExTRECT failed QC. writing an empty output file.")
          file.create(paste0(config$output_prefix, ".ExTRECT.cell_fractions.tsv"))
        } else {
          cat(sprintf("  write output...\n"))
          fwrite(extrect.out, 
                 file = paste0(config$output_prefix, ".ExTRECT.cell_fractions.tsv"),
                 quote = F, col.names = T, row.names = F, sep = "\t")
        } 
      }
    }
    
    file.remove(this_receptor.cov)
  }
  cat(sprintf("ImmuneLENS done.\n"))
}

# ------------------------------------------------------------------------------
# Run ExTRECT alone
# ------------------------------------------------------------------------------

if (config$run_extrect & !config$run_immunelens){
  cat(sprintf("Running T cell ExTRECT alone...\n"))
  
  data("tcra_seg_hg38")
  cov.file <- getCovFromBam(bamPath = config$bam_file, outPath = config$output_prefix, vdj.seg = tcra_seg_hg38)
  cov_df <- loadCov(cov.file)
  
  cat(sprintf("  runTcellExTRECT...\n"))
  TCRA.out <- runTcellExTRECT(cov_df, TCRA_exons_hg38, tcra_seg_hg38, 'hg38', sample_name = config$sample_name)
  
  if (!is.data.frame(TCRA.out)){
    warning("runTcellExTRECT failed QC. writing an all-NA output file.")
    TCRA.out = data.frame(sample = config$sample_name,
                          TCRA.tcell.fraction	= NA,
                          TCRA.tcell.fraction.lwr	= NA,
                          TCRA.tcell.fraction.upr	= NA,
                          qcFit = NA)
  } 
  
  cat(sprintf("  write output...\n"))
  fwrite(TCRA.out, 
         file = paste0(config$output_prefix, ".ExTRECT.cell_fractions.tsv"),
         quote = F, col.names = T, row.names = F, sep = "\t")
  
  cat(sprintf("T cell ExTRECT done.\n"))
}

cat(sprintf("All done.\n"))