# ---------------------------------------------------------------------------
# assign_sample_id_from_fov.R
#
# Upstream pipeline - helper sourced by merge_seu.R (step 1 of 9)
#
# Repository: Tangle_Neighbours - spatial analysis of the neuronal and glial
# microenvironment around tau-tangle-bearing neurons (CosMx 6k + IMC).
#
# Paths in this file are PLACEHOLDERS (<PROJECT_ROOT>, <RDS_ROOT>, ...). Set them
# to your own copy before running - see PLACEHOLDERS.md.
# Donor identifiers are UK Brain Banks Network (BBN) IDs, not brain-bank case IDs.
# ---------------------------------------------------------------------------
assign_sample_id_from_fov <- function(
    seu,
    fov_col = "set_fov",
    ranges,
    exceptions = NULL,
    default_label = NA_character_,
    verbose = TRUE
) {
  stopifnot(inherits(seu, "Seurat"))
  md <- seu@meta.data
  if (!fov_col %in% colnames(md)) {
    stop(sprintf("Column '%s' not found in metadata.", fov_col))
  }
  if (!all(c("sample_id", "start", "end") %in% colnames(ranges))) {
    stop("ranges must have columns: sample_id, start, end")
  }
  if (any(ranges$start > ranges$end)) stop("Each range must satisfy start <= end")
  
  # Coerce set_fov to integer (tolerant to factor/char)
  fov_vec <- md[[fov_col]]
  if (is.factor(fov_vec)) fov_vec <- as.character(fov_vec)
  suppressWarnings(fov_num <- as.integer(fov_vec))
  
  out <- rep(default_label, length(fov_num))
  
  # Assign by ranges (later rows can override earlier if overlapping)
  for (i in seq_len(nrow(ranges))) {
    idx <- !is.na(fov_num) & fov_num >= ranges$start[i] & fov_num <= ranges$end[i]
    out[idx] <- ranges$sample_id[i]
  }
  
  # Apply exceptions (explicit overrides)
  if (!is.null(exceptions) && nrow(exceptions) > 0) {
    stopifnot(all(c("fov", "sample_id") %in% colnames(exceptions)))
    if (is.factor(exceptions$fov)) exceptions$fov <- as.character(exceptions$fov)
    suppressWarnings(exc_fov <- as.integer(exceptions$fov))
    for (j in seq_len(nrow(exceptions))) {
      idx <- !is.na(fov_num) & fov_num == exc_fov[j]
      out[idx] <- exceptions$sample_id[j]
    }
  }
  
  # QA messages
  if (verbose) {
    unmapped <- sum(is.na(out))
    if (unmapped > 0) {
      warning(sprintf("Unmapped cells: %d (no matching range or NA %s).", unmapped, fov_col))
    }
    # Print quick cross-tab
    print(table(out, useNA = "ifany"))
  }
  
  seu$sample_id <- out
  return(seu)
}