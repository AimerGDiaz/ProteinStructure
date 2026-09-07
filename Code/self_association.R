# Utilities for summarizing long-range intramolecular contacts.

compute_self_extrema <- function(matrix, residue_range = NULL, seq_sep = 20,
                                 which_func, fill_value) {
  if (!is.matrix(matrix) || nrow(matrix) != ncol(matrix)) {
    stop("`matrix` must be a square matrix")
  }
  if (is.null(residue_range)) residue_range <- seq_len(nrow(matrix))
  mat <- matrix[residue_range, residue_range, drop = FALSE]
  mat[abs(row(mat) - col(mat)) <= seq_sep] <- fill_value

  extrema <- lapply(seq_len(ncol(mat)), function(column) {
    values <- mat[, column]
    if (all(values == fill_value)) return(c(value = NA, partner = NA))
    index <- which_func(values)
    c(value = values[index], partner = residue_range[index])
  })

  values <- do.call(rbind, extrema)
  data.frame(
    Residue = residue_range,
    Value = values[, "value"],
    Partner = values[, "partner"]
  )
}

compute_self_minima <- function(pae_matrix, residue_range = NULL, seq_sep = 20) {
  result <- compute_self_extrema(
    pae_matrix, residue_range, seq_sep, which.min, Inf
  )
  names(result)[2:3] <- c("MinPAE", "MinPAE_Index")
  result
}

compute_self_maxima <- function(contact_matrix, residue_range = NULL,
                                seq_sep = 20) {
  result <- compute_self_extrema(
    contact_matrix, residue_range, seq_sep, which.max, -Inf
  )
  names(result)[2:3] <- c("MaxContactProb", "MaxContactProb_Index")
  result
}

