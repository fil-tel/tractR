#' Compute the Tract Frequency Spectrum (TFS)
#'
#' Given a binary matrix the function computes the Tract Frequency Spectrum (TFS).
#'
#' @param bin_mat Binary matrix. Rows correspond to bin, columns to sampled chromosomes/individuals.
#' @param n_sample Integer. Number of sampled chromosomes/individuals. Ideally it should be equivalent to the number of columns, but it is not if some individuals do not carry any tracts.
#'
#' @returns Vector. It corresponds to the TFS.
#' @export
#'
#' @examples
compute_tfs <- function(bin_mat, n_sample=100){
  tfs <- rep(0, n_sample)
  # find how many individual share a tract
  counts <- table(rowSums(bin_mat))
  # remove the 0 counts if there are
  if(!is.na(counts["0"])) counts <- counts[-1]
  tfs[as.numeric(names(counts))] <- counts
  # normalize to get proportion of tracts
  tfs <- tfs/sum(tfs)
  tfs
}
