# List of functions to obtain tract-encoding binary matrices for the different
# tract-encoding approaches.
# Each function has a simulated data version (functionname_sim), and an empirical data one (functionname_emp)..

#' Get binary matrix of tracts for uniqueness approach - SIMULATION
#'
#' @param tracts_gr
#'
#' @returns
#' @export
#'
#' @examples
get_bin_mat_unique_sim <- function(tracts_gr){
  if (!("name" %in% colnames(GenomicRanges::mcols(tracts_gr))))
    stop(
      "There is no metadata column named 'name'. \n Note that the column storing the ID of the individual has to be called 'name'."
    )
  tracts_unique_gr <- unique(tracts_gr)
  if(length(tracts_unique_gr)>1){
    bin_mat <- sapply(unique(tracts_gr$name), function(x) get_exact_match(filter(tracts_gr, name==x), tracts_unique_gr))
  }else{
    bin_mat <- matrix(1, nrow = 1, ncol = length(unique(tracts_gr$name)))
    rownames(bin_mat) <- as.character(tracts_unique_gr)
    colnames(bin_mat) <- unique(tracts_gr$name)
  }
  return(bin_mat)
}

#' Get binary matrix of tracts for uniqueness approach - EMPIRICAL
#'
#' @param tracts_gr
#'
#' @returns
#' @export
#'
#' @examples
get_bin_mat_unique_emp <- function(tracts_gr){
  if(!("name" %in% colnames(GenomicRanges::mcols(tracts_gr))))
    stop(
      "There is no metadata column named 'name'. \n Note that the column storing the ID of the individual has to be called 'name'."

    )
  tracts_unique_gr <- unique(tracts_gr)
  if(length(tracts_unique_gr)>1){
    bin_mat <- sapply(unique(tracts_gr$name), function(x) get_exact_match(filter(tracts_gr, name==x), tracts_unique_gr))
  }else{
    bin_mat <- matrix(1, nrow = 1, ncol = length(unique(tracts_gr$name)))
    rownames(bin_mat) <- as.character(tracts_unique_gr)
    colnames(bin_mat) <- unique(tracts_gr$name)
  }
  return(bin_mat)
}

#' Get binary matrix of tracts for windows approach - SIMULATION
#'
#' @param tracts_gr
#' @param window_size
#' @param step_size
#' @param len_chr
#' @param ncores
#'
#' @returns
#' @export
#'
#' @examples
get_bin_mat_windows_sim <- function(tracts_gr, window_size = 50e3, step_size = 50e3, len_chr = NULL, ncores=1){
  if (is.null(len_chr)) {
    stop("Specify length of chromosome.")
  }
  windows_gr <- generate_windows_sim(window_size = window_size, step_size = step_size, len_chr = len_chr)
  bin_mat <- sapply(unique(tracts_gr$name), function(sample_id){
    compute_bin_vec(tracts_gr, windows_gr, sample_id, ncores=ncores)
  })
  rownames(bin_mat) <- paste0(seqnames(windows_gr), ":",ranges(windows_gr))
  bin_mat <- ifelse(bin_mat>0, 1, 0)
  bin_mat <- bin_mat[rowSums(bin_mat)!=0,]
  return(bin_mat)
}

#' Get binary matrix of tracts for windows approach - EMPIRICAL
#'
#' @param tracts_gr
#' @param window_size
#' @param step_size
#' @param ncores
#'
#' @returns
#' @export
#'
#' @examples
get_bin_mat_windows_emp <- function(tracts_gr, window_size = 50e3, step_size = 50e3, ncores=1){
  windows_gr <- generate_windows_no_gap(window_size = window_size, step_size = step_size)
  bin_mat <- sapply(unique(tracts_gr$name), function(sample_id){
    compute_bin_vec(tracts_gr, windows_gr, sample_id, ncores=ncores)
  })
  rownames(bin_mat) <- paste0(seqnames(windows_gr), ":",ranges(windows_gr))
  bin_mat <- ifelse(bin_mat>0, 1, 0)
  bin_mat <- bin_mat[rowSums(bin_mat)!=0,]
  return(bin_mat)
}

#' Get binary matrix of tracts for subtracts approach - SIMULATION
#'
#' @param tracts_gr
#'
#' @returns
#' @export
#'
#' @examples
get_bin_mat_subtracts_sim <- function(tracts_gr){
  # iterate over each chromosome
  # to get one matrix for each chromosome
  col_names <- unique(tracts_gr$name)
  bin_mat_l <- lapply(as.character(unique(seqnames(tracts_gr))),
                      function(chrom) {
                        chrom_gr <- filter(tracts_gr, seqnames==chrom)
                        mat <- get_bin_mat_subtracts_chrom(chrom_gr)
                        # This is because the some inds might not have tracts
                        # for some chromosomes and in order to create one big matrix it is important
                        # to have the dimensions consistent
                        full_matrix <- matrix(0, nrow = nrow(mat), ncol = length(col_names), dimnames = list(rownames(mat), col_names))
                        full_matrix[rownames(mat), colnames(mat)] <- mat
                        full_matrix
                      })
  bin_mat <- do.call(rbind, bin_mat_l)
  return(bin_mat)
}

#' Get binary matrix of tracts for subtracts approach - EMPIRICAL
#'
#' @param tracts_gr
#'
#' @returns
#' @export
#'
#' @examples
get_bin_mat_subtracts_emp <- function(tracts_gr){
  # iterate over each chromosome
  # to get one matrix for each chromosome
  col_names <- unique(tracts_gr$name)
  bin_mat_l <- lapply(as.character(unique(seqnames(tracts_gr))),
                      function(chrom) {
                        chrom_gr <- filter(tracts_gr, seqnames==chrom)
                        mat <- get_bin_mat_subtracts_chrom(chrom_gr)
                        # This is because the some inds might not have tracts
                        # for some chromosomes and in order to create one big matrix it is important
                        # to have the dimensions consistent
                        full_matrix <- matrix(0, nrow = nrow(mat), ncol = length(col_names), dimnames = list(rownames(mat), col_names))
                        full_matrix[rownames(mat), colnames(mat)] <- mat
                        full_matrix
                      })
  bin_mat <- do.call(rbind, bin_mat_l)
  return(bin_mat)
}

#' Get binary matrix of tracts for recombination sites approach - SIMULATION
#'
#' @param tracts_gr
#'
#' @returns Binary
#' @export
#'
#' @examples
get_bin_mat_sites_sim <- function(tracts_gr){
  bin_mat <- sapply(unique(tracts_gr$name), function(sample_id){
    get_rec_sites(tracts_gr, sample_id)
  })
  bin_mat <- bin_mat[rowSums(bin_mat)!=0,]
  return(bin_mat)
}

#' Get binary matrix of tracts for recombination sites approach - EMPIRICAL
#'
#' @param tracts_gr
#'
#' @returns
#' @export
#'
#' @examples
get_bin_mat_sites_emp <- function(tracts_gr){
  # iterate over each chromosome
  # to get one matrix for each chromosome
  col_names <- unique(tracts_gr$name)
  bin_mat_l <- lapply(as.character(unique(GenomicRanges::seqnames(tracts_gr))),
                      function(chrom) {
                        chrom_gr <- plyranges::filter(tracts_gr, seqnames==chrom)
                        # Just calling the function for the simulation cause it si desginmed to work with only one chromosome
                        mat <- get_bin_mat_sites_sim(chrom_gr)
                        # This is because the some inds might not have tracts
                        # for some chromosomes and in order to create one big matrix it is important
                        # to have the dimensions consistent
                        full_matrix <- matrix(0, nrow = nrow(mat), ncol = length(col_names), dimnames = list(rownames(mat), col_names))
                        full_matrix[rownames(mat), colnames(mat)] <- mat
                        full_matrix
                      })
  bin_mat <- do.call(rbind, bin_mat_l)
  return(bin_mat)
}

get_rec_sites <- function(tracts_gr, sample_id){
  tract_df <- tracts_gr %>% as.data.frame() %>% dplyr::select(start, end)
  # convert to matrix
  tract_mat <- tract_df %>% as.matrix()
  # get vector of starts and ends
  starts_ends <- tract_mat %>% c() %>% unique() %>% sort()
  # it is meant to be for only one chrom, so it the gr has to be only one chr
  seqname <- GenomicRanges::seqnames(tracts_gr) %>% as.character() %>% unique()
  # The start end vector is sorted, so for the start I will exclude the last element, and fot he ends of the ranges the first one
  names <- paste0(seqname, ":", as.character(starts_ends))
  vec <- rep(0, length(names))
  # get which ranges are in that sample
  samp_tracts_gr <- tracts_gr %>% plyranges::filter(name==sample_id)
  samp_start_ends <- samp_tracts_gr %>% as.data.frame() %>% dplyr::select(start, end) %>% as.matrix() %>% c() %>% unique() %>% sort()
  names(vec) <- names
  vec[paste0(seqname, ":", as.character(samp_start_ends))] <- 1
  return(vec)
}

# Helper functions that are not exported

# function to get which subtracts are contained in an individual
# NOTE: the function is meant to work with one sample at a time
get_tract_on <- function(tracts_gr, sample_id){
  tract_df <- tracts_gr %>% as.data.frame() %>% select(start, end)
  # convert to matrix
  tract_mat <- tract_df %>% as.matrix()
  # get vector of starts and ends
  starts_ends <- tract_mat %>% c() %>% unique() %>% sort()
  # it is meant to be for only one chrom, so it the gr has to be only one chr
  seqname <- seqnames(tracts_gr) %>% as.character() %>% unique()
  # The start end vector is sorted, so for the start I will exclude the last element, and fot he ends of the ranges the first one
  mask_gr <- GRanges(seqnames = seqname, ranges = IRanges(start = starts_ends[-length(starts_ends)], end = starts_ends[-1]))
  names <- as.character(mask_gr)
  vec <- rep(0, length(names))
  # get which ranges are in that sample
  samp_tracts_gr <- tracts_gr %>% filter(name==sample_id)
  hit <- subjectHits(findOverlaps(samp_tracts_gr, mask_gr, minoverlap = 2L))
  vec[hit] <- 1
  names(vec) <- names
  return(vec)
}

get_bin_mat_subtracts_chrom <- function(tracts_gr){
  bin_mat <- sapply(unique(tracts_gr$name), function(sample_id){
    get_tract_on(tracts_gr, sample_id)
  })
  # return NULL if the bin mat has nothing inisde
  # this to avoid the problem of rowSums over a NULL elements
  if(is.null(dim(bin_mat))) return(NULL)
  bin_mat <- bin_mat[rowSums(bin_mat)!=0,]
  return(bin_mat)
}

# function that generates the bins for a simulated chromosome
generate_windows_sim <- function(n_chr=1, len_chr=100e6, window_size, step_size) {

  seq_info <- Seqinfo(paste0("chr", 1:n_chr), seqlengths = len_chr, isCircular = FALSE, genome = "sim")
  gr <- GRanges(paste0("chr", 1:n_chr), IRanges(rep(1, n_chr), len_chr), seqinfo = seq_info)

  windows_grl <- slidingWindows(gr, width = window_size, step = step_size)
  for (i in seq_along(windows_grl)) {
    windows_grl[[i]]$midpoint <- (start(windows_grl[[i]]) + end(windows_grl[[i]])) / 2
  }
  unlist(windows_grl)
}

# function to generate the windows binning of empirical chromosome
generate_windows_no_gap <- function(window_size, step_size) {
  autosomes_gr <- GenomeInfoDb::getChromInfoFromUCSC("hg19") %>%
    dplyr::filter(grepl("chr\\d+$", chrom)) %>%
    dplyr::mutate(start = 1) %>%
    dplyr::rename(end = size) %>%
    makeGRangesFromDataFrame()
  # seqinfo(autosomes_gr) <- seqinfo(tracts_gr)

  windows_grl <- slidingWindows(autosomes_gr, width = window_size, step = step_size)
  for (i in seq_along(windows_grl)) {
    windows_grl[[i]]$midpoint <- (start(windows_grl[[i]]) + end(windows_grl[[i]])) / 2
  }

  unlist(windows_grl)
}

# function to obtain a binary vector given the tracts and the windows genomic ranges
compute_bin_vec <- function(tracts_gr, windows_gr, sample_id = NULL, ncores = 1) {
  if(is.null(sample_id)){
    stop("Forgot to specify sampleID.")
  }

  # Get sample
  sample_gr <- tracts_gr %>% filter(name==sample_id)
  seqinfo(sample_gr) <- seqinfo(windows_gr)
  # first compute coverage...
  cov <- coverage(sample_gr)

  # get a list where each element correspond to a chromosome with the coverage of the tracts for each bin
  tract_cov_list <- parallel::mclapply(as.character(unique(seqnames(windows_gr))),
                             function(chrom) {
                               chrom_coverage <- cov[[chrom]]
                               chrom_gr <- windows_gr[seqnames(windows_gr) == chrom]
                               # count overlaps between windows and tracts
                               # check Views object
                               cov_view <- Views(chrom_coverage, ranges(chrom_gr))
                               average_coverage_per_window <- viewMeans(cov_view)
                               average_coverage_per_window
                             }, mc.cores = ncores)

  unlist(tract_cov_list)
}

# function that returns a vector of which exact tracts are in that sample
get_exact_match <- function(target, uni_tracts_gr){
  # binary vector of results
  vector <- rep(0, length(uni_tracts_gr))
  names(vector) <- as.character(uni_tracts_gr)
  # pairs <- findOverlapPairs(target, uni_tracts_gr, type="equal")
  res <- findOverlaps(target, uni_tracts_gr, type = "equal")
  vector[subjectHits(res)] <- 1
  vector
}
