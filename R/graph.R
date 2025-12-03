# function  to get the adjacency matrix from tracts
# to build tract graph

#' Get adjacency matrix
#'
#'  Given a GRanges object containing tracts and a genomic position, the function returns the corresponding adjacency matrix for the tracts
#'  intersecting that genomic location. The matrix is symmetric and the row and column names correspond to the ID of the individuals carrying the tracts.
#'  Two individuals sharing one recombination breakpoints will have an entry equal to 1, two individuals sharing two rec breakpoints (i.e., they share the same tract) have an entry equal to 2.
#'
#' @param pos
#' @param tracts_gr
#' @param chrom
#'
#' @returns Matrix.
#' @export
#'
#' @examples
get_adj_mat <- function(tracts_gr, chrom, pos){
  # extract the tracts intersecting with that genomic location (pos)
  # my_tracts <- tracts_gr %>% plyranges::filter(seqnames==chrom, (pos-start)>=0, (pos-end)<0) %>% base::as.data.frame()
  my_tracts <- tracts_gr %>% as.data.frame() %>%  dplyr::filter(seqnames==chrom, (pos-start)>=0, (pos-end)<0)
  nodes <- my_tracts$name
  my_tracts <- my_tracts %>% dplyr::select(start, end)
  # convert to matrix
  tract_mat <- my_tracts %>% as.matrix()
  # get vector of starts and ends
  starts_ends <- tract_mat %>% c() %>% unique() %>% sort()
  tract_mat_new <- apply(tract_mat, MARGIN = 1, function(row) starts_ends %in% row %>% as.numeric) %>% t
  colnames(tract_mat_new) <- starts_ends
  # matrix multiplication
  adj_mat <- tract_mat_new%*%t(tract_mat_new)
  # since I do not want the link with itself (2) I set all the values = to 2 to 0
  # adj_mat[lower.tri(adj_mat, diag = TRUE)] <- 0
  diag(adj_mat) <- 0
  colnames(adj_mat) <- rownames(adj_mat) <- nodes
  adj_mat
}
