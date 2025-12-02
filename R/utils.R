#' Title
#'
#' @returns
#' @export
#'
#' @examples
read_metadata <- function() {
  # raw_info <- readr::read_tsv("data/neo.impute.1000g.sampleInfo_clusterInfo.txt")
  raw_info <-readr::read_tsv(system.file("extdata", "neo.impute.1000g.sampleInfo_clusterInfo.txt", package = "tractR"))
  info <-
    raw_info %>%
    dplyr::filter(groupAge != "Archaic") %>%
    dplyr::mutate(ageAverage = ifelse(groupAge == "Modern", 0, ageAverage)) %>%
    dplyr::mutate(coverage = ifelse(groupAge == "Modern", Inf, coverage))
  info
}

#' Title
#'
#' @param set
#'
#' @returns
#' @export
#'
#' @examples
read_tracts <- function(set = NULL) {
  info <- read_metadata()

  if(!is.null(set)){
    info <- dplyr::filter(info, groupAge == set)
  }

  raw_tracts <- readr::read_tsv(system.file("extdata", "Vindija33.19_raw_eurasian_wModern", package = "tractR"))

  tracts <- raw_tracts %>%
    dplyr::mutate(chrom = paste0("chr", chrom)) %>%
    dplyr::filter(ID %in% unique(info$sampleId)) %>%
    dplyr::select(ID, chrom, start, end) %>%
    dplyr::mutate(length = end - start, set = set)

  tracts
}


