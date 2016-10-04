#' samnum function
#'
#' Counts how many samples are in each assay
#'
#' @param data a multiassayexperiment object
#'
#' @import dplyr
#' @import magrittr
#' @import MultiAssayExperiment
#'
#' @return data frame with assay names and number of samples in each assay
#'
#' @examples
#'
#' samnum(data)
#'
#' @export
#'
samnum <- function(data) {
  tbl_df(sampleMap(data)) %>% group_by(assay) %>% summarise(n())
}
