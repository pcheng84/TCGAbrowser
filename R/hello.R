#' Hello function
#'
#' Outputs your name and counts characters
#'
#' @param who character(1) holds name of person
#'
#' @return character(1) displays name and number of characters
#'
#' @examples
#' hello("Martin")
#'
#' @export
hello <- function(who) {
  sprintf("Hello, %s, your name is %d characters", upper(who), nchar(who))
}
