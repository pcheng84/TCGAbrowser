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
hello <- function(who, how=c("shout", "whisper")) {
  stopifnot(
    is.character(who),
    length(who) == 1,
    !is.na(who)
  )
  how <- match.arg(how)
  fun <- switch(how, shout=shout, whisper=whisper)
  sprintf("Hello, %s, your name is %d characters", fun(who), nchar(who))
}

shout <- function(who) {
  toupper(who)
}

whisper <- function(who) {
  tolower(who)
}
