#' Estimate mode of a distribution
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}