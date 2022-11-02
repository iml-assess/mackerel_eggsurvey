# Calculate mode of numeric or character vector
#' Title
#'
#' @param v numeric or character vector
#'
#' @return value with the highest number of occurences in v
#' @export
#'
#' @examples a <- c(1,2,3,4,5,5,6);getmode(a)
#'  b <- c("dog","cat","dog","bird");getmode(b)
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
