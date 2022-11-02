#' Geometric mean
#'
#' @param x your variable
#' @param na.rm  TRUE or FALSE
#'
#' @return a value
#' @export
#'
#' @examples gm_mean(data$variable); data %>% group_by(a, b, c, ...) %>% dplyr::summarise(g_mean = gm_mean(x, na.rm = TRUE))
gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
