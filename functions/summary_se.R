#' Summarise grouped data
#' Function calculates summary statistics for grouped datasets including the mean, sd, n, se, and ci
#' @param .data 
#' @param measure_var the variable of interest
#' @param ... 
#' @param .ci specificy ci level wanted
#' @param na.rm 
#'
#' @return returns mean, sd, n, se, and ci
#' @export
#'
#' @examples mackerel_ziff %>% summary_se( pds_vif, year, div_caa)
summary_se <- function(.data, measure_var, ..., .ci = 0.95, na.rm = FALSE) {
  
  measure_var <- dplyr::enquo(measure_var)
  group_var <- dplyr::enquos(...)
  
  .data %>%
    group_by(!!! group_var) %>%
    summarise(mean = mean(!! measure_var, na.rm = na.rm),
              sd = sd(!! measure_var, na.rm = na.rm),
              n = n(),
              se = sd/sqrt(n),
              ci = se * qt(.ci/2 + 0.5, n-1)) %>%
    ungroup()
  
}

