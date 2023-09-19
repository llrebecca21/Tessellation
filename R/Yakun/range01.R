#' range01 function
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
range01 <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
