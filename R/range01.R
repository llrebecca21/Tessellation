# Function that stores the range01 function
#' Title
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
range01 <- function(x){
  (x-min(x))/(max(x)-min(x))
}
