#' Title
#'
#' @param ... Message you want to print.
#'
#' @returns Message with time.
#' @export
#'
#' @examples fun.message("Hello.")
fun.message <- function(...){
  cat(paste0("[",date(),"]"),...)
  cat("\n")
  flush.console()
}
