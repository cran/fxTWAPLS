#' Show progress bar
#'
#' @param expr R expression.
#' @inheritDotParams progressr::with_progress -handlers
#'
#' @return Return data from the function called.
#' @export
pb <- function(expr, ...) {
  progress_bar <-
    progressr::handler_progress(format = "(:current/:total) [:bar] :percent")
  progressr::with_progress(expr,
                           ...,
                           handlers = progress_bar)
}
