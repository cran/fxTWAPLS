#' Perform parallel benchmarks 
#' 
#' Perform parallel benchmarks on a function and generate a plot with execution 
#' times vs CPU count.
#'
#' @param CPUS Vector with the number of CPUs.
#' @param FUN Parallel function, MUST have a parameter called \code{cpus}.
#' @param plot Boolean flag to request a plot for the results.
#' @param quiet Boolean flag to print results of each execution.
#' @param ... Optional arguments for \code{FUN}, must be named; e.g. 
#'     \code{x = test_df}.
#'
#' @examples
#' # Define toy function that sleeps for (2/cpus) seconds
#' a <- function(cpus) {Sys.sleep(2/cpus)}
#' fxTWAPLS:::par_benchmark(c(1, 2), a)
#' \donttest{
#' fxTWAPLS:::par_benchmark(c(1, 2), a, plot = TRUE)
#' }
#' 
#' @keywords internal
#' @noRd
par_benchmark <- function(CPUS, FUN, plot = FALSE, quiet = FALSE, ...) {
  cpus <- NULL # Local binding
  tictoc::tic.clearlog()
  for (c in CPUS) {
    tictoc::tic(paste0("Using ", c, " CPUs"))
    out <- FUN(..., cpus = c)
    tictoc::toc(log = TRUE, quiet = quiet)
  }
  times <- unlist(tictoc::tic.log(format = TRUE))
  times <- gsub(" sec elapsed", "", unlist(times))
  times <- gsub(".*: ", "", unlist(times))
  times <- as.numeric(times)
  times_df <- data.frame(cpus = CPUS, times = times)
  
  if(plot) {
    print(ggplot2::qplot(cpus, times, data = times_df) + 
            ggplot2::geom_area(alpha = 0.5) + 
            ggplot2::geom_line() + 
            ggplot2::labs(x = "CPUs", y = "Execution time [seconds]") +
            ggplot2::scale_x_continuous(breaks = 1:max(CPUS)) + 
            ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(10))
    )
  }
  return(times_df)
}

#' Combine results with progress bar
#' 
#' Combine results with progress bar, to be used in combination with 
#'     \code{\link[foreach::foreach]{foreach::foreach}}.
#'     
#' @importFrom utils flush.console
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' 
#' @param iterator Number of iterations.
#' @param FUN Function to combine the results (default: \code{rbind}).
#' @param ... Optional parameters.
#' 
#' @examples
#' \donttest{
#' # Load binary operator for backend
#' `%do%` <- foreach::`%do%`
#' N <- 5
#' out <- foreach::foreach(i = 1:N, 
#'                         .combine = fxTWAPLS:::comb_pb(N)) %do% {
#'                         Sys.sleep(1)
#'                         i
#'                        }
#' }
#' 
#' @noRd
#' @keywords internal
comb_pb <- function(iterator, FUN = rbind, ...) {
  pb <- txtProgressBar(min = 1, max = iterator - 1, style = 3)
  count <- 0
  function(...) {
    count <- count + length(list(...)) - 1
    setTxtProgressBar(pb, count)
    flush.console()
    FUN(...) # this can feed into .combine option of foreach
  }
}

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
