test_that("parallel benchmark works", {
  # Define toy function that sleeps for (60/cpus) seconds
  a <- function(cpus) {Sys.sleep(2/cpus)}
  times_df <- par_benchmark(c(1, 2), a, quiet = TRUE)
  expect_equal(length(times_df$times), 2)
  expect_output(par_benchmark(c(2), a, quiet = FALSE))
  expect_output(par_benchmark(c(2), a, quiet = FALSE, plot = TRUE))
})

test_that("combine with progress bar works", {
  # Load binary operator for backend
  `%do%` <- foreach::`%do%`
  N <- 5
  out <- foreach::foreach(i = 1:N, 
                          .combine = comb_pb(N)) %do% {
                            Sys.sleep(1)
                            i
                          }
  expect_equal(c(unname(out)), 1:N)
})