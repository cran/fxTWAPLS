test_that("hexagonal logo works", {
  hex_logo(output = "hex_logo.png")
  expect_true(file.exists("hex_logo.png"))
  expect_false(dir.exists("hex_logo.png"))
  expect_gt(file.size("hex_logo.png"), 0)
  file.remove("hex_logo.png")
  expect_false(file.exists("hex_logo.png"))
})

test_that("parallel benchmark works", {
  # Define toy function that sleeps for (60/cpus) seconds
  a <- function(cpus) {Sys.sleep(4/cpus)}
  times_df <- par_benchmark(c(1, 2, 4), a, quiet = TRUE)
  expect_equal(length(times_df$times), 3)
  expect_output(par_benchmark(c(4), a, quiet = FALSE))
  expect_output(par_benchmark(c(4), a, quiet = FALSE, plot = TRUE))
  # print(list.files("."))
  # is.windows <- Sys.info()['sysname'] == "Windows"
  # # print(paste0("Is windows ", is.windows))
  # print(paste0("OS: ", Sys.info()['sysname']))
  # if (!is.windows) {
  #   expect_true(file.exists("./Rplots.pdf"))
  #   expect_false(dir.exists("./Rplots.pdf"))
  #   expect_gt(file.size("Rplots.pdf"), 0)
  # }
  # file.remove("./Rplots.pdf")
  # expect_false(file.exists("./Rplots.pdf"))
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
