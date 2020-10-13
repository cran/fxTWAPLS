# Load modern pollen data
modern_pollen <- read.csv(system.file("extdata", 
                                      "Modern_Pollen_gdd_alpha_Tmin.csv", 
                                      package = "fxTWAPLS", 
                                      mustWork = TRUE))

# Extract taxa
taxaColMin <- which(colnames(modern_pollen) == "Abies")
taxaColMax <- which(colnames(modern_pollen) == "Zygophyllaceae")
taxa <- modern_pollen[, taxaColMin:taxaColMax]

# Load reconstruction data
Holocene <- read.csv(system.file("extdata", "Holocene.csv", 
                                         package = "fxTWAPLS", 
                                         mustWork = TRUE), 
                     row.names = 1)
core <- Holocene[, -c(1:3)]

# Get the frequency of each climate variable fx
fx_Tmin <- fxTWAPLS::fx(modern_pollen$Tmin, bin = 0.02)
fx_gdd <- fxTWAPLS::fx(modern_pollen$gdd, bin = 20)
fx_alpha <- fxTWAPLS::fx(modern_pollen$alpha, bin = 0.002)

# MTCO
fit_Tmin <- fxTWAPLS::WAPLS.w(taxa, modern_pollen$Tmin, nPLS = 5)
fit_t_Tmin <- fxTWAPLS::TWAPLS.w(taxa, modern_pollen$Tmin, nPLS = 5)
fit_f_Tmin <- fxTWAPLS::WAPLS.w(taxa, 
                                modern_pollen$Tmin, 
                                nPLS = 5, 
                                usefx = TRUE, 
                                fx = fx_Tmin)
fit_tf_Tmin <- fxTWAPLS::TWAPLS.w(taxa, 
                                  modern_pollen$Tmin, 
                                  nPLS = 5, 
                                  usefx = TRUE, 
                                  fx = fx_Tmin)

test_that("fx works", {
  expect_equal(length(fx_Tmin), 6458)
  expect_equal(length(fx_gdd), 6458)
  expect_equal(length(fx_alpha), 6458)
})

test_that("WPALS training function works", {
  expected_names <- c("fit", "x", "taxon_name", "optimum", "comp", "u", 
                      "z", "s", "orth", "alpha", "meanx", "nPLS")
  expect_equal(names(fit_Tmin), expected_names)
  expect_equal(names(fit_f_Tmin), expected_names)
})

test_that("TWPALS training function works", {
  expected_names <- c("fit", "x", "taxon_name", "optimum", "comp", "u", 
                      "t", "z", "s", "orth", "alpha", "meanx", "nPLS")
  expect_equal(names(fit_t_Tmin), expected_names)
  expect_equal(names(fit_tf_Tmin), expected_names)
})

test_that("LOOCV as in rioja works", {
  # MTCO
  test_it <- 5 # Number of iterations for testing mode
  cv_Tmin <- fxTWAPLS::cv.w(taxa,
                            modern_pollen$Tmin,
                            nPLS = 5,
                            fxTWAPLS::WAPLS.w,
                            fxTWAPLS::WAPLS.predict.w,
                            cpus = 1,
                            test_mode = TRUE,
                            test_it = test_it)
  expect_equal(dim(cv_Tmin), c(test_it, 6))
})

test_that("Random t-test works", {
  # MTCO
  test_it <- 5 # Number of iterations for testing mode
  cv_Tmin <- fxTWAPLS::cv.w(taxa,
                            modern_pollen$Tmin,
                            nPLS = 5,
                            fxTWAPLS::WAPLS.w,
                            fxTWAPLS::WAPLS.predict.w,
                            cpus = 1,
                            test_mode = TRUE,
                            test_it = test_it)
  expect_output(rand_Tmin <- fxTWAPLS::rand.t.test.w(cv_Tmin, n.perm = 999))
  expect_equal(dim(cv_Tmin), c(test_it, 6))
})

test_that("Get distance between points works", {
  N <- 100 # Subset
  point <- modern_pollen[1:N, c("Long", "Lat")]
  expect_output(dist <- get_distance(point, cpus = 1))
  expect_equal(dim(dist), c(N, N))
})

test_that("Pseudo removed works", {
  N <- 100 # Subset
  point <- modern_pollen[1:N, c("Long", "Lat")]
  dist <- get_distance(point, cpus = 1)
  pseudo_Tmin <- fxTWAPLS::get_pseudo(dist, modern_pollen$Tmin[1:N], cpus = 1)
  expect_equal(length(pseudo_Tmin), N)
})

test_that("Pseudo removed LOOCV works", {
  test_it <- 5 # Number of iterations for testing mode
  point <- modern_pollen[, c("Long", "Lat")]
  dist <- get_distance(point, 
                       cpus = 1, 
                       test_mode = TRUE,
                       test_it = test_it)
  pseudo_Tmin <- fxTWAPLS::get_pseudo(dist, modern_pollen$Tmin, 
                                      cpus = 1, 
                                      test_mode = TRUE,
                                      test_it = test_it)
  # Test WAPLS
  cv_pr_Tmin <- fxTWAPLS::cv.pr.w(taxa,
                                  modern_pollen$Tmin,
                                  nPLS = 5,
                                  fxTWAPLS::TWAPLS.w,
                                  fxTWAPLS::TWAPLS.predict.w,
                                  pseudo_Tmin,
                                  cpus = 1,
                                  test_mode = TRUE,
                                  test_it = test_it)
  # Test TWAPLS
  cv_pr_t_Tmin <- fxTWAPLS::cv.pr.w(taxa,
                                    modern_pollen$Tmin,
                                    nPLS = 5,
                                    fxTWAPLS::TWAPLS.w,
                                    fxTWAPLS::TWAPLS.predict.w,
                                    pseudo_Tmin,
                                    cpus = 1, 
                                    test_mode = TRUE,
                                    test_it = test_it)
  expect_equal(dim(cv_pr_Tmin), c(test_it, 6))
  expect_equal(dim(cv_pr_t_Tmin), c(test_it, 6))
})

test_that("SSE sample works", {
  test_it <- 5 # Number of iterations for testing mode
  # MTCO
  sse_Tmin_WAPLS <- sse.sample(modern_taxa = taxa,
                               modern_climate = modern_pollen$Tmin,
                               fossil_taxa = core,
                               trainfun = fxTWAPLS::WAPLS.w,
                               predictfun = fxTWAPLS::WAPLS.predict.w,
                               nboot = 100,
                               nPLS = 5,
                               nsig = 3,
                               usefx = FALSE,
                               fx = NA,
                               cpus = 1,
                               seed = 1,
                               test_mode = TRUE,
                               test_it = test_it)
      
})

test_that("plot training data works", {
  expect_true(plot_train(fit_Tmin, 3))
  expect_true(plot_train(fit_f_Tmin, 3))
  expect_true(plot_train(fit_t_Tmin, 3))
  expect_true(plot_train(fit_tf_Tmin, 3))
})

test_that("plot residuals works", {
  expect_true(plot_residuals(fit_Tmin, 3))
  expect_true(plot_residuals(fit_f_Tmin, 3))
  expect_true(plot_residuals(fit_t_Tmin, 3))
  expect_true(plot_residuals(fit_tf_Tmin, 3))
})