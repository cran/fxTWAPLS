#' @keywords internal
"_PACKAGE"

#' Get frequency of the climate value
#' 
#' Function to get the frequency of the climate value, which will be used to 
#'     provide \code{fx} correction for WA-PLS and TWA-PLS.
#'
#' @importFrom graphics plot
#' 
#' @param x Numeric vector with the modern climate values.
#' @param bin Binwidth to get the frequency of the modern climate values.
#' @param show_plot Boolean flag to show a plot of \code{fx ~ x}.
#'
#' @return Numeric vector with the frequency of the modern climate values.
#' @export
#'
#' @examples
#' \dontrun{
#' # Load modern pollen data
#' modern_pollen <- read.csv("/path/to/modern_pollen.csv")
#' 
#' # Get the frequency of each climate variable fx
#' fx_Tmin <- fxTWAPLS::fx(modern_pollen$Tmin, bin = 0.02)
#' fx_gdd <- fxTWAPLS::fx(modern_pollen$gdd, bin = 20)
#' fx_alpha <- fxTWAPLS::fx(modern_pollen$alpha, bin = 0.002)
#' }
#' 
#' @seealso \code{\link{cv.w}}, \code{\link{cv.pr.w}}, and 
#'     \code{\link{sse.sample}}
fx <- function(x, bin, show_plot = FALSE) {
  pbin <- round((max(x) - min(x)) / bin, digits = 0)
  bin <- (max(x) - min(x)) / pbin
  hist <- hist(x, breaks = seq(min(x), max(x), by = bin), plot = show_plot)
  xbin <- seq(min(x) + bin / 2, max(x) - bin / 2, by = bin)
  counts <- hist[["counts"]]
  fx <- rep(NA, length(x))
  for (i in seq_len(length(x))) {
    fx[i] <- counts[which.min(abs(x[i] - xbin))]
  }
  if (any(fx == 0)) {
    print("Some x have a count of 0!")
  }
  if (show_plot)
    plot(fx ~ x)
  return(fx)
}

#' WA-PLS training function
#' 
#' WA-PLS training function, which can perform \code{fx} correction.
#' 
#' @importFrom stats lm
#' 
#' @param modern_taxa The modern taxa abundance data, each row represents a 
#'     sampling site, each column represents a taxon.
#' @param modern_climate The modern climate value at each sampling site.
#' @param nPLS The number of components to be extracted.
#' @param usefx Boolean flag on whether or not use \code{fx} correction.
#' @param fx The frequency of the climate value for \code{fx} correction: if 
#'     \code{usefx = FALSE}, this should be \code{NA}; otherwise, this should 
#'     be obtained from the \code{\link{fx}} function.
#'
#' @return A list of the training results, which will be used by the predict 
#'     function. Each element in the list is described below:
#'     \describe{
#'     \item{\code{fit}}{the fitted values using each number of components.}
#'     \item{\code{x}}{the observed modern climate values.}
#'     \item{\code{taxon_name}}{the name of each taxon.}
#'     \item{\code{optimum}}{the updated taxon optimum (u* in the WA-PLS 
#'     paper).}
#'     \item{\code{comp}}{each component extracted (will be used in step 7 
#'     regression).}
#'     \item{\code{u}}{taxon optimum for each component (step 2).}
#'     \item{\code{z}}{a parameter used in standardization for each component 
#'     (step 5).}
#'     \item{\code{s}}{a parameter used in standardization for each component 
#'     (step 5).}
#'     \item{\code{orth}}{a list that stores orthogonalization parameters 
#'     (step 4).}
#'     \item{\code{alpha}}{a list that stores regression coefficients (step 7).}
#'     \item{\code{meanx}}{mean value of the observed modern climate values.}
#'     \item{\code{nPLS}}{the total number of components extracted.}
#'     }
#'     
#' @export
#'
#' @examples
#' \dontrun{
#' # Load modern pollen data
#' modern_pollen <- read.csv("/path/to/modern_pollen.csv")
#'                                       
#' # Extract taxa
#' taxaColMin <- which(colnames(modern_pollen) == "taxa0")
#' taxaColMax <- which(colnames(modern_pollen) == "taxaN")
#' taxa <- modern_pollen[, taxaColMin:taxaColMax]
#' 
#' # Get the frequency of each climate variable fx
#' fx_Tmin <- fxTWAPLS::fx(modern_pollen$Tmin, bin = 0.02)
#' fx_gdd <- fxTWAPLS::fx(modern_pollen$gdd, bin = 20)
#' fx_alpha <- fxTWAPLS::fx(modern_pollen$alpha, bin = 0.002)
#' 
#' # MTCO
#' fit_Tmin <- fxTWAPLS::WAPLS.w(taxa, modern_pollen$Tmin, nPLS = 5)
#' fit_f_Tmin <- fxTWAPLS::WAPLS.w(taxa, 
#'                                 modern_pollen$Tmin, 
#'                                 nPLS = 5, 
#'                                 usefx = TRUE, 
#'                                 fx = fx_Tmin)
#' }
#' 
#' @seealso \code{\link{fx}}, \code{\link{TWAPLS.w}}, and
#'     \code{\link{WAPLS.predict.w}}
WAPLS.w <- function(modern_taxa, 
                    modern_climate, 
                    nPLS = 5, 
                    usefx = FALSE, 
                    fx = NA) {
  # Step 0. Centre the environmental variable by subtracting the weighted mean
  x <- modern_climate
  y <- modern_taxa
  y <- as.matrix(y)
  nc <- ncol(modern_taxa)
  nr <- nrow(modern_taxa)
  Ytottot <- sum(y)
  sumk_yik <- rowSums(y)
  sumi_yik <- colSums(y)
  
  # Define some matrix to store the values
  u <- matrix(NA, nc, nPLS) # u of each component
  # u of each component, standardized the same way as r
  u_sd <- matrix(NA, nc, nPLS) 
  optimum <- matrix(NA, nc, nPLS) # u updated
  r <- matrix(NA, nr, nPLS) # site score
  z <- matrix(NA, 1, nPLS) # standardize
  s <- matrix(NA, 1, nPLS) # standardize
  orth <- list() # store orthogonalization parameters
  alpha <- list() # store regression coefficients
  comp <- matrix(NA, nr, nPLS) # each component
  fit <- matrix(NA, nr, nPLS) # current estimate
  
  pls <- 1
  # Step 1. Take the centred environmental variable(xi) as initial site 
  # scores (ri). 
  r[, pls] <- x - mean(x)
  
  # Step 2. Calculate new species scores (uk* by weighted averaging of the 
  # site scores)
  u[, pls] <- t(y) %*% x / sumi_yik # uk = sumi_yik*xi/sumi_yik; 1*nmodern_taxa
  
  # Step 3. Calculate new site scores (ri) by weighted averaging of the species 
  # scores
  r[, pls] <- y %*% u[, pls] / sumk_yik # xi=sumk_yik*uk/sumk_yik; 1*nsite
  
  # Step 4. For the first axis go to Step 5.
  
  # Step 5. Standardize the new site scores (ri) ter braak 1987 5.2.c
  z[, pls] <- mean(r[, pls], na.rm = TRUE)
  s[, pls] <- sqrt(sum((r[, pls] - z[, pls]) ^ 2, na.rm = TRUE) / Ytottot)
  r[, pls] <- (r[, pls] - z[, pls]) / s[, pls]
  
  # Step 6. Take the standardized score as the new component
  comp[, pls] <- r[, pls]
  
  # Step 7. Regress the environmental variable (xJ on the components obtained 
  # so far using weights and take the fitted values as current estimates
  if(usefx == FALSE) {
    lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], 
                    weights = sumk_yik / Ytottot)
  } else{
    lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], weights = 1 / fx ^ 2)
  }
  
  fit[, pls] <- lm[["fitted.values"]]
  alpha[[pls]] <- lm[["coefficients"]]
  u_sd[, pls] <- (u[, pls] - z[, pls]) / s[, pls]
  optimum[, pls] <- alpha[[pls]][1] + u_sd[, pls] * alpha[[pls]][2]
  
  for(pls in 2:nPLS) {
    # Go to Step 2 with the residuals of the regression as the new site 
    # scores (rJ.
    r[, pls] <- lm[["residuals"]]
    
    # Step 2. Calculate new species scores (uk* by weighted averaging of the 
    # site scores)
    # uk=sumi_yik*xi/sumi_yik; 1*nmodern_taxa
    u[, pls] <- t(y) %*% r[, pls] / sumi_yik
    
    # Step 3. Calculate new site scores (ri) by weighted averaging of the 
    # species scores
    # xi=sumk_yik*uk/sumk_yik; 1*nsite
    r[, pls] <- y %*% u[, pls] / sumk_yik 
    
    # Step 4. For second and higher components, make the new site scores (ri) 
    # uncorrelated with the previous components by orthogonalization 
    # (Ter Braak, 1987 : Table 5 .2b)
    v <- rep(NA, pls - 1)
    for (j in 1:(pls - 1)) {
      fi <- r[, pls - j]
      xi <- r[, pls]
      v[pls - j] <- sum(sumk_yik * fi * xi) / Ytottot
      xinew <- xi - v[pls - j] * fi
    }
    orth[[pls]] <- v
    r[, pls] <- xinew
    
    # Step 5. Standardize the new site scores (ri) ter braak 1987 5.2.c
    z[, pls] <- mean(r[, pls], na.rm = TRUE)
    s[, pls] <- sqrt(sum((r[, pls] - z[, pls]) ^ 2, na.rm = TRUE) / Ytottot)
    r[, pls] <- (r[, pls] - z[, pls]) / s[, pls]
    
    # Step 6. Take the standardized score as the new component
    comp[, pls] <- r[, pls]
    
    # Step 7. Regress the environmental variable on the components obtained so 
    # far using weights and take the fitted values as current estimates 
    if(usefx == FALSE) {
      lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], 
                      weights = sumk_yik / Ytottot)
    } else{
      lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], weights = 1 / fx ^ 2)
    }
    
    fit[, pls] <- lm[["fitted.values"]]
    alpha[[pls]] <- lm[["coefficients"]]
    
    u_sd[, pls] <- (u[, pls] - z[, pls]) / s[, pls]
    optimum[, pls] <-
      alpha[[pls]][1] + u_sd[, 1:pls] %*% as.matrix(alpha[[pls]][2:(pls + 1)])
  }
  
  list <- list(fit, 
               modern_climate, 
               colnames(modern_taxa), 
               optimum,comp, 
               u, 
               z, 
               s, 
               orth, 
               alpha, 
               mean(modern_climate), 
               nPLS)
  names(list) <- c("fit", 
                   "x", 
                   "taxon_name", 
                   "optimum", 
                   "comp", 
                   "u", 
                   "z", 
                   "s", 
                   "orth", 
                   "alpha", 
                   "meanx", 
                   "nPLS")
  return(list)
}

#' TWA-PLS training function
#' 
#' TWA-PLS training function, which can perform \code{fx} correction.
#' 
#' @importFrom stats lm
#' 
#' @inheritParams WAPLS.w
#'
#' @return A list of the training results, which will be used by the predict 
#'     function. Each element in the list is described below:
#'     \describe{
#'     \item{\code{fit}}{the fitted values using each number of components.}
#'     \item{\code{x}}{the observed modern climate values.}
#'     \item{\code{taxon_name}}{the name of each taxon.}
#'     \item{\code{optimum}}{the updated taxon optimum}
#'     \item{\code{comp}}{each component extracted (will be used in step 7 
#'     regression).}
#'     \item{\code{u}}{taxon optimum for each component (step 2).}
#'     \item{\code{t}}{taxon tolerance for each component (step 2).}
#'     \item{\code{z}}{a parameter used in standardization for each component 
#'     (step 5).}
#'     \item{\code{s}}{a parameter used in standardization for each component 
#'     (step 5).}
#'     \item{\code{orth}}{a list that stores orthogonalization parameters 
#'     (step 4).}
#'     \item{\code{alpha}}{a list that stores regression coefficients (step 7).}
#'     \item{\code{meanx}}{mean value of the observed modern climate values.}
#'     \item{\code{nPLS}}{the total number of components extracted.}
#'     }
#'     
#' @export
#'
#' @examples
#' \dontrun{
#' # Load modern pollen data
#' modern_pollen <- read.csv("/path/to/modern_pollen.csv")
#'                                       
#' # Extract taxa
#' taxaColMin <- which(colnames(modern_pollen) == "taxa0")
#' taxaColMax <- which(colnames(modern_pollen) == "taxaN")
#' taxa <- modern_pollen[, taxaColMin:taxaColMax]
#' 
#' # Get the frequency of each climate variable fx
#' fx_Tmin <- fxTWAPLS::fx(modern_pollen$Tmin, bin = 0.02)
#' fx_gdd <- fxTWAPLS::fx(modern_pollen$gdd, bin = 20)
#' fx_alpha <- fxTWAPLS::fx(modern_pollen$alpha, bin = 0.002)
#' 
#' # MTCO
#' fit_t_Tmin <- fxTWAPLS::TWAPLS.w(taxa, modern_pollen$Tmin, nPLS = 5)
#' fit_tf_Tmin <- fxTWAPLS::TWAPLS.w(taxa, 
#'                                   modern_pollen$Tmin, 
#'                                   nPLS = 5, 
#'                                   usefx = TRUE, 
#'                                   fx = fx_Tmin)
#' }
#' 
#' @seealso \code{\link{fx}}, \code{\link{TWAPLS.predict.w}}, and
#'     \code{\link{WAPLS.w}}
TWAPLS.w <- function(modern_taxa,
                     modern_climate,
                     nPLS = 5,
                     usefx = FALSE,
                     fx = NA){
  # Step 0. Centre the environmental variable by subtracting the weighted mean
  x <- modern_climate
  y <- modern_taxa
  y <- as.matrix(y)
  nc <- ncol(modern_taxa)
  nr <- nrow(modern_taxa)
  Ytottot <- sum(y)
  sumk_yik <- rowSums(y)
  sumi_yik <- colSums(y)
  
  #Define some matrix to store the values
  u <- matrix(NA, nc, nPLS) # u of each component
  # u of each component, standardized the same way as r
  u_sd <- matrix(NA, nc, nPLS)
  optimum <- matrix(NA, nc, nPLS) # u updated
  t <- matrix(NA, nc, nPLS) # tolerance
  r <- matrix(NA, nr, nPLS) # site score
  z <- matrix(NA, 1, nPLS) # standardize
  s <- matrix(NA, 1, nPLS) # standardize
  orth <- list() # store orthogonalization parameters
  alpha <- list() # store regression coefficients
  comp <- matrix(NA, nr, nPLS) # each component
  fit <- matrix(NA, nr, nPLS) # current estimate
  
  pls <- 1
  # Step 1. Take the centred environmental variable (xi) as initial site 
  # scores (ri). 
  r[, pls] <- x - mean(x)
  
  # Step 2. Calculate uk and tk
  u[, pls] <- t(y) %*% x / sumi_yik # uk=sumi_yik*xi/sumi_yik; 1*nmodern_taxa
  n2 <- matrix(NA, nc, 1)
  for (k in 1:nc) {
    t[k, pls] <- sqrt(sum(y[, k] * (x - u[k, pls]) ^ 2) / sumi_yik[k])
    n2[k] <- 1 / sum((y[, k] / sum(y[, k])) ^ 2)
    t[k, pls] <- t[k, pls] / sqrt(1 - 1 / n2[k])
  }
  
  # Step 3. Calculate new site scores (ri)
  #xi; 1*nsite
  r[, pls] <- (y %*% (u[, pls] / t[, pls] ^ 2)) / (y %*% (1 / t[, pls] ^ 2))
  
  # Step 4. For the first axis go to Step 5.
  
  # Step 5. Standardize the new site scores (ri) ter braak 1987 5.2.c
  z[, pls] <- mean(r[, pls], na.rm = TRUE)
  s[, pls] <- sqrt(sum((r[, pls] - z[, pls]) ^ 2, na.rm = TRUE) / Ytottot)
  r[, pls] <- (r[, pls] - z[, pls]) / s[, pls]
  
  # Step 6. Take the standardized score as the new component
  comp[, pls] <- r[, pls]
  
  # Step 7. Regress the environmental variable on the components obtained so far
  # using weights and take the fitted values as current estimates 
  if (usefx == FALSE) {
    lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], 
                    weights = sumk_yik / Ytottot)
  } else {
    lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], weights = 1 / fx ^ 2)
  }
  
  fit[, pls] <- lm[["fitted.values"]]
  alpha[[pls]] <- lm[["coefficients"]]
  u_sd[, pls] <- (u[, pls] - z[, pls]) / s[, pls]
  optimum[, pls] <- alpha[[pls]][1] + u_sd[, pls] * alpha[[pls]][2]
  
  for (pls in 2:nPLS) {
    # Go to Step 2 with the residuals of the regression as the new site 
    # scores (ri).
    r[, pls] <- lm[["residuals"]]
    
    # Step 2. Calculate new uk and tk
    # uk=sumi_yik*xi/sumi_yik; 1*nmodern_taxa
    u[, pls] <- t(y) %*% r[, pls] / sumi_yik
    n2 <- matrix(NA, nc, 1)
    for (k in 1:nc) {
      t[k, pls] <- sqrt(sum(y[, k] * (r[, pls] - u[k, pls]) ^ 2) / sumi_yik[k])
      n2[k] <- 1 / sum((y[, k] / sum(y[, k])) ^ 2)
      t[k, pls] <- t[k, pls] / sqrt(1 - 1 / n2[k])
    }
    
    # Step 3. Calculate new site scores (r;) by weighted averaging of the 
    # species scores, i.e. new
    r[,pls] <- (y%*%(u[,pls]/t[,pls]^2))/(y%*%(1/t[,pls]^2)); #xi; 1*nsite
    
    # Step 4. For second and higher components, make the new site scores (r;) 
    # uncorrelated with the previous components by orthogonalization 
    # (Ter Braak, 1987 : Table 5 .2b)
    v <- rep(NA, pls - 1)
    for (j in 1:(pls - 1)) {
      fi <- r[, pls - j]
      xi <- r[, pls]
      v[pls - j] <- sum(sumk_yik * fi * xi) / Ytottot
      xinew <- xi - v[pls - j] * fi
    }
    orth[[pls]] <- v
    # plot(xinew~r[,pls]);abline(0,1)
    r[, pls] <- xinew
    
    # Step 5. Standardize the new site scores (ri) ter braak 1987 5.2.c
    z[, pls] <- mean(r[, pls], na.rm = TRUE)
    s[, pls] <- sqrt(sum((r[, pls] - z[, pls]) ^ 2, na.rm = TRUE) / Ytottot)
    r[, pls] <- (r[, pls] - z[, pls]) / s[, pls]
    
    # Step 6. Take the standardized scores as the new component
    comp[, pls] <- r[, pls]
    
    # Step 7. Regress the environmental variable (xJ on the components obtained 
    # so far using weights and take the fitted values as current estimates 
    if (usefx == FALSE) {
      lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], 
                      weights = sumk_yik / Ytottot)
    } else {
      lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], weights = 1 / fx ^ 2)
    }
    
    fit[, pls] <- lm[["fitted.values"]]
    alpha[[pls]] <- lm[["coefficients"]]
    
    u_sd[, pls] <- (u[, pls] - z[, pls]) / s[, pls]
    optimum[, pls] <-
      alpha[[pls]][1] + u_sd[, 1:pls] %*% as.matrix(alpha[[pls]][2:(pls + 1)])
  }
  
  list <- list(fit, 
               modern_climate, 
               colnames(modern_taxa), 
               optimum, 
               comp, 
               u, 
               t, 
               z, 
               s, 
               orth, 
               alpha, 
               mean(modern_climate), 
               nPLS)
  names(list) <- c("fit", 
                   "x", 
                   "taxon_name", 
                   "optimum", 
                   "comp", 
                   "u", 
                   "t", 
                   "z", 
                   "s", 
                   "orth", 
                   "alpha", 
                   "meanx", 
                   "nPLS")
  return(list)
}

#' WA-PLS predict function
#'
#' @param WAPLSoutput The output of the \code{\link{WAPLS.w}} training function, 
#'     either with or without \code{fx} correction.
#' @param fossil_taxa Fossil taxa abundance data to reconstruct past climates, 
#'     each row represents a site to be reconstructed, each column represents a 
#'     taxon.
#'
#' @return A list of the reconstruction results. Each element in the list is 
#'     described below:
#'     \describe{
#'     \item{\code{fit}}{The fitted values using each number of components.}
#'     \item{\code{nPLS}}{The total number of components extracted.}
#'     }
#'     
#' @export
#'
#' @examples
#' \dontrun{
#' # Load modern pollen data
#' modern_pollen <- read.csv("/path/to/modern_pollen.csv")
#'                                       
#' # Extract taxa
#' taxaColMin <- which(colnames(modern_pollen) == "taxa0")
#' taxaColMax <- which(colnames(modern_pollen) == "taxaN")
#' taxa <- modern_pollen[, taxaColMin:taxaColMax]
#' 
#' # Load reconstruction data
#' Holocene <- read.csv("/path/to/Holocene.csv")
#' taxaColMin <- which(colnames(Holocene) == "taxa0")
#' taxaColMax <- which(colnames(Holocene) == "taxaN")
#' core <- Holocene[, taxaColMin:taxaColMax]
#' 
#' # Get the frequency of each climate variable fx
#' fx_Tmin <- fxTWAPLS::fx(modern_pollen$Tmin, bin = 0.02)
#' fx_gdd <- fxTWAPLS::fx(modern_pollen$gdd, bin = 20)
#' fx_alpha <- fxTWAPLS::fx(modern_pollen$alpha, bin = 0.002)
#' 
#' # MTCO
#' ## Train
#' fit_Tmin <- fxTWAPLS::WAPLS.w(modern_taxa = taxa, 
#'                               modern_climate = modern_pollen$Tmin, 
#'                               nPLS = 5)
#' fit_f_Tmin <- fxTWAPLS::WAPLS.w(modern_taxa = taxa, 
#'                                 modern_climate = modern_pollen$Tmin,
#'                                 nPLS = 5, 
#'                                 usefx = TRUE, 
#'                                 fx = fx_Tmin)
#' ## Predict
#' fossil_Tmin <- fxTWAPLS::WAPLS.predict.w(fit_Tmin, core)
#' fossil_f_Tmin <- fxTWAPLS::WAPLS.predict.w(fit_f_Tmin, core)
#' }
#' 
#' @seealso \code{\link{WAPLS.w}}
WAPLS.predict.w <- function(WAPLSoutput, fossil_taxa) {
  y <- fossil_taxa
  y <- as.matrix(y)
  nc <- ncol(fossil_taxa)
  nr <- nrow(fossil_taxa)
  Ytottot <- sum(y)
  sumk_yik <- rowSums(y)
  sumi_yik <- colSums(y)
  
  nPLS <- WAPLSoutput[["nPLS"]]
  meanx <- WAPLSoutput[["meanx"]]
  u <- WAPLSoutput[["u"]]
  z <- WAPLSoutput[["z"]]
  s <- WAPLSoutput[["s"]]
  orth <- WAPLSoutput[["orth"]]
  alpha <- WAPLSoutput[["alpha"]]
  
  if (nc != nrow(u)) {
    print("Number of taxa doesn't match!")
  }
  if (all(colnames(fossil_taxa) == WAPLSoutput[["taxon_name"]]) == FALSE) {
    print("Taxa don't match!")
  }
  
  # Define some matrix to store the values
  fit <- matrix(NA, nr, nPLS)
  r <- matrix(NA, nr, nPLS)
  comp <- matrix(NA, nr, nPLS)
  
  pls <- 1
  # xi=sumk_yik*uk/sumk_yik; 1*nsite
  r[, pls] <- y %*% u[, pls] / sumk_yik # xi=sumk_yik*uk/sumk_yik; 1*nsite
  # standardize the same way
  r[, pls] <- (r[, pls] - z[, pls]) / s[, pls]
  # multiply the same regression coefficients
  comp[, pls] <- r[, pls]
  fit[, 1] <- alpha[[pls]][1] + comp[, pls] * alpha[[pls]][2]
  
  for(pls in 2:nPLS) {
    # xi = sumk_yik*uk/sumk_yik; 1*nsite
    r[, pls] <- y %*% u[, pls] / sumk_yik # xi = sumk_yik*uk/sumk_yik; 1*nsite
    # orthoganlization the same way
    for (j in 1:(pls - 1)) {
      fi <- r[, pls - j]
      xi <- r[, pls]
      xinew <- xi - orth[[pls]][pls - j] * fi
    }
    r[, pls] <- xinew
    
    # standardize the same way
    r[, pls] <- (r[, pls] - z[, pls]) / s[, pls]
    # multiply the same regression coefficients
    comp[, pls] <- r[, pls]
    fit[, pls] <-
      alpha[[pls]][1] + comp[, 1:pls] %*% as.matrix(alpha[[pls]][2:(pls + 1)])
  }
  
  list <- list(fit, nPLS)
  names(list) <- c(c("fit", "nPLS"))
  return(list)
}

#' TWA-PLS predict function
#'
#' @param TWAPLSoutput The output of the \code{\link{TWAPLS.w}} training 
#'     function, either with or without \code{fx} correction.
#' @param fossil_taxa Fossil taxa abundance data to reconstruct past climates, 
#'     each row represents a site to be reconstructed, each column represents 
#'     a taxon.
#'
#' @return A list of the reconstruction results. Each element in the list is 
#'     described below:
#'     \describe{
#'     \item{\code{fit}}{the fitted values using each number of components.}
#'     \item{\code{nPLS}}{the total number of components extracted.}
#'     }
#'     
#' @export
#'
#' @examples
#' \dontrun{
#' # Load modern pollen data
#' modern_pollen <- read.csv("/path/to/modern_pollen.csv")
#'                                       
#' # Extract taxa
#' taxaColMin <- which(colnames(modern_pollen) == "taxa0")
#' taxaColMax <- which(colnames(modern_pollen) == "taxaN")
#' taxa <- modern_pollen[, taxaColMin:taxaColMax]
#' 
#' # Load reconstruction data
#' Holocene <- read.csv("/path/to/Holocene.csv")
#' taxaColMin <- which(colnames(Holocene) == "taxa0")
#' taxaColMax <- which(colnames(Holocene) == "taxaN")
#' core <- Holocene[, taxaColMin:taxaColMax]
#' 
#' # Get the frequency of each climate variable fx
#' fx_Tmin <- fxTWAPLS::fx(modern_pollen$Tmin, bin = 0.02)
#' fx_gdd <- fxTWAPLS::fx(modern_pollen$gdd, bin = 20)
#' fx_alpha <- fxTWAPLS::fx(modern_pollen$alpha, bin = 0.002)
#' 
#' # MTCO
#' ## Train
#' fit_t_Tmin <- fxTWAPLS::TWAPLS.w(modern_taxa = taxa, 
#'                                  modern_climate = modern_pollen$Tmin, 
#'                                  nPLS = 5)
#' fit_tf_Tmin <- fxTWAPLS::TWAPLS.w(modern_taxa = taxa, 
#'                                   modern_climate = modern_pollen$Tmin,
#'                                   nPLS = 5, 
#'                                   usefx = TRUE, 
#'                                   fx = fx_Tmin)
#'     
#' ## Predict
#' fossil_t_Tmin <- fxTWAPLS::TWAPLS.predict.w(fit_t_Tmin, core)
#' fossil_tf_Tmin <- fxTWAPLS::TWAPLS.predict.w(fit_tf_Tmin, core)
#' }
#' 
#' @seealso \code{\link{TWAPLS.w}}
TWAPLS.predict.w <- function(TWAPLSoutput, fossil_taxa) {
  y <- fossil_taxa
  y <- as.matrix(y)
  nc <- ncol(fossil_taxa)
  nr <- nrow(fossil_taxa)
  Ytottot <- sum(y)
  sumk_yik <- rowSums(y)
  sumi_yik <- colSums(y)
  
  nPLS <- TWAPLSoutput[["nPLS"]]
  meanx <- TWAPLSoutput[["meanx"]]
  u <- TWAPLSoutput[["u"]]
  t <- TWAPLSoutput[["t"]]
  z <- TWAPLSoutput[["z"]]
  s <- TWAPLSoutput[["s"]]
  orth <- TWAPLSoutput[["orth"]]
  alpha <- TWAPLSoutput[["alpha"]]
  
  if (nc != nrow(u)) {
    print("Number of taxa doesn't match!")
  }
  if (all(colnames(fossil_taxa) == TWAPLSoutput[["taxon_name"]]) == FALSE) {
    print("Taxa don't match!")
  }
  
  # Define some matrix to store the values
  fit <- matrix(NA, nr, nPLS)
  r <- matrix(NA, nr, nPLS)
  comp <- matrix(NA, nr, nPLS)
  
  pls <- 1
  # xi=sumk_yik*uk/sumk_yik; 1*nsite
  # xi; 1*nsite
  r[, pls] <- (y %*% (u[, pls] / t[, pls] ^ 2)) / (y %*% (1 / t[, pls] ^ 2))
  # standardize the same way
  r[, pls] <- (r[, pls] - z[, pls]) / s[, pls]
  # multiply the same regression coefficients
  comp[, pls] <- r[, pls]
  fit[, 1] <- alpha[[pls]][1] + comp[, pls] * alpha[[pls]][2]
  
  for(pls in 2:nPLS) {
    # xi=sumk_yik*uk/sumk_yik; 1*nsite
    # xi; 1*nsite
    r[, pls] <- (y %*% (u[, pls] / t[, pls] ^ 2)) / (y %*% (1 / t[, pls] ^ 2))
    # orthoganlization the same way
    for (j in 1:(pls - 1)) {
      fi <- r[, pls - j]
      xi <- r[, pls]
      xinew <- xi - orth[[pls]][pls - j] * fi
    }
    r[, pls] <- xinew
    
    # standardize the same way
    r[, pls] <- (r[, pls] - z[, pls]) / s[, pls]
    # multiply the same regression coefficients
    comp[, pls] <- r[, pls]
    fit[, pls] <-
      alpha[[pls]][1] + comp[, 1:pls] %*% as.matrix(alpha[[pls]][2:(pls + 1)])
  }
  
  list <- list(fit, nPLS)
  names(list) <- c(c("fit", "nPLS"))
  return(list)
}

#' Calculate Sample Specific Errors
#'
#' @param modern_taxa The modern taxa abundance data, each row represents a 
#'     sampling site, each column represents a taxon.
#' @param modern_climate The modern climate value at each sampling site
#' @param fossil_taxa Fossil taxa abundance data to reconstruct past climates, 
#'     each row represents a site to be reconstructed, each column represents a 
#'     taxon.
#' @param trainfun Training function you want to use, either 
#'     \code{\link{WAPLS.w}} or \code{\link{TWAPLS.w}}.
#' @param predictfun Predict function you want to use: if \code{trainfun} is 
#'     \code{\link{WAPLS.w}}, then this should be \code{\link{WAPLS.predict.w}}; 
#'     if \code{trainfun} is \code{\link{TWAPLS.w}}, then this should be 
#'     \code{\link{TWAPLS.predict.w}}.
#' @param nboot The number of bootstrap cycles you want to use.
#' @param nPLS The number of components to be extracted.
#' @param nsig The significant number of components to use to reconstruct past 
#'     climates, this can be obtained from the cross-validation results.
#' @param usefx Boolean flag on whether or not use \code{fx} correction.
#' @param fx The frequency of the climate value for \code{fx} correction: if 
#'     \code{usefx = FALSE}, this should be \code{NA}; otherwise, this should 
#'     be obtained from the \code{\link{fx}} function.
#' @param cpus Number of CPUs for simultaneous iterations to execute, check
#'     \code{parallel::detectCores()} for available CPUs on your machine.
#' @param seed Seed for reproducibility.
#' @param test_mode Boolean flag to execute the function with a limited number
#'     of iterations, \code{test_it}, for testing purposes only.
#' @param test_it Number of iterations to use in the test mode.
#'
#' @return The bootstrapped standard error for each site.
#' @export
#'
#' @examples
#' \dontrun{
#' # Load modern pollen data
#' modern_pollen <- read.csv("/path/to/modern_pollen.csv")
#'                                       
#' # Extract taxa
#' taxaColMin <- which(colnames(modern_pollen) == "taxa0")
#' taxaColMax <- which(colnames(modern_pollen) == "taxaN")
#' taxa <- modern_pollen[, taxaColMin:taxaColMax]
#' 
#' # Get the frequency of each climate variable fx
#' fx_Tmin <- fxTWAPLS::fx(modern_pollen$Tmin, bin = 0.02)
#' fx_gdd <- fxTWAPLS::fx(modern_pollen$gdd, bin = 20)
#' fx_alpha <- fxTWAPLS::fx(modern_pollen$alpha, bin = 0.002)
#' 
#' # Load reconstruction data
#' Holocene <- read.csv("/path/to/Holocene.csv")
#' taxaColMin <- which(colnames(Holocene) == "taxa0")
#' taxaColMax <- which(colnames(Holocene) == "taxaN")
#' core <- Holocene[, taxaColMin:taxaColMax]
#' 
#' # MTCO
#' ## fx
#' fit_Tmin <- fxTWAPLS::WAPLS.w(taxa, modern_pollen$Tmin, nPLS = 5)
#' ## SSE
#' nboot <- 5 # Recommended 1000
#' ### without fx
#' sse_Tmin_WAPLS <- fxTWAPLS::sse.sample(modern_taxa = taxa,
#'                                        modern_climate = modern_pollen$Tmin,
#'                                        fossil_taxa = core,
#'                                        trainfun = fxTWAPLS::WAPLS.w,
#'                                        predictfun = 
#'                                          fxTWAPLS::WAPLS.predict.w,
#'                                        nboot = nboot,
#'                                        nPLS = 5,
#'                                        nsig = 3,
#'                                        usefx = FALSE,
#'                                        fx = NA,
#'                                        cpus = 2,
#'                                        seed = 1)
#' ### with fx
#' sse_f_Tmin_WAPLS <- fxTWAPLS::sse.sample(modern_taxa = taxa,
#'                                          modern_climate = 
#'                                            modern_pollen$Tmin,
#'                                          fossil_taxa = core,
#'                                          trainfun = fxTWAPLS::WAPLS.w,
#'                                          predictfun = 
#'                                            fxTWAPLS::WAPLS.predict.w,
#'                                          nboot = nboot,
#'                                          nPLS = 5,
#'                                          nsig = 3,
#'                                          usefx = TRUE,
#'                                          fx = fx_Tmin,
#'                                          cpus = 2,
#'                                          seed = 1)
#' }
#' 
#' @seealso \code{\link{fx}}, \code{\link{TWAPLS.w}}, 
#'     \code{\link{TWAPLS.predict.w}}, \code{\link{WAPLS.w}}, and 
#'     \code{\link{WAPLS.predict.w}}
sse.sample <- function(modern_taxa,
                       modern_climate,
                       fossil_taxa,
                       trainfun,
                       predictfun,
                       nboot,
                       nPLS,
                       nsig,
                       usefx,
                       fx,
                       cpus = 4,
                       seed = NULL,
                       test_mode = FALSE,
                       test_it = 5) {
  i <- NULL # Local binding
  # Check the number of CPUs does not exceed the availability
  avail_cpus <- parallel::detectCores() - 1
  cpus <- ifelse(cpus > avail_cpus, avail_cpus, cpus)
  
  # Start parallel backend
  cl <- parallel::makeCluster(cpus)
  doParallel::registerDoParallel(cl)
  
  # Load binary operator for backend
  `%dopar%` <- foreach::`%dopar%`
  # `%dorng%` <- doRNG::`%dorng%`
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Make list of row numbers by sampling with 
  # replacement
  k_samples <- replicate(nboot, sample(1:nrow(modern_taxa),
                                       size = nrow(modern_taxa),
                                       replace = TRUE))
  
  # Create list of indices to loop through
  idx <- 1:nboot
  # Reduce the list of indices, if test_mode = TRUE
  if (test_mode) {
    idx <- 1:test_it
  }
  xboot <- foreach::foreach(i = idx,
                            .combine = cbind) %dopar% {
                              tryCatch({
                                # Extract list of row numbers by sampling with 
                                # replacement
                                k <- k_samples[, i]
                                
                                # Reorganise modern_taxa obs in k order
                                modern_taxak <- modern_taxa[k, ] 
                                modern_climatek <- modern_climate[k]
                                fxk<-fx[k]
                                col_not0 <- which(colSums(modern_taxak) > 0)
                                # Strip out zero-sum cols
                                modern_taxak <- modern_taxak[, col_not0]
                                # Apply train function, with modern_climate also 
                                # in k order
                                if (usefx == FALSE) {
                                  mod <- trainfun(modern_taxak, 
                                                  modern_climatek, 
                                                  nPLS = nPLS)
                                } else {
                                  mod <- trainfun(modern_taxak, 
                                                  modern_climatek, 
                                                  nPLS = nPLS, 
                                                  usefx = TRUE, 
                                                  fx = fxk)
                                }
                                
                                # Make reconstruction
                                predictfun(mod, 
                                           fossil_taxa[, col_not0])$fit[, nsig]
                              }, error=function(e){})
                            }
  parallel::stopCluster(cl) # Stop cluster
  
  avg.xboot <- rowMeans(xboot, na.rm = TRUE)
  v1 <- boot.mean.square <- rowMeans((xboot - avg.xboot) ^ 2 , na.rm = TRUE)
  return(sqrt(v1))
}

#' Leave-one-out cross-validation
#' 
#' Leave-one-out cross-validation as 
#'     \code{rioja} (\url{https://cran.r-project.org/package=rioja}).
#' 
#' @importFrom foreach `%dopar%`
#' 
#' @param modern_taxa The modern taxa abundance data, each row represents a 
#'     sampling site, each column represents a taxon.
#' @param modern_climate The modern climate value at each sampling site.
#' @param nPLS The number of components to be extracted.
#' @param trainfun Training function you want to use, either 
#'     \code{\link{WAPLS.w}} or \code{\link{TWAPLS.w}}.
#' @param predictfun Predict function you want to use: if \code{trainfun} is 
#'     \code{\link{WAPLS.w}}, then this should be \code{\link{WAPLS.predict.w}}; 
#'     if \code{trainfun} is \code{\link{TWAPLS.w}}, then this should be 
#'     \code{\link{TWAPLS.predict.w}}.
#' @param usefx Boolean flag on whether or not use \code{fx} correction.
#' @param fx The frequency of the climate value for \code{fx} correction: if 
#'     \code{usefx = FALSE}, this should be \code{NA}; otherwise, this should 
#'     be obtained from the \code{\link{fx}} function.
#' @param cpus Number of CPUs for simultaneous iterations to execute, check
#'     \code{parallel::detectCores()} for available CPUs on your machine.
#' @param test_mode boolean flag to execute the function with a limited number
#'     of iterations, \code{test_it}, for testing purposes only.
#' @param test_it number of iterations to use in the test mode.
#'
#' @return leave-one-out cross validation results
#' @export
#'
#' @examples
#' \dontrun{
#' # Load modern pollen data
#' modern_pollen <- read.csv("/path/to/modern_pollen.csv")
#'                                       
#' # Extract taxa
#' taxaColMin <- which(colnames(modern_pollen) == "taxa0")
#' taxaColMax <- which(colnames(modern_pollen) == "taxaN")
#' taxa <- modern_pollen[, taxaColMin:taxaColMax]
#' 
#' # Get the frequency of each climate variable fx
#' fx_Tmin <- fxTWAPLS::fx(modern_pollen$Tmin, bin = 0.02)
#' fx_gdd <- fxTWAPLS::fx(modern_pollen$gdd, bin = 20)
#' fx_alpha <- fxTWAPLS::fx(modern_pollen$alpha, bin = 0.002)
#' 
#' # MTCO
#' ## fx
#' fit_Tmin <- fxTWAPLS::WAPLS.w(taxa, modern_pollen$Tmin, nPLS = 5)
#' ## LOOCV
#' test_mode <- TRUE # It should be set to FALSE before running
#' ### without fx
#' cv_Tmin <- fxTWAPLS::cv.w(taxa,
#'                           modern_pollen$Tmin,
#'                           nPLS = 5,
#'                           fxTWAPLS::WAPLS.w,
#'                           fxTWAPLS::WAPLS.predict.w,
#'                           cpus = 2, # Remove the following line
#'                           test_mode = test_mode)
#' ### with fx
#' cv_f_Tmin <- fxTWAPLS::cv.w(taxa,
#'                             modern_pollen$Tmin,
#'                             nPLS = 5,
#'                             fxTWAPLS::WAPLS.w,
#'                             fxTWAPLS::WAPLS.predict.w,
#'                             usefx = TRUE,
#'                             fx = fx_Tmin,
#'                             cpus = 2, # Remove the following line
#'                             test_mode = test_mode)  
#' }
#' @seealso \code{\link{fx}}, \code{\link{TWAPLS.w}}, 
#'     \code{\link{TWAPLS.predict.w}}, \code{\link{WAPLS.w}}, and 
#'     \code{\link{WAPLS.predict.w}}
cv.w <- function(modern_taxa,
                 modern_climate,
                 nPLS = 5,
                 trainfun,
                 predictfun,
                 usefx = FALSE,
                 fx = NA,
                 cpus = 4,
                 test_mode = FALSE,
                 test_it = 5) {
  i <- NULL # Local binding
  x <- modern_climate
  y <- modern_taxa
  
  # Check the number of CPUs does not exceed the availability
  avail_cpus <- parallel::detectCores() - 1
  cpus <- ifelse(cpus > avail_cpus, avail_cpus, cpus)
  
  # Start parallel backend
  cl <- parallel::makeCluster(cpus, setup_strategy = "sequential")
  doParallel::registerDoParallel(cl)

  # Load binary operator for backend
  `%dopar%` <- foreach::`%dopar%`
  
  # Create list of indices to loop through
  idx <- seq_len(length(x))
  # Reduce the list of indices, if test_mode = TRUE
  if (test_mode) {
    idx <- 1:test_it
  }
  all.cv.out <- foreach::foreach(i = idx,
                                 .combine = rbind, #comb_pb(max(idx)),
                                 .verbose = FALSE) %dopar% {
                                   fit <- trainfun(y[-i, ], 
                                                   x[-i], 
                                                   nPLS, 
                                                   usefx, 
                                                   fx[-i])
                                   xnew <- predictfun(fit, y[i, ])[["fit"]]
                                   data.frame(x[i], xnew)
                                 }
  parallel::stopCluster(cl) # Stop cluster
  colnames(all.cv.out) <- c("test.x", paste0("comp", 1:nPLS))
  return(all.cv.out)
}

#' Get the distance between points
#' 
#' Get the distance between points, the output will be used in 
#'     \code{\link{get_pseudo}}.
#' 
#' @importFrom foreach `%dopar%` 
#' 
#' @param point Each row represents a sampling site, the first column is 
#'     longitude and the second column is latitude, both in decimal format.
#' @param cpus Number of CPUs for simultaneous iterations to execute, check
#'     \code{parallel::detectCores()} for available CPUs on your machine.
#' @param test_mode Boolean flag to execute the function with a limited number
#'     of iterations, \code{test_it}, for testing purposes only.
#' @param test_it Number of iterations to use in the test mode.
#'    
#' @return Distance matrix, the value at the \code{i-th} row, means the distance 
#'     between the \code{i-th} sampling site and the whole sampling sites.
#' @export
#' 
#' @examples
#' \dontrun{
#' # Load modern pollen data
#' modern_pollen <- read.csv("/path/to/modern_pollen.csv")
#' 
#' point <- modern_pollen[, c("Long", "Lat")]
#' test_mode <- TRUE # It should be set to FALSE before running
#' dist <- fxTWAPLS::get_distance(point, 
#'                                cpus = 2, # Remove the following line
#'                                test_mode = test_mode)
#' }
#' 
#' @seealso \code{\link{get_pseudo}}
get_distance <- function(point, cpus = 4, test_mode = FALSE, test_it = 5) {
  i <- NULL # Local binding
  colnames(point) <- c("Long", "Lat")
  tictoc::tic("Distance between points")
  
  # Check the number of CPUs does not exceed the availability
  avail_cpus <- parallel::detectCores() - 1
  cpus <- ifelse(cpus > avail_cpus, avail_cpus, cpus)
  
  # Start parallel backend
  cl <- parallel::makeCluster(cpus, setup_strategy = "sequential")
  doParallel::registerDoParallel(cl)
  
  # Load binary operator for backend
  `%dopar%` <- foreach::`%dopar%`
  
  # Create list of indices to loop through
  idx <- seq_len(nrow(point))
  # Reduce the list of indices, if test_mode = TRUE
  if (test_mode) {
    idx <- 1:test_it
  }
  dist <- foreach::foreach(i = idx,
                           .combine = rbind) %dopar% {
                             tmp <- rep(0, nrow(point))
                             lon1 <- point[i, "Long"]
                             lat1 <- point[i, "Lat"]
                             for (j in seq_len(nrow(point))) {
                               lon2 <- point[j, "Long"]
                               lat2 <- point[j, "Lat"]
                               tmp[j] <- geosphere::distm(c(lon1, lat1),
                                                          c(lon2, lat2),
                                                          fun = geosphere::distHaversine)
                             }
                             tmp
                           }
  
  parallel::stopCluster(cl) # Stop cluster
  tictoc::toc()
  return(dist)
}

#' Get geographically and climatically close sites
#' 
#' Get the sites which are both geographically and climatically close to the 
#'     test site, which could result in pseudo-replication and inflate the 
#'     cross-validation statistics. The output will be used in 
#'     \code{\link{cv.pr.w}}.
#'
#' @param dist Distance matrix which contains the distance from other sites.
#' @param x The modern climate values.
#' @param cpus Number of CPUs for simultaneous iterations to execute, check
#'     \code{parallel::detectCores()} for available CPUs on your machine.
#' @param test_mode Boolean flag to execute the function with a limited number
#'     of iterations, \code{test_it}, for testing purposes only.
#' @param test_it Number of iterations to use in the test mode.
#' 
#' @return The geographically and climatically close sites to each test site.
#' @export
#'
#' @examples
#' \dontrun{
#' # Load modern pollen data
#' modern_pollen <- read.csv("/path/to/modern_pollen.csv")
#' 
#' point <- modern_pollen[, c("Long", "Lat")]
#' test_mode <- TRUE # It should be set to FALSE before running
#' dist <- fxTWAPLS::get_distance(point, 
#'                                cpus = 2, # Remove the following line
#'                                test_mode = test_mode)
#' pseudo_Tmin <- fxTWAPLS::get_pseudo(dist, 
#'                                     modern_pollen$Tmin, 
#'                                     cpus = 2, # Remove the following line
#'                                     test_mode = test_mode)
#' }
#' @seealso \code{\link{get_distance}}
get_pseudo <- function(dist, x, cpus = 4, test_mode = FALSE, test_it = 5) {
  i <- NULL # Local binding
  # Check the number of CPUs does not exceed the availability
  avail_cpus <- parallel::detectCores() - 1
  cpus <- ifelse(cpus > avail_cpus, avail_cpus, cpus)
  
  # Start parallel backend
  cl <- parallel::makeCluster(cpus, setup_strategy = "sequential")
  doParallel::registerDoParallel(cl)
  
  # Load binary operator for backend
  `%dopar%` <- foreach::`%dopar%`
  
  # Create list of indices to loop through
  idx <- seq_len(length(x))
  # Reduce the list of indices, if test_mode = TRUE
  if (test_mode) {
    idx <- 1:test_it
  }
  pseudo <- foreach::foreach(i = idx) %dopar% {
    which(dist[i, ] < 50000 & abs(x - x[i]) < 0.02 * (max(x) - min(x)))
  }
  parallel::stopCluster(cl) # Stop cluster
  return(pseudo)
}

#' Pseudo-removed leave-out cross-validation
#'
#' @param modern_taxa The modern taxa abundance data, each row represents a 
#'     sampling site, each column represents a taxon.
#' @param modern_climate The modern climate value at each sampling site.
#' @param nPLS The number of components to be extracted.
#' @param trainfun Training function you want to use, either 
#'     \code{\link{WAPLS.w}} or \code{\link{TWAPLS.w}}.
#' @param predictfun Predict function you want to use: if \code{trainfun} is 
#'     \code{\link{WAPLS.w}}, then this should be \code{\link{WAPLS.predict.w}}; 
#'     if \code{trainfun} is \code{\link{TWAPLS.w}}, then this should be 
#'     \code{\link{TWAPLS.predict.w}}.
#' @param pseudo The geographically and climatically close sites to each test 
#'     site, obtained from \code{\link{get_pseudo}} function.
#' @param usefx Boolean flag on whether or not use \code{fx} correction.
#' @param fx The frequency of the climate value for \code{fx} correction: if 
#'     \code{usefx} is FALSE, this should be \code{NA}; otherwise, this should 
#'     be obtained from the \code{\link{fx}} function.
#' @param cpus Number of CPUs for simultaneous iterations to execute, check
#'     \code{parallel::detectCores()} for available CPUs on your machine.
#' @param test_mode Boolean flag to execute the function with a limited number
#'     of iterations, \code{test_it}, for testing purposes only.
#' @param test_it Number of iterations to use in the test mode.
#'
#' @return Leave-one-out cross validation results.
#' @export
#'
#' @examples
#' \dontrun{
#' # Load modern pollen data
#' modern_pollen <- read.csv("/path/to/modern_pollen.csv")
#'                                       
#' # Extract taxa
#' taxaColMin <- which(colnames(modern_pollen) == "taxa0")
#' taxaColMax <- which(colnames(modern_pollen) == "taxaN")
#' taxa <- modern_pollen[, taxaColMin:taxaColMax]
#' 
#' point <- modern_pollen[, c("Long", "Lat")]
#' test_mode <- TRUE # It should be set to FALSE before running
#' dist <- fxTWAPLS::get_distance(point, 
#'                                cpus = 2, # Remove the following line
#'                                test_mode = test_mode)
#' pseudo_Tmin <- fxTWAPLS::get_pseudo(dist, 
#'                                     modern_pollen$Tmin, 
#'                                     cpus = 2, # Remove the following line
#'                                     test_mode = test_mode)
#' # Test WAPLS
#' cv_pr_Tmin <- fxTWAPLS::cv.pr.w(taxa,
#'                                 modern_pollen$Tmin,
#'                                 nPLS = 5,
#'                                 fxTWAPLS::WAPLS.w,
#'                                 fxTWAPLS::WAPLS.predict.w,
#'                                 pseudo_Tmin,
#'                                 cpus = 2, # Remove the following line
#'                                 test_mode = test_mode)
#' # Test TWAPLS
#' cv_pr_Tmin2 <- fxTWAPLS::cv.pr.w(taxa,
#'                                  modern_pollen$Tmin,
#'                                  nPLS = 5,
#'                                  fxTWAPLS::TWAPLS.w,
#'                                  fxTWAPLS::TWAPLS.predict.w,
#'                                  pseudo_Tmin,
#'                                  cpus = 2, # Remove the following line
#'                                  test_mode = test_mode)
#' }
#' 
#' @seealso \code{\link{fx}}, \code{\link{TWAPLS.w}}, 
#'     \code{\link{TWAPLS.predict.w}}, \code{\link{WAPLS.w}}, and 
#'     \code{\link{WAPLS.predict.w}}
cv.pr.w <- function(modern_taxa,
                    modern_climate,
                    nPLS = 5,
                    trainfun,
                    predictfun,
                    pseudo,
                    usefx = FALSE,
                    fx = NA,
                    cpus = 4,
                    test_mode = TRUE,
                    test_it = 5) {
  i <- NULL # Local binding
  x <- modern_climate
  y <- modern_taxa
  
  # Check the number of CPUs does not exceed the availability
  avail_cpus <- parallel::detectCores() - 1
  cpus <- ifelse(cpus > avail_cpus, avail_cpus, cpus)
  
  # Start parallel backend
  cl <- parallel::makeCluster(cpus, setup_strategy = "sequential")
  doParallel::registerDoParallel(cl)
  
  # Load binary operator for backend
  `%dopar%` <- foreach::`%dopar%`
  
  # Create list of indices to loop through
  idx <- seq_len(length(x))
  # Reduce the list of indices, if test_mode = TRUE
  if (test_mode) {
    idx <- 1:test_it
  }
  all.cv.out <- foreach::foreach(i = idx,
                                 .combine = rbind) %dopar% {
                                   leave <- unlist(pseudo[i])
                                   fit <- trainfun(y[-leave, ], 
                                                   x[-leave], 
                                                   nPLS, 
                                                   usefx, 
                                                   fx[-leave])
                                   xnew <- predictfun(fit, y[i, ])[["fit"]]
                                   data.frame(x[i], xnew)
                                 }
  parallel::stopCluster(cl) # Stop cluster
  
  # assign column names to all.cv.out
  colnames(all.cv.out) <- c("test.x", paste0("comp", 1:nPLS))
  
  return(all.cv.out)
}

#' Random t-test
#' 
#' Do a random t-test to the cross-validation results.
#' 
#' @importFrom stats cor lm rbinom
#' 
#' @param cvoutput Cross-validation output either from \code{\link{cv.w}} or 
#'     \code{\link{cv.pr.w}}.
#' @param n.perm The number of permutation times to get the p value, which 
#'     assesses whether using the current number of components is significantly 
#'     different from using one less.
#'
#' @return A matrix of the statistics of the cross-validation results. Each 
#'     component is described below:
#'     \describe{
#'     \item{\code{R2}}{the coefficient of determination (the larger, the 
#'     better the fit).}
#'     \item{\code{Avg.Bias}}{average bias.}
#'     \item{\code{Max.Bias}}{maximum bias.}
#'     \item{\code{Min.Bias}}{minimum bias.}
#'     \item{\code{RMSEP}}{root-mean-square error of prediction (the smaller, 
#'     the better the fit).}
#'     \item{\code{delta.RMSEP}}{the percent change of RMSEP using the current 
#'     number of components than using one component less.}
#'     \item{\code{p}}{assesses whether using the current number of components 
#'     is significantly different from using one component less, which is used 
#'     to choose the last significant number of components to avoid 
#'     over-fitting.}
#'     \item{\code{-}}{The degree of overall compression is assessed by doing 
#'     linear regression to the cross-validation result and the observed 
#'     climate values.
#'         \itemize{
#'         \item \code{Compre.b0}: the intercept.
#'         \item \code{Compre.b1}: the slope (the closer to 1, the less the 
#'         overall compression).
#'         \item \code{Compre.b0.se}: the standard error of the intercept.
#'         \item \code{Compre.b1.se}: the standard error of the slope.
#'         }
#'     }
#'     }
#'         
#' @export
#'
#' @examples
#' \dontrun{
#' # Load modern pollen data
#' modern_pollen <- read.csv("/path/to/modern_pollen.csv")
#'                                       
#' # Extract taxa
#' taxaColMin <- which(colnames(modern_pollen) == "taxa0")
#' taxaColMax <- which(colnames(modern_pollen) == "taxaN")
#' taxa <- modern_pollen[, taxaColMin:taxaColMax]
#' 
#' # Get the frequency of each climate variable fx
#' fx_Tmin <- fxTWAPLS::fx(modern_pollen$Tmin, bin = 0.02)
#' fx_gdd <- fxTWAPLS::fx(modern_pollen$gdd, bin = 20)
#' fx_alpha <- fxTWAPLS::fx(modern_pollen$alpha, bin = 0.002)
#' 
#' # MTCO
#' ## fx
#' fit_Tmin <- fxTWAPLS::WAPLS.w(taxa, modern_pollen$Tmin, nPLS = 5)
#' ## LOOCV
#' test_mode <- TRUE # It should be set to FALSE before running
#' ### without fx
#' cv_Tmin <- fxTWAPLS::cv.w(taxa,
#'                           modern_pollen$Tmin,
#'                           nPLS = 5,
#'                           fxTWAPLS::WAPLS.w,
#'                           fxTWAPLS::WAPLS.predict.w,
#'                           cpus = 2, # Remove the following line
#'                           test_mode = test_mode)
#' ### with fx
#' cv_f_Tmin <- fxTWAPLS::cv.w(taxa,
#'                             modern_pollen$Tmin,
#'                             nPLS = 5,
#'                             fxTWAPLS::WAPLS.w,
#'                             fxTWAPLS::WAPLS.predict.w,
#'                             usefx = TRUE,
#'                             fx = fx_Tmin,
#'                             cpus = 2, # Remove the following line
#'                             test_mode = test_mode)
#'                             
#' ## Random t-test
#' rand_Tmin <- fxTWAPLS::rand.t.test.w(cv_Tmin, n.perm = 999)
#' rand_f_Tmin <- fxTWAPLS::rand.t.test.w(cv_f_Tmin, n.perm = 999)
#' }
#' 
#' @seealso \code{\link{cv.w}} and \code{\link{cv.pr.w}}
rand.t.test.w <- function(cvoutput, n.perm = 999) {
  ncomp <- ncol(cvoutput) - 1
  output <- matrix(NA, ncomp, 11)
  colnames(output) <- c("R2",
                        "Avg.Bias",
                        "Max.Bias",
                        "Min.Bias",
                        "RMSEP",
                        "delta.RMSEP",
                        "p",
                        "Compre.b0",
                        "Compre.b1",
                        "Compre.b0.se",
                        "Compre.b1.se")
  
  for (i in 1:ncomp) {
    cv.x <- cvoutput[, 1]
    cv.i <- cvoutput[, 1 + i]
    output[i, "RMSEP"] <- sqrt(mean((cv.i - cv.x) ^ 2))
    output[i, "R2"] <- cor(cv.i, cv.x) ^ 2
    output[i, "Avg.Bias"] <- mean(cv.i - cv.x)
    output[i, "Max.Bias"] <- max(abs(cv.i - cv.x))
    output[i, "Min.Bias"] <- min(abs(cv.i - cv.x))
    output[i, c("Compre.b0", "Compre.b1")] <-
      summary(lm(cv.i ~ cv.x))[["coefficients"]][, "Estimate"]
    output[i, c("Compre.b0.se", "Compre.b1.se")] <-
      summary(lm(cv.i ~ cv.x))[["coefficients"]][, "Std. Error"]
  }
  # get delta.RMSEP
  for(i in 1:ncomp) {
    if (i == 1) {
      rmsep.null <- sqrt(mean((cv.x - mean(cv.x)) ^ 2))
      output[i, "delta.RMSEP"] <-
        (output[i, "RMSEP"] - rmsep.null) * 100 / rmsep.null
    } else {
      output[i, "delta.RMSEP"] <-
        (output[i, "RMSEP"] - output[i - 1, "RMSEP"]) * 100 / output[i - 1, "RMSEP"]
    }
  }
  # get p-value, which describes whether using the number of components now has 
  # a significant difference than using one less
  e0 <- cv.x - mean(cv.x)
  e <- cbind(e0, cvoutput[, 2:ncol(cvoutput)] - cv.x)
  
  t.res <- vector("numeric", ncomp)
  t <- vector("numeric", n.perm + 1)
  t.res[] <- NA
  n <- nrow(e)
  for (i in 1:ncomp) {
    d <- e[, i] ^ 2 - e[, i + 1] ^ 2
    t[1] <- mean(d, na.rm = TRUE)
    for (j in 1:n.perm) {
      sig <- 2 * rbinom(n, 1, 0.5) - 1
      t[j + 1] <- mean(d * sig, na.rm = TRUE)
    }
    t.res[i] <- sum(t >= t[1]) / (n.perm + 1)
  }
  output[, "p"] <- t.res
  
  print(output)
  return(output)
}

#' Plot the training results
#' 
#' Plot the training results, the black line is the 1:1 line, the red line is 
#'     the linear regression line to fitted and \code{x}, which shows the degree 
#'     of overall compression.
#' 
#' @param train_output Training output, can be the output of WA-PLS, WA-PLS with 
#'     \code{fx} correction, TWA-PLS, or TWA-PLS with \code{fx} correction.
#' @param col Choose which column of the fitted value to plot, in other words, 
#'     how many number of components you want to use.
#'     
#' @export
#' 
#' @return Plotting status.
#' 
#' @examples
#' \dontrun{
#' # Load modern pollen data
#' modern_pollen <- read.csv("/path/to/modern_pollen.csv")
#'                                       
#' # Extract taxa
#' taxaColMin <- which(colnames(modern_pollen) == "taxa0")
#' taxaColMax <- which(colnames(modern_pollen) == "taxaN")
#' taxa <- modern_pollen[, taxaColMin:taxaColMax]
#' 
#' # Get the frequency of each climate variable fx
#' fx_Tmin <- fxTWAPLS::fx(modern_pollen$Tmin, bin = 0.02)
#' fx_gdd <- fxTWAPLS::fx(modern_pollen$gdd, bin = 20)
#' fx_alpha <- fxTWAPLS::fx(modern_pollen$alpha, bin = 0.002)
#' 
#' # MTCO
#' ## WAPLS
#' fit_Tmin <- fxTWAPLS::WAPLS.w(taxa, modern_pollen$Tmin, nPLS = 5)
#' fit_f_Tmin <- fxTWAPLS::WAPLS.w(taxa, 
#'                                 modern_pollen$Tmin, 
#'                                 nPLS = 5, 
#'                                 usefx = TRUE, 
#'                                 fx = fx_Tmin)
#' ## TWAPLS
#' fit_t_Tmin <- fxTWAPLS::TWAPLS.w(taxa, modern_pollen$Tmin, nPLS = 5)
#' fit_tf_Tmin <- fxTWAPLS::TWAPLS.w(taxa, 
#'                                   modern_pollen$Tmin, 
#'                                   nPLS = 5, 
#'                                   usefx = TRUE, 
#'                                   fx = fx_Tmin)
#' fxTWAPLS::plot_train(fit_Tmin, 3)
#' fxTWAPLS::plot_train(fit_f_Tmin, 3)
#' fxTWAPLS::plot_train(fit_t_Tmin, 3)
#' fxTWAPLS::plot_train(fit_tf_Tmin, 3)
#' }
#' 
#' @seealso \code{\link{TWAPLS.w}} and \code{\link{WAPLS.w}}
plot_train <- function(train_output, col) {
  x <- train_output[["x"]]
  fitted <- train_output[["fit"]][, col]
  plotdata <- cbind.data.frame(x, fitted)
  
  max <- max(fitted, x)
  min <- min(fitted, x)

  # plot the fitted curve, the black line is the 1:1 line, the red line is the 
  # linear regression line to fitted and x, which shows the overall compression
  ggplot2::ggplot(plotdata, ggplot2::aes(x, fitted)) + 
    ggplot2::geom_point(size = 0.4) + ggplot2::theme_bw() +
    ggplot2::geom_abline(slope = 1, intercept = 0) + 
    ggplot2::xlim(min, max) + ggplot2::ylim(min, max) +
    ggplot2::geom_smooth(method = 'lm',
                         formula = y ~ x,
                         color = 'red')
  return(TRUE)
}

#' Plot the residuals
#' 
#' Plot the residuals, the black line is 0 line, the red line is the locally 
#'     estimated scatterplot smoothing, which shows the degree of local 
#'     compression.
#' 
#' @param train_output Training output, can be the output of WA-PLS, WA-PLS with 
#'     \code{fx} correction, TWA-PLS, or TWA-PLS with \code{fx} correction
#' @param col Choose which column of the fitted value to plot, in other words, 
#'     how many number of components you want to use.
#' 
#' @export
#' 
#' @return Plotting status.
#' 
#' @examples
#' \dontrun{
#' # Load modern pollen data
#' modern_pollen <- read.csv("/path/to/modern_pollen.csv")
#'                                       
#' # Extract taxa
#' taxaColMin <- which(colnames(modern_pollen) == "taxa0")
#' taxaColMax <- which(colnames(modern_pollen) == "taxaN")
#' taxa <- modern_pollen[, taxaColMin:taxaColMax]
#' 
#' # Get the frequency of each climate variable fx
#' fx_Tmin <- fxTWAPLS::fx(modern_pollen$Tmin, bin = 0.02)
#' fx_gdd <- fxTWAPLS::fx(modern_pollen$gdd, bin = 20)
#' fx_alpha <- fxTWAPLS::fx(modern_pollen$alpha, bin = 0.002)
#' 
#' # MTCO
#' ## WAPLS
#' fit_Tmin <- fxTWAPLS::WAPLS.w(taxa, modern_pollen$Tmin, nPLS = 5)
#' fit_f_Tmin <- fxTWAPLS::WAPLS.w(taxa, 
#'                                 modern_pollen$Tmin, 
#'                                 nPLS = 5, 
#'                                 usefx = TRUE, 
#'                                 fx = fx_Tmin)
#' ## TWAPLS
#' fit_t_Tmin <- fxTWAPLS::TWAPLS.w(taxa, modern_pollen$Tmin, nPLS = 5)
#' fit_tf_Tmin <- fxTWAPLS::TWAPLS.w(taxa, 
#'                                   modern_pollen$Tmin, 
#'                                   nPLS = 5, 
#'                                   usefx = TRUE, 
#'                                   fx = fx_Tmin)
#' fxTWAPLS::plot_residuals(fit_Tmin, 3)
#' fxTWAPLS::plot_residuals(fit_f_Tmin, 3)
#' fxTWAPLS::plot_residuals(fit_t_Tmin, 3)
#' fxTWAPLS::plot_residuals(fit_tf_Tmin, 3)
#' }
#' 
#' @seealso \code{\link{TWAPLS.w}} and \code{\link{WAPLS.w}}
plot_residuals <- function(train_output, col) {
  x <- train_output[["x"]]
  residuals <- train_output[["fit"]][, col] - train_output[["x"]]
  plotdata <- cbind.data.frame(x, residuals)
  
  maxr <- max(abs(residuals))

  # plot the residuals, the black line is 0 line, the red line is the locally 
  # estimated scatterplot smoothing, which shows the degree of local compression
  ggplot2::ggplot(plotdata, ggplot2::aes(x, residuals)) + 
    ggplot2::geom_point(size = 0.4) + ggplot2::theme_bw() +
    ggplot2::geom_abline(slope = 0, intercept = 0) + 
    ggplot2::xlim(min(x), max(x)) + ggplot2::ylim(-maxr, maxr) +
    ggplot2::geom_smooth(method = 'loess',
                         color = 'red',
                         se = FALSE)
  return(TRUE)
}