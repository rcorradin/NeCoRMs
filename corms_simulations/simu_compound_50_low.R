library(ggplot2)
Rcpp::sourceCpp('utilities.cpp')
Rcpp::sourceCpp('main_compound_CPP.cpp')

#-----------------------------------------------
dgp1 <- function(n, t){
  means_temp <- c(sample(c(6, 10, 15) / t, size = n, prob = c(0.5, 0.25, 0.25), replace = TRUE),
                  sample(c(10, 15, 20) / t, size = n, prob = c(0.25, 0.5, 0.25), replace = TRUE),
                  sample(c(6, 10) / t, size = n, prob = c(0.5, 0.5), replace = TRUE),
                  sample(c(0, 3, 15, 20) / t, size = n, prob = c(0.2, 0.4, 0.2, 0.2), replace = TRUE),
                  sample(c(0, 15) / t, size = n, prob = c(0.5, 0.5), replace = TRUE),
                  sample(c(0, 6, 10, 15) / t, size = n, prob = c(0.25, 0.25, 0.25, 0.25), replace = TRUE))
  Y <- rnorm(n * 6, sd = sqrt(0.6)) + means_temp
  true_part <- as.numeric(as.factor(means_temp))
  list(Y, true_part, rep(0:5, each = n))
}

true_dens <- function(grid, t){
  rbind(
    0.5 * dnorm(grid, 6 / t, sqrt(0.6)) + 0.25 * dnorm(grid, 10 / t, sqrt(0.6)) + 0.25 * dnorm(grid, 15 / t, sqrt(0.6)),
    0.25 * dnorm(grid, 10 / t, sqrt(0.6)) + 0.5 * dnorm(grid, 15 / t, sqrt(0.6)) + 0.25 * dnorm(grid, 20 / t, sqrt(0.6)),
    0.5 * dnorm(grid, 6 / t, sqrt(0.6)) + 0.5 * dnorm(grid, 10 / t, sqrt(0.6)),
    0.2 * dnorm(grid, 0 / t, sqrt(0.6)) + 0.4 * dnorm(grid, 3 / t, sqrt(0.6)) + 0.2 * dnorm(grid, 15 / t, sqrt(0.6)) + 0.2 * dnorm(grid, 20 / t, sqrt(0.6)),
    0.5 * dnorm(grid, 0 / t, sqrt(0.6)) + 0.5 * dnorm(grid, 15 / t, sqrt(0.6)),
    0.25 * dnorm(grid, 0 / t, sqrt(0.6)) + 0.25 * dnorm(grid, 6 / t, sqrt(0.6)) + 0.25 * dnorm(grid, 10 / t, sqrt(0.6)) + 0.25 * dnorm(grid, 15 / t, sqrt(0.6))
  )
}


VI <- function(c1, c2){
  n <- length(c1)
  
  T1 <- table(c1) / n
  T2 <- table(c2) / n
  
  H1 <- - sum(T1 * log(T1))
  H2 <- - sum(T2 * log(T2))
  
  T12 <- table(c1, c2) / n
  I12 <- sum(T12 * log(T12 / T1 %*% t(T2)), na.rm = T)
  
  return((H1 + H2 - 2 * I12) / log(n))
}

#-----------------------------------------------
sample_size <- 50
scaling_factor <- 1.5
nrep <- 100
niter <- 15000
nburn <- 10000
grid <- seq(-5, 25, length.out = 250)
writeLines(c(""), "output_50_low.txt")
out_list_50_low <- list()

set.seed(42)
for(i in 1:nrep){
  dg <- dgp1(sample_size, scaling_factor)
  res <- c()
  mod <- main_compound(Y = dg[[1]], group = dg[[3]], d = 6, niter = niter, nburn = nburn, thin = 1,
                       m0 = 0, k0 = 0.01, a0 = 1.1, b0 = 0.1, m1 = 0, s21 = 10, tau1 = .1, tau2 = .1,
                       a1 = .1, b1 = .1, sigma = 0.2, phi = 0.5, a_sigma = 2, b_sigma = 2, grid = grid, 
                       M = 100, tol = 0.0001, max_iter = 50, epsilon = 0.005, eval_density = T, 
                       IS_nval = 1000, nupd = 50, MH_var = 0.01)
  opt_part  <- mod$clust[which.min(VI_LB(clean_partition(mod$clust), psm_mat = psm(mod$clust))),] + 1
  res[1] <- VI(opt_part, dg[[2]])
  res[2] <- mean(
    apply(cbind(apply(mod$dens, c(1,2), mean),true_dens(grid, scaling_factor)), 1, function(x) 
      (sum(x[1:(length(grid))] * log(x[1:(length(grid))] / x[(length(grid) + 1):(2 * length(grid))])) + 
        sum(x[(length(grid) + 1):(2 * length(grid))] * log(x[(length(grid) + 1):(2 * length(grid))] / x[1:(length(grid))]))) / 2)
  )
  cat(paste("completed: ", rep / 100, "%\n"), file = "output_50_low.txt", append = TRUE)
  out_list_50_low[[i]] <- res
}
save.image(file = "simu_compound_50_low.Rdata")

