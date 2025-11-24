# read the data
load(file="data_list.RData")
niter <- 15000
nburn <- 10000
grid <- seq(-5, 25, length.out = 250)

#--------------------
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

true_dens <- function(grid, g, t){
  rbind(
    matrix(rep(0.5 * dnorm(grid, 0 / t, sqrt(0.6)) + 0.25 * dnorm(grid, 3 / t, sqrt(0.6)) + 0.25 * dnorm(grid, 10 / t, sqrt(0.6)), g), 
           ncol = length(grid), byrow = T),
    matrix(rep(0.5 * dnorm(grid, 0 / t, sqrt(0.6)) + 0.5 * dnorm(grid, 6 / t, sqrt(0.6)), g), 
           ncol = length(grid), byrow = T), 
    matrix(rep(0.4 * dnorm(grid, 0 / t, sqrt(0.6)) + 0.2 * dnorm(grid, 3 / t, sqrt(0.6)) + 0.2 * dnorm(grid, 10 / t, sqrt(0.6)) + 0.2 * dnorm(grid, 15 / t, sqrt(0.6)), g), 
           ncol = length(grid), byrow = T))
}

scaling_factor <- 1
  
#--------------------
nrep <- length(data_list)

#--------------------
writeLines(c(""), "output_high_group.txt")
Rcpp::sourceCpp('utilities.cpp')
Rcpp::sourceCpp('main_nested_CPP.cpp')
Rcpp::sourceCpp('main_CAM_CPP.cpp')
out_list_100_HG_HS <- list()

set.seed(42)
for(i in 1:100){
  res <- c()
  
  #-------- NESTED
  est_nested <- main_nested(Y = data_list[[i]][[1]], group = data_list[[i]][[2]], d = 30, niter = niter, 
                            nburn = nburn, thin = 1,  m0 = 0, k0 = 0.01, a0 = 1.1, b0 = 0.1, m1 = 0, s21 = 10, 
                            tau1 = .1, tau2 = .1, a1 = .1, b1 = .1, sigma = 0.2, phi = 0.5, a_sigma = 2, 
                            b_sigma = 2, grid = grid, M = 100, tol = 0.005, max_iter = 50, epsilon = 0.005, 
                            eval_density = T, IS_nval = 1000, nupd = 500, MH_var = 0.01, beta = 1)
  opt_part  <- est_nested$clust[which.min(VI_LB(clean_partition(est_nested$clust), 
                                                psm_mat = psm(est_nested$clust))),] + 1
  res[1] <- VI(opt_part, data_list[[i]][[3]])
  res[2] <- mean(
    apply(cbind(apply(est_nested$dens, c(1,2), mean),true_dens(grid, 10, scaling_factor)), 1, function(x) 
      (sum(x[1:(length(grid))] * log(x[1:(length(grid))] / x[(length(grid) + 1):(2 * length(grid))])) + 
         sum(x[(length(grid) + 1):(2 * length(grid))] * log(x[(length(grid) + 1):(2 * length(grid))] / x[1:(length(grid))]))) / 2)
  )
  opt_part_nested  <- est_nested$nested[which.min(VI_LB(clean_partition(est_nested$nested), 
                                                        psm_mat = psm(est_nested$nested))),] + 1
  res[3] <- VI(opt_part_nested, rep(1:3, each = 10))
  
  #--------- CAM
  est_CAM <- main_CAM(Y = data_list[[i]][[1]], group = data_list[[i]][[2]], J = 30, niter = niter, 
                      nburn = nburn, thin = 1,  m0 = 0, k0 = 0.01, a0 = 1.1, b0 = 0.1, 
                      p_alpha_1 = .1, p_alpha_2 = .1, p_beta_1 = .1, p_beta_2 = .1, kappa = 0.5, 
                      grid = grid, eval_density = T, nupd = 500)
  opt_part  <- est_CAM$clust[which.min(VI_LB(clean_partition(est_CAM$clust), 
                                             psm_mat = psm(est_CAM$clust))),] + 1
  res[4] <- VI(opt_part, data_list[[i]][[3]])
  res[5] <- mean(
    apply(cbind(apply(est_CAM$dens, c(1,2), mean),true_dens(grid, 10, scaling_factor)), 1, function(x) 
      (sum(x[1:(length(grid))] * log(x[1:(length(grid))] / x[(length(grid) + 1):(2 * length(grid))])) + 
         sum(x[(length(grid) + 1):(2 * length(grid))] * log(x[(length(grid) + 1):(2 * length(grid))] / x[1:(length(grid))]))) / 2)
  )
  opt_part_nested  <- est_CAM$nested[which.min(VI_LB(clean_partition(est_CAM$nested), 
                                                     psm_mat = psm(est_CAM$nested))),] + 1
  res[6] <- VI(opt_part_nested, rep(1:3, each = 10))
  
  cat(paste("completed: ", i / 100, "%\n"), file = "output_high_group.txt", append = TRUE)
  out_list_100_HG_HS[[i]] <- res
}
save.image(file = "simu_nested_high_group_high_sep.Rdata")