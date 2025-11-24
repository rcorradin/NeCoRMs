# read the data
library(mvtnorm)
dims <- c(2, 5, 10)
file_name <- c("data_list_d2.RData", "data_list_d5.RData", "data_list_d10.RData")

# initialize
niter <- 15000
nburn <- 10000

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

true_dens <- function(grid, t, d){
  
  means <- t(sapply(c(0, 3, 6, 10, 15), function(x) rep(x, d)))
  vcovs <- array(0, dim = c(d, d, 5))
  vcovs[,,1] <- vcovs[,,2] <- diag(0.6, d)
  vcovs[,,3] <- vcovs[,,4] <- vcovs[,,5] <- 0.6 * (matrix(0.75, d, d) + diag(0.25, d))
  probs <- rbind(c(0.5, 0.25, 0, 0.25, 0), 
                 c(0.5, 0, 0.5, 0, 0), 
                 c(0.5, 0.25, 0, 0.25, 0), 
                 c(0.5, 0, 0.5, 0, 0), 
                 c(0.4, 0.2, 0, 0.2, 0.2), 
                 c(0.5, 0.25, 0, 0.25, 0))
  dens_out <- matrix(0, ncol = nrow(grid), nrow = 6)
  for(i in 1:6){
    for(j in 1:5){
      dens_out[i,] = dens_out[i,] + probs[i,j] * dmvnorm(grid, means[j,] / t, vcovs[,,j])
    }
  }
  return(dens_out)
}


scaling_factor <- 1
  
#--------------------
nrep <- length(data_list)

#--------------------
writeLines(c(""), "output_low_sep.txt")
Rcpp::sourceCpp('utilities.cpp')
Rcpp::sourceCpp('main_nested_CPP.cpp')
out_list_diff_dim <- list()

set.seed(42)
for(j in 1:3){
  out_list_diff_dim[[j]] <- list()
  load(file = file_name[j])
  for(i in 1:nrep){
    res <- c()
    d <- dims[j]
    est_nested <- main_nested(Y = data_list[[i]][[1]], group = data_list[[i]][[2]], d = 6, niter = niter, 
                              nburn = nburn, thin = 1,  m0 = rep(0, d), k0 = 0.01, S0 = diag(1, d), n0 = d + 2,
                              sigma = 0.2, phi = 0.5, a_sigma = 2, 
                              b_sigma = 2, grid = data_list[[i]][[1]], M = 100, tol = 0.005, max_iter = 50, epsilon = 0.005, 
                              eval_density = T, IS_nval = 1000, nupd = 500, MH_var = 0.01, beta = 1)
    
    opt_part  <- est_nested$clust[which.min(VI_LB(clean_partition(est_nested$clust), 
                                                  psm_mat = psm(est_nested$clust))),] + 1
    res[1] <- VI(opt_part, data_list[[i]][[3]])
    res[2] <- mean(
      apply(cbind(apply(est_nested$dens, c(1,2), mean),true_dens(data_list[[i]][[1]], scaling_factor, d)), 1, function(x) 
        (sum(x[1:(nrow(data_list[[i]][[1]]))] * log(x[1:(nrow(data_list[[i]][[1]]))] / x[(nrow(data_list[[i]][[1]]) + 1):(2 * nrow(data_list[[i]][[1]]))])) + 
           sum(x[(nrow(data_list[[i]][[1]]) + 1):(2 * nrow(data_list[[i]][[1]]))] * log(x[(nrow(data_list[[i]][[1]]) + 1):(2 * nrow(data_list[[i]][[1]]))] / x[1:(nrow(data_list[[i]][[1]]))]))) / 2)
    )
    opt_part_nested  <- est_nested$nested[which.min(VI_LB(clean_partition(est_nested$nested), 
                                                          psm_mat = psm(est_nested$nested))),] + 1
    res[3] <- VI(opt_part_nested, c(1,2,1,2,3,1))
    res[4] <- est_nested$time
    
    cat(paste("completed: ", j, " - ", i, " %\n"), file = "output_low_sep.txt", append = TRUE)
    out_list_diff_dim[[j]][[i]] <- res
  }
}
save.image(file = "simu_nested_dims.Rdata")
