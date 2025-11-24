library(mvtnorm)

#---------
# generate all the datasets
#---------

dgp1 <- function(n, t, d){
  
  means <- c(0, 3, 6, 10, 15)
  vcovs <- array(0, dim = c(d, d, 5))
  vcovs[,,1] <- vcovs[,,2] <- diag(0.6, d)
  vcovs[,,3] <- vcovs[,,4] <- vcovs[,,5] <- 0.6 * (matrix(0.75, d, d) + diag(0.25, d))
  
  idx_temp <- c(sample(c(1, 2, 4), size = n, prob = c(0.5, 0.25, 0.25), replace = TRUE),
                  sample(c(1, 3), size = n, prob = c(0.5, 0.5), replace = TRUE),
                  sample(c(1, 2, 4), size = n, prob = c(0.5, 0.25, 0.25), replace = TRUE),
                  sample(c(1, 3), size = n, prob = c(0.5, 0.5), replace = TRUE),
                  sample(c(1, 2, 4, 5), size = n, prob = c(0.4, 0.2, 0.2, 0.2), replace = TRUE),
                  sample(c(1, 2, 4), size = n, prob = c(0.5, 0.25, 0.25), replace = TRUE))
  Y <- t(sapply(idx_temp, function(x) rmvnorm(1, mean = rep(means[x], d) / t, sigma = vcovs[,,x])))
  list(Y, idx_temp, rep(0:5, each = n))
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

#--------

set.seed(42)
data_list <- list()
nrep <- 100
for(i in 1:nrep){
  temp <- dgp1(100, 1, 5)
  Y <- temp[[1]]
  part_true <- temp[[2]]
  group <- temp[[3]]
  
  #---------
  
  data_list[[i]] <- list(Y, group, part_true)
  
  #---------
  write.table(x = cbind(group, Y), file = paste0("data/rep_", i, ".csv"), row.names = F, col.names = F, sep = ",")
  write.table(x = cbind(part_true), file = paste0("true_part/rep_", i, ".csv"), sep = ",", row.names = F, col.names = F)
}

#--------

save(data_list, file="data_list.RData")

