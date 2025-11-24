#---------
# generate all the datasets
#---------

dgp1 <- function(n, t){
  means_temp <- c(sample(c(0, 3, 10) / t, size = n, prob = c(0.5, 0.25, 0.25), replace = TRUE),
                  sample(c(0, 6) / t, size = n, prob = c(0.5, 0.5), replace = TRUE),
                  sample(c(0, 3, 10) / t, size = n, prob = c(0.5, 0.25, 0.25), replace = TRUE),
                  sample(c(0, 6) / t, size = n, prob = c(0.5, 0.5), replace = TRUE),
                  sample(c(0, 3, 10, 15) / t, size = n, prob = c(0.4, 0.2, 0.2, 0.2), replace = TRUE),
                  sample(c(0, 3, 10) / t, size = n, prob = c(0.5, 0.25, 0.25), replace = TRUE))
  Y <- rnorm(n * 6, sd = sqrt(0.6)) + means_temp
  true_part <- as.numeric(as.factor(means_temp))
  list(Y, true_part, rep(0:5, each = n))
}

true_dens <- function(grid, t){
  rbind(
    0.5 * dnorm(grid, 0 / t, sqrt(0.6)) + 0.25 * dnorm(grid, 3 / t, sqrt(0.6)) + 0.25 * dnorm(grid, 10 / t, sqrt(0.6)),
    0.5 * dnorm(grid, 0 / t, sqrt(0.6)) + 0.5 * dnorm(grid, 6 / t, sqrt(0.6)),
    0.5 * dnorm(grid, 0 / t, sqrt(0.6)) + 0.25 * dnorm(grid, 3 / t, sqrt(0.6)) + 0.25 * dnorm(grid, 10 / t, sqrt(0.6)),
    0.5 * dnorm(grid, 0 / t, sqrt(0.6)) + 0.5 * dnorm(grid, 6 / t, sqrt(0.6)),
    0.4 * dnorm(grid, 0 / t, sqrt(0.6)) + 0.2 * dnorm(grid, 3 / t, sqrt(0.6)) + 0.2 * dnorm(grid, 10 / t, sqrt(0.6)) + 0.2 * dnorm(grid, 15 / t, sqrt(0.6)),
    0.5 * dnorm(grid, 0 / t, sqrt(0.6)) + 0.25 * dnorm(grid, 3 / t, sqrt(0.6)) + 0.25 * dnorm(grid, 10 / t, sqrt(0.6))
  )
}

#--------

set.seed(42)
data_list <- list()
nrep <- 100
for(i in 1:nrep){
  temp <- dgp1(100, 1)
  Y <- temp[[1]]
  part_true <- temp[[2]]
  group <- temp[[3]]
  
  #---------
  
  data_list[[i]] <- list(Y, group, part_true)
  
  #---------
  write.table(x = cbind(group, Y), file = paste0("/home/riccardocorradin/Dropbox/PROJECTS/Andrea&Federico/CoRMs/01Code/simulations/nested/low_dim_high_sep/data/rep_", i, ".csv"), row.names = F, col.names = F, sep = ",")
  write.table(x = cbind(part_true), file = paste0("/home/riccardocorradin/Dropbox/PROJECTS/Andrea&Federico/CoRMs/01Code/simulations/nested/low_dim_high_sep/true_part/rep_", i, ".csv"), sep = ",", row.names = F, col.names = F)
}

#--------

save(data_list, file="data_list.RData")

