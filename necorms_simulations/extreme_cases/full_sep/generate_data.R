#---------
# generate all the datasets
#---------

dgp1 <- function(n, t){
  means_temp <- rep(c(0,3,6,10,15,20), each = n)
  Y <- rnorm(n * 6, sd = sqrt(0.6)) + means_temp
  true_part <- as.numeric(as.factor(means_temp))
  list(Y, true_part, rep(0:5, each = n))
}

true_dens <- function(grid, t){
  rbind(
    dnorm(grid, 0 / t, sqrt(0.6)),
    dnorm(grid, 3 / t, sqrt(0.6)),
    dnorm(grid, 6 / t, sqrt(0.6)),
    dnorm(grid, 10 / t, sqrt(0.6)),
    dnorm(grid, 15 / t, sqrt(0.6)),
    dnorm(grid, 20 / t, sqrt(0.6))
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
  write.table(x = cbind(group, Y), file = paste0("data/rep_", i, ".csv"), row.names = F, col.names = F, sep = ",")
  write.table(x = cbind(part_true), file = paste0("true_part/rep_", i, ".csv"), sep = ",", row.names = F, col.names = F)
}

#--------

save(data_list, file="data_list.RData")

