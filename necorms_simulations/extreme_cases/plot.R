load(file = "necorms_simulations/extreme_cases/full_sep/simu_nested_full_sep.Rdata")
out_mat <- matrix(0, ncol = 5, nrow = 400)
out_temp_VI <- rep(0, 400)
out_temp_VI_nest <- rep(0, 400)
out_temp_sKL <- rep(0, 400)

for(i in 1:100){
  out_temp_VI[i] <- out_list_full[[i]][1]
  out_temp_sKL[i] <- out_list_full[[i]][2]
  out_temp_VI_nest[i] <- out_list_full[[i]][3]
  out_mat[i,4] <- "FS"
  out_mat[i,5] <- "nCoRMs"
  
  out_temp_VI[i+100] <- out_list_full[[i]][4]
  out_temp_sKL[i+100] <- out_list_full[[i]][5]
  out_temp_VI_nest[i+100] <- out_list_full[[i]][6]
  out_mat[i+100,4] <- "FS"
  out_mat[i+100,5] <- "CAM"
  
}

load(file = "necorms_simulations/extreme_cases/no_sep/simu_nested_no_sep.Rdata")
for(i in 1:100){
  out_temp_VI[i+200] <- out_list_no[[i]][1]
  out_temp_sKL[i+200] <- out_list_no[[i]][2] 
  out_temp_VI_nest[i+200] <- out_list_no[[i]][3]
  out_mat[i+200,4] <- "NS"
  out_mat[i+200,5] <- "nCoRMs"
  
  out_temp_VI[i+300] <- out_list_no[[i]][4]
  out_temp_sKL[i+300] <- out_list_no[[i]][5] 
  out_temp_VI_nest[i+300] <- out_list_no[[i]][6]
  out_mat[i+300,4] <- "NS"
  out_mat[i+300,5] <- "CAM"
}
out_mat <- as.data.frame(out_mat)
out_mat[,1] <- out_temp_VI
out_mat[,2] <- out_temp_sKL
out_mat[,3] <- out_temp_VI_nest

colnames(out_mat) <- c("VI", "sKL", "VInest", "scenario", "model")
out_mat$model <- factor(out_mat$model, levels = c("nCoRMs", "CAM"))
out_mat$scenario <- factor(out_mat$scenario, levels = c("FS", "NS"))

#-------------------------------

library(ggplot2)
library(ggpubr)
blue <- "#100f7a"
red <- "#7a0f10"
dgray <- "#999999"

p1 <- ggplot(out_mat) + 
  geom_boxplot(mapping = aes(y = VI, x = scenario, fill = model), alpha = 0.5) +
  theme_classic() +
  ylab("variation of information") + 
  scale_fill_manual(values = c(dgray, blue))
p2 <- ggplot(out_mat) + 
  geom_boxplot(mapping = aes(y = sKL, x = scenario, fill = model), alpha = 0.5)+
  theme_classic() + 
  ylab("J-divergence") + 
  scale_fill_manual(values = c(dgray, blue)) +
  scale_y_log10()

pdf(file = "box_simu_extremes.pdf", width = 7, height = 3.5, onefile = F)
ggarrange(p1, p2, ncol = 2, common.legend = T, legend = "bottom")
dev.off()

aggregate(out_mat[,1:3], list(out_mat$scenario, out_mat$model), function(x) mean(x))# sum(x > 0)/100)

