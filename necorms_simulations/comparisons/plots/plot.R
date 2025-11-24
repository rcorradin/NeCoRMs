out_mat <- matrix(0, ncol = 5, nrow = 800)
out_temp_VI <- rep(0, 800)
out_temp_VI_nest <- rep(0, 800)
out_temp_sKL <- rep(0, 800)

load(file = "necorms_simulations/comparisons/low_groups_low_sep/simu_nested_low_group_low_sep.Rdata")
for(i in 1:100){
  out_temp_VI[i] <- out_list_100_LG_LS[[i]][1]
  out_temp_sKL[i] <- out_list_100_LG_LS[[i]][2]
  out_temp_VI_nest[i] <- out_list_100_LG_LS[[i]][3]
  out_mat[i,4] <- "LG-HS"
  out_mat[i,5] <- "nCoRMs"
  
  out_temp_VI[i+100] <- out_list_100_LG_LS[[i]][4]
  out_temp_sKL[i+100] <- out_list_100_LG_LS[[i]][5]
  out_temp_VI_nest[i+100] <- out_list_100_LG_LS[[i]][6]
  out_mat[i+100,4] <- "LG-HS"
  out_mat[i+100,5] <- "CAM"
  
}

load(file = "necorms_simulations/comparisons/low_groups_high_sep/simu_nested_low_group_high_sep.Rdata")
for(i in 1:100){
  out_temp_VI[i+200] <- out_list_100_LG_HS[[i]][1]
  out_temp_sKL[i+200] <- out_list_100_LG_HS[[i]][2] 
  out_temp_VI_nest[i+200] <- out_list_100_LG_HS[[i]][3]
  out_mat[i+200,4] <- "LG-HS"
  out_mat[i+200,5] <- "nCoRMs"
  
  out_temp_VI[i+300] <- out_list_100_LG_HS[[i]][4]
  out_temp_sKL[i+300] <- out_list_100_LG_HS[[i]][5] 
  out_temp_VI_nest[i+300] <- out_list_100_LG_HS[[i]][6]
  out_mat[i+300,4] <- "LG-HS"
  out_mat[i+300,5] <- "CAM"
}
out_mat <- as.data.frame(out_mat)
out_mat[,1] <- out_temp_VI
out_mat[,2] <- out_temp_sKL
out_mat[,3] <- out_temp_VI_nest

load(file = "necorms_simulations/comparisons/high_groups_low_sep/simu_nested_high_group_low_sep.Rdata")
for(i in 1:100){
  out_temp_VI[i+400] <- out_list_100_HG_LS[[i]][1]
  out_temp_sKL[i+400] <- out_list_100_HG_LS[[i]][2] 
  out_temp_VI_nest[i+400] <- out_list_100_HG_LS[[i]][3]
  out_mat[i+400,4] <- "HG-LS"
  out_mat[i+400,5] <- "nCoRMs"
  
  out_temp_VI[i+500] <- out_list_100_HG_LS[[i]][5]
  out_temp_sKL[i+500] <- out_list_100_HG_LS[[i]][6] 
  out_temp_VI_nest[i+500] <- out_list_100_HG_LS[[i]][7]
  out_mat[i+500,4] <- "HG-LS"
  out_mat[i+500,5] <- "CAM"
}
out_mat <- as.data.frame(out_mat)
out_mat[,1] <- out_temp_VI
out_mat[,2] <- out_temp_sKL
out_mat[,3] <- out_temp_VI_nest

load(file = "necorms_simulations/comparisons/high_groups_high_sep/simu_nested_high_group_high_sep.Rdata")
for(i in 1:100){
  out_temp_VI[i+600] <- out_list_100_HG_HS[[i]][1]
  out_temp_sKL[i+600] <- out_list_100_HG_HS[[i]][2] 
  out_temp_VI_nest[i+600] <- out_list_100_HG_HS[[i]][3]
  out_mat[i+600,4] <- "HG-HS"
  out_mat[i+600,5] <- "nCoRMs"
  
  out_temp_VI[i+700] <- out_list_100_HG_HS[[i]][4]
  out_temp_sKL[i+700] <- out_list_100_HG_HS[[i]][5] 
  out_temp_VI_nest[i+700] <- out_list_100_HG_HS[[i]][6]
  out_mat[i+700,4] <- "HG-HS"
  out_mat[i+700,5] <- "CAM"
}
out_mat <- as.data.frame(out_mat)
out_mat[,1] <- out_temp_VI
out_mat[,2] <- out_temp_sKL
out_mat[,3] <- out_temp_VI_nest

colnames(out_mat) <- c("VI", "sKL", "VInest", "scenario", "model")
out_mat$model <- factor(out_mat$model, levels = c("nCoRMs", "CAM"))
out_mat$scenario <- factor(out_mat$scenario, levels = c("LG-HS", "LG-LS", "HG-HS", "HG-LS"))

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

pdf(file = "box_simu_nested.pdf", width = 8, height = 4, onefile = F)
ggarrange(p1, p2, ncol = 2, common.legend = T, legend = "bottom")
dev.off()

aggregate(out_mat[,1:3], list(out_mat$scenario, out_mat$model), function(x) mean(x))# sum(x > 0)/100)

