res_mat <- matrix(0, ncol = 2, nrow = 400)

load("simu_compound_50_high.Rdata")
for(i in 1:100){
  res_mat[i, 1:2] <- out_list_50_high[[i]]
}

load("simu_compound_200_high.Rdata")
for(i in 1:100){
  res_mat[i+200, 1:2] <- out_list_200_high[[i]]
}

load("simu_compound_50_low.Rdata")
for(i in 1:100){
  res_mat[i+100, 1:2] <- out_list_50_low[[i]]
}

load("simu_compound_200_low.Rdata")
for(i in 1:100){
  res_mat[i+300, 1:2] <- out_list_200_low[[i]]
}

#---

library(ggplot2)
library(ggpubr)
blue <- "#100f7a"
red <- "#7a0f10"
dgray <- "#999999"

df_pl <- data.frame(VI = res_mat[,1], sKL = res_mat[,2], size = factor(rep(c("50", "50", "200", "200"), each = 100), 
                                                                       levels = c("50", "200")), 
                    separation = factor(rep(c("HIGH", "LOW", "HIGH", "LOW"), each = 100), levels = c("LOW", "HIGH")))
p1 <- ggplot(df_pl) + 
  geom_boxplot(mapping = aes(y = VI, x = size, fill = separation), alpha = 0.5) +
  theme_classic() +
  ylab("variation of information") + 
  scale_fill_manual(values = c(dgray, red))
p2 <- ggplot(df_pl) + 
  geom_boxplot(mapping = aes(y = sKL, x = size, fill = separation), alpha = 0.5) +
  theme_classic() + 
  ylab("J-divergence") + 
  scale_fill_manual(values = c(dgray, red))  
  # geom_hline(aes(yintercept = 1))
  # scale_fill_grey()

pdf(file = "box_simu_compound.pdf", width = 6.5, height = 3, onefile = F)
ggarrange(p1, p2, ncol = 2, common.legend = T, legend = "bottom")
dev.off()
