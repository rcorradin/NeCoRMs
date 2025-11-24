load(file="simu_nested_diff_params.RData")

time_mat <- dist_mat <- dens_mat <- part_mat <- matrix(0, ncol = 12, nrow = 100)
for(i in 1:100){
  for(j in 1:12){
    time_mat[i,j] <- out_list[[i]][[j]][4]
    dist_mat[i,j] <- out_list[[i]][[j]][3]
    dens_mat[i,j] <- out_list[[i]][[j]][2]
    part_mat[i,j] <- out_list[[i]][[j]][1]
  }
}

library(ggplot2)
library(ggpubr)
blue <- "#0a4f16"
mblue <- "#519c5f"
lblue <- "#97dba4"
red <- "#7a0f10"
dgray <- "#999999"

data_plot <- data.frame(x = as.vector(part_mat), 
                        y = as.factor(rep(rep(c(0.5, 1, 5), each = 100), times = 4)), 
                        z = as.factor(rep(c(0.01, 0.2, 0.5, 0.75), each = 300)))

part_plot <- ggplot(data_plot) + 
  geom_boxplot(mapping = aes(y = x, x = y, fill = z), alpha = 0.75) +
  theme_classic() + 
  xlab(expression(phi)) +
  ylab("variation of information") + 
  theme(legend.position = "null") +
  scale_fill_manual(values = c("white", lblue, mblue, blue), name = expression(sigma))

data_plot <- data.frame(x = as.vector(dens_mat), 
                        y = as.factor(rep(rep(c(0.5, 1, 5), each = 100), times = 4)), 
                        z = as.factor(rep(c(0.01, 0.2, 0.5, 0.75), each = 300)))
idx <- which(data_plot[,2] == "0.5" & data_plot[,3] == "0.75" & data_plot[,1] >= 1)
data_plot <- data_plot[-idx,]

dens_plot <- ggplot(data_plot) + 
  geom_boxplot(mapping = aes(y = x, x = y, fill = z), alpha = 0.75) +
  theme_classic() +
  xlab(expression(phi)) +
  ylab("J-divergence") +
  theme(legend.position = "null") +
  scale_fill_manual(values = c("white", lblue, mblue, blue), name = expression(sigma))

pdf(file = "box_simu_nested.pdf", width = 8, height = 3.5, onefile = F)
ggarrange(part_plot, dens_plot, ncol = 2, common.legend = T, legend = "bottom")
dev.off()
