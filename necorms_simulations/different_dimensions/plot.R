load(file="simu_nested_dims.RData")

time_mat <- dist_mat <- dens_mat <- part_mat <- matrix(0, ncol = 3, nrow = 100)
for(j in 1:3){
  for(i in 1:100){
    time_mat[i,j] <- out_list[[i]][[j]][4]
    dist_mat[i,j] <- out_list[[i]][[j]][3]
    dens_mat[i,j] <- out_list[[i]][[j]][2]
    part_mat[i,j] <- out_list[[i]][[j]][1]
  }
}

library(ggplot2)
library(ggpubr)
blue <- "#992826"
mblue <- "#d68281"
red <- "#7a0f10"
dgray <- "#999999"

data_plot <- data.frame(x = as.vector(part_mat), 
                        y = factor(rep(rep(c("p = 2", "p = 5", "p = 10"), each = 100), times = 4), levels = c("p = 2", "p = 5", "p = 10")))

part_plot <- ggplot(data_plot) + 
  geom_boxplot(mapping = aes(y = x, x = y, fill = y), alpha = 0.75) +
  theme_classic() + 
  xlab("") +
  ylab("variation of information") + 
  theme(legend.position = "null") +
  scale_fill_manual(values = c("white", mblue, blue))

data_plot <- data.frame(x = as.vector(dens_mat), 
                        y = factor(rep(rep(c("d = 2", "d = 5", "d = 10"), each = 100), times = 4), levels = c("d = 2", "d = 5", "d = 10")))

dens_plot <- ggplot(data_plot) + 
  geom_boxplot(mapping = aes(y = x, x = y, fill = y), alpha = 0.75) +
  theme_classic() + 
  xlab("") +
  ylab("J-diver") + 
  theme(legend.position = "null") +
  scale_fill_manual(values = c("white", mblue, blue))

pdf(file = "box_simu_nested.pdf", width = 6, height = 3, onefile = F)
ggarrange(part_plot, dens_plot, ncol = 2)
dev.off()
