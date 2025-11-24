#--------------------------------------------------------------------------
# ANALISI dati galaxy_new -------------------------------------------------
#--------------------------------------------------------------------------

Rcpp::sourceCpp('utilities.cpp')
Rcpp::sourceCpp('main_compound_CPP.cpp')
library(ggplot2)

data <- read.csv("galaxy_data.csv")
data <- data[(data[,1] > 0 & data[,1] < 4),]
clas_lum <- as.numeric(as.factor(as.numeric(cut(10^data[,3], c(0, 0.2, 0.5, 1.4, 6.5, 900))))) - 1
clas_den <- as.numeric(as.factor(as.numeric(cut(data[,2], c(-26, -22, -21, -20, -19, 0))))) - 1
group_tot <- as.numeric(as.factor(apply(cbind(clas_lum, clas_den), 1, function(x) paste(x[1], "x", x[2])))) - 1

data_use  <- data[data[,1] > 0 & data[,1] < 4 & clas_lum == 3,1]
group_use <- clas_den[data[,1] > 0 & data[,1] < 4 & clas_lum == 3]
data_tot  <- data[data[,1] > 0 & data[,1] < 4,1]

grid <- seq(-2, 6, length.out = 100)

# label per plot ----------------------------------------------------------

cl_lum <- as.factor(as.numeric(cut(10^data[,3], c(0, 0.2, 0.5, 1.4, 6.5, 900))))
cl_den <- as.factor(as.numeric(cut(data[,2], c(-26, -22, -21, -20, -19, 0))))

levels(cl_den) <- c("Luminosity = H", "Luminosity = H/M", "Luminosity = M", "Luminosity = M/L", "Luminosity = L")
levels(cl_lum) <- c("Environment = L", "Environment = L/M", "Environment = M", "Environment = M/H", "Environment = H")

temp <- data.frame(group_tot, cl_lum, cl_den)
label_match <- unique(temp)

# DDP model ---------------------------------------------------------------
sum(mod$phi[-1,1] == mod$phi[-nrow(mod$phi),1]) / nrow(mod$phi)

set.seed(42)
system.time(mod <- main_compound(Y = data_tot, group = group_tot, d = 25, niter = 15000, nburn = 10000, thin = 1,
                                 m0 = 0, k0 = 0.1, a0 = 2, b0 = 0.2, m1 = 0, s21 = 10, tau1 = 2, tau2 = 0.2,
                                 a1 = 2, b1 = 0.2, sigma = 0.2, phi = 0.5, a_sigma = 2, b_sigma = 2, 
                                 a_phi = 2, b_phi = 2, grid = grid, 
                                 M = 1000, tol = 0.0001, max_iter = 100, epsilon = 0.00075, eval_density = T, 
                                 IS_nval = 1000, nupd = 50, MH_var = 0.1, MH_var2 = 0.25))
# save.image(file = "mod_galaxy_est.RData")
load(file = "mod_galaxy_est2.RData")
opt_part  <- mod$clust[which.min(VI_LB(clean_partition(mod$clust), psm_mat = psm(mod$clust))),]+1
table(opt_part)[table(opt_part) > 1000]
opt_part[!(opt_part %in% c(1,2,31,33,35))] <- 200
opt_part <- as.numeric(as.factor(opt_part))

plot_data <- data.frame(
  x = rep(grid, 25), 
  y = as.vector(sapply(1:25, function(x) apply(mod$dens[x,,], 1, mean))), 
  ylow = as.vector(sapply(1:25, function(x) apply(mod$dens[x,,], 1, quantile, p = 0.025))), 
  yup = as.vector(sapply(1:25, function(x) apply(mod$dens[x,,], 1, quantile, p = 0.975))),
  group = as.factor(rep(1:25, each = length(grid)))
)

ann_data <- matrix(NA, ncol = 3, nrow = 25)

names(plot_data) <- c("grid", "V2", "Q1", "Q2", "V3")
# names(mar_data) <- c("grid", "V2", "Q1", "Q2", "V3")
plot_data$V3 <- as.factor(plot_data$V3)
data_hist <- as.data.frame(cbind(data_tot, as.factor(as.numeric(group_tot))))
names(data_hist) <-  c("V1", "V3")
data_hist$V3 <- as.factor(group_tot)

ann_df <- as.data.frame(cbind(ann_data, 0:24))
names(ann_df) <- c("V1", "V3")

labs1 <- label_match[unlist(sapply(as.numeric(plot_data$V3), function(x) which(label_match[,1] + 1 == x))),c(2,3)]
labs2 <- label_match[unlist(sapply(as.numeric(mar_data$V3), function(x) which(label_match[,1] == x))),c(2,3)]
labs3 <- label_match[unlist(sapply(as.numeric(ann_df$V3), function(x) which(label_match[,1] + 1 == x))),c(2,3)]

plot_data <- cbind(plot_data, labs1)
# mar_data <- cbind(mar_data, labs2)
ann_df <- cbind(ann_df, labs3)
data_hist <- data.frame(V1 = data_hist[,1], cl_lum, cl_den)

grey <- '#3c3c3c'
lgray <- '#EFEFEF'
blue <- "#100f7a"
red <- "#7a0f10"

drsimonj_colors <- c("#d11141", "#00b159", "#00aedb", "#f37735",
    "#ffc425", "#cccccc", "#8c8c8c")

pl <- ggplot() +
  geom_histogram(data = data_hist, mapping = aes(x = V1, stat(density)), color = grey, alpha=0.0, bins = 40) +
  geom_line(data = plot_data, aes(x = grid, y = V2),  lwd = 0.4,show.legend=F, color = blue) +
  # geom_line(data = mar_data, aes(x = grid, y = V2, color = 1),  lwd = 0.4,show.legend=F, color = 1, lty = 2) +
  # geom_text(data = ann_df, mapping = aes(x = 1.5, y = 3, label = V1), size = 2.2) +
  theme_minimal() +
  # theme_classic() +
  # geom_line(data = plot_data, aes(x = grid,y = Q1), lty = 1, col = "#100f7a") +
  # geom_line(data = plot_data, aes(x = grid,y = Q2), lty = 1, col = "#100f7a") +
  geom_ribbon(data = plot_data, aes(x = grid,ymin = Q1, ymax = Q2), alpha = 0.5, fill = blue) +
  # geom_ribbon(data = mar_data, aes(x = grid,ymin = Q1, ymax = Q2), alpha = 0.3, fill = 1) +
  scale_x_continuous(limits = c(-0.0, 3.5)) +
  theme(strip.background = element_blank()) + #, strip.text.y = element_blank(), strip.text.x = element_blank()) +
  # theme(strip.background = element_rect(fill=lgray)) +
  facet_grid(cl_den~cl_lum) +
  # labs(x = "Luminosity - Low to High", y = "Density - Low to High") 
  labs(x = "y", y = "density")  + 
  theme(legend.position = "null") + 
  scale_colour_manual(values = drsimonj_colors)

pdf(file = "dens_compound.pdf", width = 9, height = 6)
pl
dev.off()

pdf(file = "dens_compound2.pdf", width = 9, height = 6)
pl +
  geom_segment(data = data.frame(x = data_tot, cval = as.factor(opt_part), cl_den = cl_den, cl_lum = cl_lum), 
               mapping = aes(x = x, y= -0.05, xend = x, yend = 0.05, color = cval))
dev.off()

# plot_ancestor

plt_data_ancestor <- data.frame(x = grid, 
                                y = colMeans(mod$ancestor),
                                ylow = apply(mod$ancestor, 2, quantile, p = 0.025),
                                yup = apply(mod$ancestor, 2, quantile, p = 0.975))
pla <- ggplot() + 
  geom_histogram(data = data_hist, mapping = aes(x = V1, stat(density)), color = grey, alpha=0.0, bins = 50) +
  geom_line(data = plt_data_ancestor, aes(x = x, y = y), color = red) +
  # geom_line(data = plt_data_ancestor, aes(x = x,y = ylow), lty = 1, col = "#100f7a") +
  # geom_line(data = plt_data_ancestor, aes(x = x,y = yup), lty = 1, col = "#100f7a") +
  geom_ribbon(data = plt_data_ancestor, aes(x = x, ymin = ylow, ymax = yup), alpha = 0.5, fill = red) +
  scale_x_continuous(limits = c(-0.0, 3.5)) +
  labs(x = "y", y = "density") +
  theme_minimal()

pdf(file = "ancestor_compound.pdf", width = 4, height = 3)
pla
dev.off()

plt_data_diag <- data.frame(x = 1:length(mod$n_eta) + 10000, 
                                y = mod$phi,
                                z = mod$sigma)
pd1 <- ggplot(plt_data_diag) + 
  geom_line(aes(x = x, y = y)) +
  labs(x = "iter", y = expression(phi)) +
  theme_minimal() +
  geom_hline(mapping = aes(yintercept = mean(y)), col = red, lty = 2) + 
  ggtitle(label = expression("Traceplot of " * phi))

pd2 <- ggplot(plt_data_diag) + 
  geom_line(aes(x = x, y = z)) +
  labs(x = "iter", y = expression(sigma)) +
  theme_minimal() +
  geom_hline(mapping = aes(yintercept = mean(z)), col = red, lty = 2) + 
  ggtitle(label = expression("Traceplot of " * sigma))

pdf(file = "diag_galaxy1.pdf", width = 4, height = 3)
pd1
dev.off()

pdf(file = "diag_galaxy2.pdf", width = 4, height = 3)
pd2
dev.off()

