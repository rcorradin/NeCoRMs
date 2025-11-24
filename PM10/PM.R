library(ggplot2)
library(latex2exp)

Rcpp::sourceCpp('main_nested_CPP_with_quantile.cpp')
Rcpp::sourceCpp('utilities.cpp')

# --------
# IMPORT THE DATA

filenames <- list.files(path = "PM10/data", pattern="*.csv", full.names=TRUE)
data_list <- lapply(lapply(filenames, function(x) read.csv(x, skip = 4)[,2]), function(y) y[y > 0])
data <- as.numeric(unlist(data_list))
group <- unlist(apply(cbind(1:12, unlist(lapply(data_list, length))), 1, function(x) rep(x[1], x[2]))) - 1
grid <- seq(0.1, 120, length.out = 100)
grid_std <- seq(-4, 6, length.out = 100)
sdata = scale(data)
ldata = log(data)
gridl = log(grid)
# --------

set.seed(42)
mod <- main_nested(Y = data, group = group, d = 12, niter = 7500, nburn = 2500, thin = 1, beta = 1,
                   m0 = 0, k0 = 0.1, a0 = 2, b0 = 0.2, m1 = 0, s21 = 10, tau1 = 2, tau2 = 2,
                   a1 = 2, b1 = 2, sigma = 0.2, phi = 1, a_sigma = 2, b_sigma = 2, 
                   a_phi = 1, b_phi = 1, grid = grid, MH_var2 = 0.15, 
                   M = 100, tol = 0.0001, max_iter = 100, epsilon = 0.01, eval_density = T, 
                   IS_nval = 1000, nupd = 100, quantile = 50)

opt_part_nested  <- as.numeric(as.factor(mod$nested[which.min(VI_LB(clean_partition(mod$nested), psm_mat = psm(mod$nested))),]+1))
table(opt_part_nested)

# --------

grey <- '#3c3c3c'
lgray <- '#EFEFEF'
blue <- "#100f7a"
red <- "#7a0f10"

df_plot <- data.frame(
  x = rep(grid, 12),
  y = as.vector(sapply(1:12, function(x) apply(mod$dens[x,,], 1, mean))),
  ylow = as.vector(sapply(1:12, function(x) apply(mod$dens[x,,], 1, quantile, p = 0.05))),
  yup = as.vector(sapply(1:12, function(x) apply(mod$dens[x,,], 1, quantile, p = 0.95))),
  group_label = factor(rep(1:12, each = length(grid)))
)

hist_data <- data.frame(x = data, group_label = factor(group + 1))

qmat <- 1 - mod$quantile
temp_vals <- round(colMeans(qmat), digits = 2)
Aval <- sapply(temp_vals, function(x) expression(paste0("bar(Q)[i,50] =" , x)))
data_ann <- data.frame(city = c("Bergamo", "Brescia", "Como", "Cremona", "Lecco",
        "Lodi", "Milano", "Mantova", "Monza", "Pavia", "Sondrio", "Varese"),
        # labels_ann = sapply(round(apply(qmat, 2 , mean), digits = 2), function(y) paste0("bar(Q)[i,50] = ", y)),
        labels_ann = as.character(paste0(expression(bar(Q)[50]))),
        labels_ann1 = round(apply(qmat, 2 , mean), digits = 2),
        labels_ann2 = paste0("(", sprintf("%.2f",apply(qmat, 2 , quantile, p = c(0.025))), " - ",
                             sprintf("%.2f",apply(qmat, 2 , quantile, p = c(0.975))), ")"),
        group_label = factor(1:12))
# --------
sort_ord <- c(which(opt_part_nested == 1), which(opt_part_nested == 2), which(opt_part_nested == 3), which(opt_part_nested == 4),
              which(opt_part_nested == 5), which(opt_part_nested == 6), which(opt_part_nested == 7))
psm_nested <- psm(mod$nested[,sort_ord])

plot_all <- ggplot() +
  geom_line(data = df_plot, aes(x = x, y = y), colour = blue) +
  geom_histogram(data = hist_data, mapping = aes(x = x, stat(density)), color = grey, alpha=0.0, bins = 40) +
  geom_ribbon(data = df_plot, mapping = aes(x = x, ymin = ylow, ymax = yup), alpha = 0.3, fill = blue) +
  facet_wrap(~factor(group_label, levels = sort_ord), ncol = 4, nrow = 3, labeller = as_labeller(
    c("1" = "Bergamo",
      "2" = "Brescia",
      "3" = "Como",
      "4" = "Cremona",
      "5" = "Lecco",
      "6" = "Lodi",
      "7" = "Milano",
      "8" = "Mantova",
      "9" = "Monza",
      "10" = "Pavia",
      "11" = "Sondrio",
      "12" = "Varese")
  )) +
  geom_vline(xintercept = 50, color = red) +
  geom_text(data = data_ann, mapping = aes(x = 88, y = 0.04, label = labels_ann), size = 4, parse = T) +
  geom_text(data = data_ann, mapping = aes(x = 110, y = 0.04, label = paste0(" = ", labels_ann1)), size = 4) +
  geom_text(data = data_ann, mapping = aes(x = 103, y = 0.035, label = paste0("CI = ",labels_ann2)), size = 2.25) +
  theme_minimal() +
  xlab(expression(PM[1][0])) +
  ylab("") +
  xlim(c(0, 150)) +
  ylab("density")

pdf(file = "all_dens.pdf", width = 12, height = 7)
plot_all
dev.off()

# ----------------
matrix_nested <- reshape2::melt(psm_nested[1:nrow(psm_nested), , drop = FALSE])
plot_heat <- ggplot(matrix_nested) +
  geom_tile(aes(x = Var2, y = Var1, fill = value), na.rm = TRUE, size = 0.0,show.legend = FALSE) +
  scale_fill_gradient(
    low = "transparent",
    high = "black"
  ) +
  labs(x = "", y = "", fill = "") +
  theme_linedraw() +
  coord_cartesian(xlim = c(0, ncol(psm_nested) + 1), ylim = c(0, ncol(psm_nested) + 1), expand = FALSE) +
  scale_x_continuous(breaks = seq(1, 12, by = 1),
                     labels = c("Bergamo", "Brescia", "Como", "Cremona", "Lecco",
                                "Lodi", "Milano", "Mantova", "Monza", "Pavia", "Sondrio", "Varese")[sort_ord]) +
  scale_y_continuous(breaks = seq(1, 12, by = 1),
                     labels = c("Bergamo", "Brescia", "Como", "Cremona", "Lecco",
                                "Lodi", "Milano", "Mantova", "Monza", "Pavia", "Sondrio", "Varese")[sort_ord]) +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(size=0.1))

pdf(file = "heat_PM.pdf", width = 3.25, height = 3)
plot_heat
dev.off()

# ----------------
plt_data_diag <- data.frame(x = (1:length(mod$n_eta)) + 10000,
                            y = log(mod$phi),
                            z = log(mod$sigma))
pd2 <- ggplot(plt_data_diag) +
  geom_line(aes(x = x, y = z)) +
  labs(x = "iter", y = expression(sigma)) +
  theme_minimal() +
  geom_hline(mapping = aes(yintercept = mean(z)), col = red, lty = 2) +
  ggtitle(label = expression("Traceplot of " * sigma))

pd3 <- ggplot(plt_data_diag) +
  geom_line(aes(x = x, y = y)) +
  labs(x = "iter", y = expression(phi)) +
  theme_minimal() +
  geom_hline(mapping = aes(yintercept = mean(y)), col = red, lty = 2) +
  ggtitle(label = expression("Traceplot of " * phi))

pdf(file = "diag_PM.pdf", width = 4, height = 3)
pd2
dev.off()

pdf(file = "diag_PM2.pdf", width = 4, height = 3)
pd3
dev.off()
