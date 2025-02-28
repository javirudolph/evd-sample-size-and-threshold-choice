
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(purrr)
library(MetBrewer)
library(stringr)
library(latex2exp)


source("functions.R")


# FIGURE 1 ---------

# load("data/allinone_original/mixturesamples.RData")

desired_means <- c(28, 32, 40, 50)
desired_sds <- c(49.7, 39.9, 33.3, 31.1)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)
# These are the weights for the distribution
pis <- c( 0.1, 0.2, 0.3, 0.4)

# The color associated to each mixture component
dens_cols <- c("#264653", "#2a9d8f", "#f4a261", "#e76f51")
# dens_cols <- pal_mixture

# Build the densities for plotting
weighted_densities <- purrr::map(1:4, function(y) stat_function(fun = dlnorm, args = list(meanlog = pars$meanlog[y], sdlog = pars$sdlog[y]), color = dens_cols[y], size=1, alpha = 0.8))

# Plot density components of the mixture
ggplot() +
   weighted_densities +
   theme_bw(base_size = 7) +
   labs(y = "Component Densities") +
   lims(x = c(0, 150)) -> densities_plot

# Plot the means and sd for each of the curves
data.frame(desired_means, desired_sds, gID = factor(1:4)) %>%
   mutate(lo = desired_means - desired_sds,
          hi = desired_means + desired_sds) %>%
   ggplot(., aes(x = gID, y = desired_means, color = gID)) +
   geom_point(size = 2) +
   geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.1, size = 1) +
   scale_color_manual(values = dens_cols) +
   labs(y = "Mean +/- SD") +
   theme_bw(base_size = 7) +
   theme(axis.title.x = element_blank(),
         axis.text.x = element_blank(),
         legend.position = "none") -> means_sd_plot


tru_tail <- data.frame(values = truth_df$x_samps, y = 100) %>% arrange(desc(values)) %>% filter(values >=50)

# Histogram plus the points in the tail.
thresh_tests <- c(50, 100, 150, 250, 500, 750, 1000)
truth_df %>%
   ggplot(., aes(x = x_samps)) +
   geom_histogram(bins = 100) +
   geom_point(data = tru_tail[1:50,], aes(x = values, y = y), color = "grey47", alpha = 0.9, size = 1) +
   geom_vline(xintercept = thresh_tests, color = dens_cols[4], linetype = "dashed", size = 0.3) +
   labs(y = "Frequency", x = "Distance") +
   #lims(x = c(0, 150)) +
   theme_bw(base_size = 7) -> truth_hist
# truth_hist

# Density curve for the mixture based on the 50k samples.
truth_df %>%
   ggplot(., aes(x = x_samps)) +
   geom_density() +
   labs(y = "Mixture Density", x = "Distance") +
   lims(x = c(0, 150), y = c(0,0.05)) +
   theme_bw(base_size = 7) +
   theme(axis.title.x = element_blank()) -> truth_density
# truth_density

# Visualize both in one frame
# plot_grid(truth_hist, truth_density)

## Assemble -----------------------------


top_row <- plot_grid(means_sd_plot, densities_plot, truth_density, nrow = 1, rel_widths = c(1,2,2), labels = "AUTO", label_size = 8)

plot_grid(top_row, truth_hist, nrow = 2, labels = c("", "D"), label_size = 8)

ggsave(filename = "Figures/Figure1.png", width = 6, height = 3, scale = 1.5, bg = "white")



# FIGURE 2 ------------------------
# See the estimations compared to original data before the correction

# load("data/allinone_original/nreps_mles_df.RData")
# load("data/allinone_original/GP_runs.RData")
pal_sampsize <- met.brewer("Hokusai3", n = 6)

nreps_mles_df %>%     
   mutate(samp_n_tests = factor(samp_n_tests),
          thresh_tests = factor(thresh_tests)) %>% 
   ggplot(., aes(x = thresh_tests, y = log_St_hat - log(theta), fill = samp_n_tests, color = samp_n_tests)) +
   # facet_wrap(~samp_n_tests, nrow = 1) +
   geom_boxplot(alpha = 0.9) +
   geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
   scale_fill_manual(values = pal_sampsize) +
   scale_color_manual(values = pal_sampsize) +
   labs(x = " ", y = TeX(r'($log(\hat{\theta} / \theta)$)'), subtitle = "Simple Lomax") +
   theme_bw(base_size = 8) +
   coord_cartesian(ylim = c(-23, 10)) +
   theme(legend.position = "none") -> a

nreps_mles_df %>% 
   mutate(samp_n_tests = factor(samp_n_tests),
          thresh_tests = factor(thresh_tests)) %>% 
   ggplot(., aes(x = thresh_tests, y = log_St_hat2 - log(theta), fill = samp_n_tests, color = samp_n_tests)) +
   # facet_wrap(~samp_n_tests, nrow = 1) +
   geom_boxplot(alpha = 0.9) +
   geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
   scale_fill_manual(values = pal_sampsize, name = "Sample Size") +
   scale_color_manual(values = pal_sampsize, name = "Sample Size") +
   labs(x = "", y = TeX(r'($log(\hat{\theta} / \theta)$)'), subtitle = "GLM Lomax") +
   theme_bw(base_size = 8) +
   coord_cartesian(ylim = c(-23, 10)) +
   theme(legend.position = "right") -> b

ableg <- get_legend(b)

nreps_mles_df %>% 
   mutate(samp_n_tests = factor(samp_n_tests),
          thresh_tests = factor(thresh_tests)) %>% 
   ggplot(., aes(x = thresh_tests, y = log(GP_theta_hat) - log(theta), fill = samp_n_tests, color = samp_n_tests)) +
   # facet_wrap(~samp_n_tests, nrow = 1) +
   geom_boxplot(alpha = 0.9) +
   geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
   scale_fill_manual(values = pal_sampsize, name = "Sample Size") +
   scale_color_manual(values = pal_sampsize, name = "Sample Size") +
   labs(x = "", y = TeX(r'($log(\hat{\theta} / \theta)$)'), subtitle = "Generalized Pareto") +
   theme_bw(base_size = 8) +
   coord_cartesian(ylim = c(-23, 10)) +
   theme(legend.position = "none") -> c

# GP_mles %>% 
#    mutate(samp_n_tests = factor(samp_n_tests),
#           thresh_tests = factor(thresh_tests)) %>% 
#    ggplot(., aes(x = thresh_tests, y = log(quant_theta/theta), fill = samp_n_tests, color = samp_n_tests)) +
#    # facet_wrap(~samp_n_tests, nrow = 1) +
#    geom_boxplot(alpha = 0.9) +
#    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
#    scale_fill_manual(values = pal_sampsize) +
#    scale_color_manual(values = pal_sampsize) +
#    labs(x = "Distance level", y = TeX(r'($log(\hat{\theta} / \theta)$)'), subtitle = "Generalized Pareto") +
#    theme_bw(base_size = 8) +
#    theme(legend.position = "none") -> c

abplot <- plot_grid(a, b + theme(legend.position = "none"), c, nrow = 3, labels = "AUTO", label_size = 8)

plot_grid(abplot, ableg, ncol = 2, rel_widths = c(1, 0.2))
ggsave("Figures/Figure2.png", width = 6, height = 3, scale = 1.5, bg = "white")


# FIGURE 3 -----------------------
# Introduce the numerical estimation with the simple lomax

nreps_mles_df %>% 
   ggplot()+
   geom_point(aes(x = alpha_star, y = k_star, color = "Simple Lomax"), size = 4, alpha = 0.5) +
   geom_point(aes(x = alpha_glm, y = k_glm, color = "Lomax GLM"), size = 4, alpha = 0.5) +
   scale_y_log10() +
   scale_x_log10() +
   scale_color_manual(values = c("grey20", "dodgerblue4"), name = "") +
   labs(x = "Log(alpha)", y = "Log(k)") +
   theme_bw() +
   theme(legend.position = "bottom") -> a

nreps_mles_df %>% 
   mutate(labelsamp = paste("n =", samp_n_tests),
          labelsamp = factor(labelsamp, levels = c("n = 80", "n = 200", "n = 500", "n = 800", "n = 1000", "n = 1600"))) %>% 
   ggplot()+
   facet_wrap(~labelsamp) +
   geom_point(aes(x = alpha_star, y = k_star, color = "Simple Lomax"), size = 4, alpha = 0.5) +
   geom_point(aes(x = alpha_glm, y = k_glm, color = "Lomax GLM"), size = 4, alpha = 0.5) +
   scale_y_log10() +
   scale_x_log10() +
   scale_color_manual(values = c("grey20", "dodgerblue4"), name = "") +
   labs(x = "Log(alpha)", y = "Log(k)") +
   theme_bw() +
   theme(legend.position = "none",
         axis.text = element_blank()) -> b

ableg <- get_legend(a)

top <- plot_grid(a + theme(legend.position = "none"), b, ncol = 2, labels = "AUTO", label_size = 8)

plot_grid(top, ableg, nrow = 2, rel_heights = c(1, 0.1))
ggsave("Figures/Figure3.png", width = 4, height = 3, scale = 1.5, bg = "white")


# FIGURE 4 -------------------------

# load("data/allinone_original/allin1df.RData")
allin1df %>%
   dplyr::select(., thresh_tests, samp_n_tests, gp_theta_hat, gp_theta_bar, theta) %>%
   mutate(original = log(gp_theta_hat/theta),
          corrected = log(gp_theta_bar/theta)) %>%
   pivot_longer(., cols = c(original, corrected), names_to = "type", values_to = "theta_val") %>%
   drop_na() %>%
   mutate(thresh_tests = factor(thresh_tests),
          samp_n_tests = factor(samp_n_tests)) %>% 
   mutate(labelsamp = paste("n =", samp_n_tests),
          labelsamp = factor(labelsamp, levels = c("n = 80", "n = 200", "n = 500", "n = 800", "n = 1000", "n = 1600"))) %>%
   ggplot(., aes(x = thresh_tests, y = theta_val, color = type, fill = type)) +
   facet_wrap(~labelsamp, nrow = 3) +
   geom_boxplot() +
   scale_color_manual(values = dens_cols[c(1,4)], name = "Bias") +
   scale_fill_manual(values = dens_cols[c(2,3)], name = "Bias") +
   geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
   labs(y = TeX(r'($log(\hat{\theta} / \theta)$)'), x = "Distance Level Estimated") +
   theme_bw()
ggsave("Figures/Figure4.png", width = 4, height = 4, scale = 1.5, bg = "white")



# FIGURE 5 --------------------------

mixtures_plot_fx <- function(pars = pars, dens_cols = dens_cols, pis = pis, desired_means = desired_means, desired_sds = desired_sds, weighted = FALSE){
   
   if(weighted == TRUE){
      lnorm_densities <- purrr::map(1:4, function(y) stat_function(fun = dlnorm, args = list(meanlog = pars$meanlog[y], sdlog = pars$sdlog[y]), color = dens_cols[y], size=pis[y]*5, alpha = 0.8))
   } else{
      lnorm_densities <- purrr::map(1:4, function(y) stat_function(fun = dlnorm, args = list(meanlog = pars$meanlog[y], sdlog = pars$sdlog[y]), color = dens_cols[y], size=1, alpha = 0.8))
   }
   
   # Plot density components of the mixture
   ggplot() +
      lnorm_densities +
      theme_minimal(base_size = 6) +
      labs(y = "Density") +
      lims(x = c(0, 150)) -> densities_plot
   # densities_plot
   
   # Plot the means and sd for each of the curves
   data.frame(desired_means, desired_sds, gID = factor(1:4)) %>%
      mutate(lo = desired_means - desired_sds,
             hi = desired_means + desired_sds) %>%
      ggplot(., aes(x = gID, y = desired_means, color = gID)) +
      geom_point(size = 2) +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.1, size = 1) +
      scale_color_manual(values = dens_cols) +
      labs(y = "Mean +/- SD") +
      theme_minimal(base_size = 6) +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            legend.position = "none") -> means_sd_plot
   # means_sd_plot
   
   # Align these plots to have multiple panels lates
   plot_grid(NULL, means_sd_plot, densities_plot, rel_widths = c(0.1, 1, 3), nrow = 1)
}

bias_bxplot_fx <- function(data){
   data %>% 
      dplyr::select(., thresh_tests, samp_n_tests, gp_theta_hat, gp_theta_bar, theta) %>%
      mutate(original = log(gp_theta_hat/theta),
             corrected = log(gp_theta_bar/theta)) %>%
      pivot_longer(., cols = c(original, corrected), names_to = "type", values_to = "theta_val") %>%
      drop_na() %>%
      mutate(thresh_tests = factor(thresh_tests),
             samp_n_tests = factor(samp_n_tests)) %>% 
      mutate(labelsamp = paste("n =", samp_n_tests),
             labelsamp = factor(labelsamp, levels = c("n = 80", "n = 200", "n = 500", "n = 800", "n = 1000", "n = 1600"))) %>%
      ggplot(., aes(x = thresh_tests, y = theta_val, color = type, fill = type)) +
      # facet_wrap(~labelsamp, nrow = 1) +
      geom_boxplot() +
      scale_color_manual(values = dens_cols[c(1,4)], name = "Bias") +
      scale_fill_manual(values = dens_cols[c(2,3)], name = "Bias") +
      geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
      labs(y = TeX(r'($log(\hat{\theta} / \theta)$)'), x = " ") +
      theme_bw(base_size = 8) +
      theme(legend.position = "none")
}

# Scenarios:
pis <- c( 0.1, 0.2, 0.3, 0.4)
dens_cols <- pal_mixture

# SAME SD INCREASING MEAN ------------------------------
dir_scenario <- "ain1_samesd_upmean"

desired_means <- c(30, 40, 50, 60)
desired_sds <- c(25, 25, 25, 25)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)
A1 <- mixtures_plot_fx(pars = pars, dens_cols = dens_cols, pis = pis, desired_means = desired_means, desired_sds = desired_sds)

load(paste0("data/", dir_scenario, "/allin1df.RData"))
B1 <- bias_bxplot_fx(allin1df)

# SAME MEAN INCREASING SD ------------------------------

dir_scenario <- "ain1_samemean_upsd"

desired_means <- c(30, 30, 30, 30)
desired_sds <- c(20, 30, 40, 50)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)

A2 <- mixtures_plot_fx(pars = pars, dens_cols = dens_cols, pis = pis, desired_means = desired_means, desired_sds = desired_sds)

load(paste0("data/", dir_scenario, "/allin1df.RData"))
B2 <- bias_bxplot_fx(allin1df)


# INCREASING SD INCREASING MEAN ------------------------------
dir_scenario <- "ain1_upsd_upmean"

desired_means <- c(30, 40, 50, 60)
desired_sds <- c(20, 30, 40, 50)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)
A3 <- mixtures_plot_fx(pars = pars, dens_cols = dens_cols, pis = pis, desired_means = desired_means, desired_sds = desired_sds)

load(paste0("data/", dir_scenario, "/allin1df.RData"))
B3 <- bias_bxplot_fx(allin1df)


# DECREASING BOTH ------------------------------
dir_scenario <- "ain1_downsd_downmean"

desired_means <- c(60, 50, 40, 30)
desired_sds <- c(50, 40, 30, 20)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)
A4 <- mixtures_plot_fx(pars = pars, dens_cols = dens_cols, pis = pis, desired_means = desired_means, desired_sds = desired_sds)

load(paste0("data/", dir_scenario, "/allin1df.RData"))
B4 <- allin1df %>%
   dplyr::select(., thresh_tests, samp_n_tests, gp_theta_hat, gp_theta_bar, theta) %>%
   mutate(original = log(gp_theta_hat/theta),
          corrected = log(gp_theta_bar/theta)) %>%
   pivot_longer(., cols = c(original, corrected), names_to = "type", values_to = "theta_val") %>%
   drop_na() %>%
   mutate(thresh_tests = factor(thresh_tests),
          samp_n_tests = factor(samp_n_tests)) %>% 
   mutate(labelsamp = paste("n =", samp_n_tests),
          labelsamp = factor(labelsamp, levels = c("n = 80", "n = 200", "n = 500", "n = 800", "n = 1000", "n = 1600"))) %>%
   ggplot(., aes(x = thresh_tests, y = theta_val, color = type, fill = type)) +
   # facet_wrap(~labelsamp, nrow = 1) +
   geom_boxplot() +
   scale_color_manual(values = dens_cols[c(1,4)], name = "Bias") +
   scale_fill_manual(values = dens_cols[c(2,3)], name = "Bias") +
   geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
   labs(y = TeX(r'($log(\hat{\theta} / \theta)$)'), x = "Distance level estimated") +
   theme_bw(base_size = 8) +
   theme(legend.position = "none")

#INCREASING MEANS DECRESING SD --------------------------------------
dir_scenario <- "ain1_upmean_downsd"

desired_means <- c(30, 35, 40, 50)
desired_sds <- c(50, 40, 35, 30)
pars <- desired_mean_sd(mu_x = desired_means, sd_x = desired_sds)
A5 <- mixtures_plot_fx(pars = pars, dens_cols = dens_cols, pis = pis, desired_means = desired_means, desired_sds = desired_sds)

load(paste0("data/", dir_scenario, "/allin1df.RData"))
B5 <- bias_bxplot_fx(allin1df)

# PLOT ----------------------------------

nolegBs <- plot_grid(B1, B2, B3, B4, ncol = 1)
legend_Bs <- get_legend(B4+ theme(legend.position = "right")) 
Bs <- plot_grid(nolegBs, legend_Bs, ncol = 2, rel_widths = c(1, 0.2))

As <- plot_grid(A1, A2, A3, A4, ncol = 1, labels = "AUTO", label_size = 8)

plot_grid(As, Bs, ncol = 2, rel_widths = c(1, 1.5))

ggsave("Figures/Figure5.png", width = 6, height = 4, scale = 1.5, bg = "white")
