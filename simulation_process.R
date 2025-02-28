
## Simulate data -------------------------------------------------------

# Using purrr to get samples instead of a for loop. More efficient.

# We are getting 50k samples from the mixture
# These are the weights for the distribution
pis <- c( 0.1, 0.2, 0.3, 0.4)
tru_n <- 50000
w_samp_sizes <- pis*tru_n
# Using purrr to pull values from each component of the mixture according to the weights
# This instead of using a for loop
purrrsampls <- tibble(gID = c(1:length(w_samp_sizes)), pars) %>%
   mutate(x_samps = purrr::map(gID, function(y) rlnorm(w_samp_sizes[y], pars$meanlog[y], pars$sdlog[y])))

# Making it an easier to read data frame with only the group ID (mixture components) and the sampled data
purrrsampls %>%
   dplyr::select(., gID, x_samps) %>%
   unnest(cols = c(x_samps)) -> truth_df

save(truth_df, file = paste0("data/", dir_scenario, "/mixturesamples.RData"))

# This is the data that we consider our "truth"
tru_n <- nrow(truth_df)

thresh_tests <- c(50, 100, 150, 250, 500, 750, 1000)

tibble(thresh_tests) %>%
   mutate(n_overthresh = map_dbl(1:length(thresh_tests),
                                 function(y) length(which(truth_df$x_samps >= thresh_tests[y]))),
          theta = n_overthresh/tru_n) -> thetas
thetas


# Estimation --------------------------------------------------------------

# This section gets the estimation for the thetas, the tails
# using the lomax, glm lomax, and the pareto with a threshold of 0 (GEV)
# nreps_mles_df <- bxplt_data_fx(points_for_boxplots) # can use function but changing threshold below
# save(nreps_mles_df, file = paste0("data/", dir_scenario, "/nreps_mles_df.RData"))

nreps <- points_for_boxplots
out <- data.frame()

for(j in 1:nreps){
   # tests for the number of samples
   samp_n_tests <- c(80, 200, 500, 800, 1000, 1600)
   # combinations for number of samples and estimating beyond a given threshold
   mles_df <- data.frame(expand.grid(thresh_tests = thresh_tests, samp_n_tests = samp_n_tests))
   mles_df$nrep <- j
   
   for(i in 1:nrow(mles_df)){
      ith_n <- mles_df$samp_n_tests[i]
      ith_thresh <- mles_df$thresh_tests[i]
      ith_samps <- data.frame(x_star = sample(truth_df$x_samps, ith_n))
      ith_tail <- length(which(ith_samps$x_star >= ith_thresh))
      mles_df$tail_p[i] <- ith_tail/ith_n
      
      # Fit Lomax
      optim_out <- optim(par=log(c(1.5, 1.5)), fn=nllike.simp, method="BFGS", Y=ith_samps$x_star)
      mles_star <- exp(optim_out$par)
      alpha_star <- mles_star[1]
      k_star <- mles_star[2]
      mles_df$alpha_star[i] <-alpha_star
      mles_df$k_star[i] <- k_star
      
      # Check with Lomax GLM
      glm_out <- lomax.glm(formula=~1, my.dataf=ith_samps, response=ith_samps$x_star)
      alpha.2 <- glm_out$alphas.hat[1]
      k.2     <- glm_out$k.hat
      mles_df$alpha_glm[i] <-alpha.2
      mles_df$k_glm[i] <- k.2
      
      # Check with GP fit
      ith_evd <- fevd(ith_samps$x_star, threshold = median(ith_samps$x_star), type = "GP")
      gp_scale <- summary(ith_evd)$par[1]
      gp_shape <- summary(ith_evd)$par[2]
      mles_df$alpha_GP[i] <- gp_scale*gp_shape
      mles_df$k_GP[i] <- 1/gp_shape
      mles_df$scale[i] <- gp_scale
      mles_df$shape[i] <- gp_shape
      
      
      
      # Estimate the tail
      #Lomax simple
      log_St_star <- lomax.St(x = ith_thresh,alpha = alpha_star,k = k_star,log.scale=TRUE)
      mles_df$log_St_hat[i] <- log_St_star
      
      # Lomax glm
      log.st.star2 <- lomax.St(x=ith_thresh,alpha=alpha.2,k=k.2,log.scale=TRUE)
      mles_df$log_St_hat2[i] <- log.st.star2
      
      # GP tail
      # mles_df$log_GP[i] <- log(pextRemes(ith_evd, ith_thresh, lower.tail = FALSE))
      mles_df$GP_theta_hat[i] <- pextRemes(ith_evd, ith_thresh, lower.tail = FALSE)
      
      
   }
   
   mles_df %>%
      right_join(., thetas[, c(1,3)]) -> mles_df
   out <- rbind.data.frame(out, mles_df)
}

nreps_mles_df <- out

save(nreps_mles_df, file = paste0("data/", dir_scenario, "/nreps_mles_df.RData"))

## GP threshold choices --------------------------------
# Here we incorporate the threshold
# Compare to a fit using a Generalized Pareto with two different approaches
# One approach sets the threshold at the median of the sample
# The other sets the threshold based on the SD of a nonparametric bootstrap of the data

out <- data.frame()

for(j in 1:points_for_boxplots){
   samp_n_tests <- c(80, 200, 500, 800, 1000, 1600)
   mles_df <- data.frame(expand.grid(thresh_tests = thresh_tests, samp_n_tests = samp_n_tests))
   mles_df$nrep <- j
   
   for(i in 1:nrow(mles_df)){
      ith_n <- mles_df$samp_n_tests[i]
      ith_thresh <- mles_df$thresh_tests[i]
      ith_samps <- data.frame(x_star = sample(truth_df$x_samps, ith_n))
      
      
      # nonparam bootstrap to define the threshold
      k_out <- NULL
      for(k in 1:1000){
         kth_samps <- data.frame(x_star = sample(ith_samps$x_star, ith_n, replace = TRUE))
         # k_quantiles <- round(as.numeric(quantile(kth_samps$x_star, 0.3)))
         k_sd <- quantile(kth_samps$x_star, 0.5)
         
         k_out <- rbind(k_out, k_sd)
      }
      kth_thresh <- max(k_out)
      
      ith_tail <- length(which(ith_samps$x_star >= ith_thresh))
      mles_df$tail_p[i] <- ith_tail/ith_n
      
      ith_quant <- round(as.numeric(quantile(ith_samps$x_star, 0.5)))
      mles_df$kth_thresh[i] <- kth_thresh
      mles_df$ith_quant[i] <- ith_quant
      
      
      # GP fit and tail based on the quantile - median
      ith_evd <- fevd(ith_samps$x_star, threshold = kth_thresh, type = "GP")
      mles_df$quant_theta[i] <- pextRemes(ith_evd, ith_thresh, lower.tail = FALSE)
      
      # GP fit and tail based on the standard deviation
      kth_evd <- fevd(ith_samps$x_star, threshold = kth_thresh, type = "GP")
      mles_df$sd_theta[i] <- pextRemes(kth_evd, kth_thresh, lower.tail = FALSE)
      
   }
   
   mles_df %>%
      right_join(., thetas[, c(1,3)]) -> mles_df
   out <- rbind.data.frame(out, mles_df)
}

GP_mles <- out
save(GP_mles, file = paste0("data/", dir_scenario, "/GP_runs.RData"))



# Bias correction ---------------------------------------------------------

## Sample sizes --------------------
samp_n_tests <- c(80, 200, 500, 800, 1000, 1600)

all_test_combos <- data.frame(expand.grid(thresh_tests = thresh_tests, samp_n_tests = samp_n_tests))


nruns <- paste0("nrun_", 1:points_for_boxplots)
combo_tests <- data.frame(expand.grid(thresh_tests = thresh_tests, samp_n_tests = samp_n_tests, nruns = nruns))
allin1df <- combo_tests
boots_df <- NULL

for(i in 1:nrow(combo_tests)){
   # Get the sample
   ith_n <- combo_tests$samp_n_tests[i]
   ith_thresh <- combo_tests$thresh_tests[i]
   ith_samps <- data.frame(x_star = sample(truth_df$x_samps, ith_n))
   
   # Lomax GLM
   ith_glm <- lomax.glm(formula=~1, my.dataf=ith_samps, response=ith_samps$x_star)
   ith_alpha_hat <- ith_glm$alphas.hat[1]
   ith_k_hat <- ith_glm$k.hat
   
   # Estimate tail
   ith_theta_hat <- lomax.St(x = ith_thresh, alpha = ith_alpha_hat, k = ith_k_hat, log.scale=TRUE)
   allin1df$lomax_theta_hat[i] <- ith_theta_hat
   
   # Generalized Pareto
   ith_quant <- quantile(ith_samps$x_star, 0.5)
   allin1df$sampl_quant[i] <- ith_quant
   ith_evd <- fevd(ith_samps$x_star, threshold = ith_quant, type = "GP")
   allin1df$gp_theta_hat[i] <- pextRemes(ith_evd, ith_thresh, lower.tail = FALSE)
   
   # Nonparametric bootstrapping
   bth_df <- data.frame(n_B = 1:B)
   
   for(b in 1:B){
      # Sample with replacement from the given sample
      bth_samps <- data.frame(x_star = sample(ith_samps$x_star, ith_n, replace = TRUE))
      
      # Estimate params using glm Lomax
      bth_glm <- lomax.glm(formula=~1, my.dataf=bth_samps, response=bth_samps$x_star)
      bth_alpha_star <- bth_glm$alphas.hat[1]
      bth_k_star <- bth_glm$k.hat
      bth_df$lomax_theta_star[b] <- lomax.St(x = ith_thresh, alpha = bth_alpha_star, k = bth_k_star, log.scale=TRUE)
      
      # Generalized Pareto with threshold set to the median of the sample
      bth_quant <- quantile(bth_samps$x_star, 0.5)
      bth_evd <- fevd(bth_samps$x_star, threshold = bth_quant, type = "GP")
      bth_df$quant_star[b] <- bth_quant
      bth_df$gp_theta_star[b] <- pextRemes(bth_evd, ith_thresh, lower.tail = FALSE)
      
      bth_df$ith_row[b] <- paste(i)
      bth_df$ith_n[b] <- ith_n
      bth_df$ith_thresh <- ith_thresh
      
   }
   
   boots_df[[i]] <- bth_df
   
   allin1df$lomax_theta_stars[i] <- mean(bth_df$lomax_theta_star)
   allin1df$gp_theta_stars[i] <- mean(bth_df$gp_theta_star)
   
}

allin1df %>%
   right_join(., thetas[,c(1,3)]) %>%
   mutate(lomax_theta_bar = 2*lomax_theta_hat - lomax_theta_stars,
          gp_theta_bar = 2*gp_theta_hat - gp_theta_stars,
          lomax_ratio = lomax_theta_bar - log(theta),
          gp_ratio = log(gp_theta_bar/theta)) -> allin1df

save(allin1df, file = paste0("data/", dir_scenario, "/allin1df.RData"))



