#-------------------------------------------------------------------------------
#                                 Functions
#-------------------------------------------------------------------------------

Simulate_Data <- function(nSubj, nGroups, dimVals, M, H, WM, WP, Noise, fileName) {
  # This function simulates augmented Gaussian gradients for two groups of N=nSubj
  # M, H, WM, WP, and Noise are vectors of length nGroups
  
  df <- list()
  
  for (g in 1:nGroups) {
    
    # generate individual parameters
    m <- rnorm(nSubj, M[g], .1)
    h = rnorm(nSubj, H[g], 1)
    wm = abs(rnorm(nSubj, WM[g], .1))
    wp = abs(rnorm(nSubj, WP[g], .1))
    noise = rep(Noise, nSubj)
    
    # generate individual gradients
    simData <- matrix(nrow = nSubj, ncol = length(dimVals))
    
    for (i in 1:nSubj) {
      gausLeft <- h[i] * exp(1)^-(((dimVals-m[i])^2) / (2 * wm[i]^2)) + 
        rnorm(length(dimVals), 0, noise)
      gausRight <- h[i] * exp(1)^-(((dimVals-m[i])^2) / (2 * wp[i]^2)) + 
        rnorm(length(dimVals), 0, noise)
      mIdx <- which.min(abs(dimVals-m[i]))
      if (mIdx == 1) {
        simData[i,] <- c(gausLeft[1], gausRight[(mIdx+1):length(dimVals)])
      } else {
        simData[i,] <- c(gausLeft[1:(mIdx-1)], gausRight[mIdx:length(dimVals)])
      }
    }
    subj <- rep(1:nSubj + (g-1)*nSubj, each = length(dimVals))
    x <- rep(dimVals, times = nSubj)
    y <- as.vector(t(simData))
    y[y > 100] <- 100
    y[y < 0] <- 0
    group <- rep(paste0("group", g), nSubj)
    df[[g]] <- data.frame(cbind(subj, group, x, y))
  }
  
  out <- do.call("rbind", df) 
  out <- out %>%
    mutate(x = as.numeric(as.character(x)),
           y = as.numeric(as.character(y)),
           subj = as.numeric(as.character(subj)))
  
  write_csv(out, paste0("data/", fileName ,".csv"))
}

#_______________________________________________________________________________

Read_Gen_Data <- function(fileName, dimVals, groupNames, xBreaks, xLab) {
  
  # This function reads data in long format and prepares the data list for stan
  # Notes
  # - Responses should be labelled as "y", subject ID as "subj", group names 
  # labelled as "group" (and should match groupNames)
  # - The dimVals argument does not have to match the labelling of the "x" column
  # - The CS+ should have a dimVal of 0
  # - Choose appropriate breaks and labels for the x-axis
  # See demo_data.csv for an example
  
  data <- read.csv(fileName, header = TRUE)
  
  out <- vector("list", length(groupNames))
  
  for (i in 1:length(groupNames)) {
    subset_data <- data %>% 
      filter(group == groupNames[i]) %>%
      arrange(subj, x)
    out[[i]] <- list(subj = subset_data[["subj"]],
                     responses = matrix(subset_data[["y"]], ncol = length(dimVals), 
                                        byrow = TRUE),
                     nSubj = n_distinct(subset_data[["subj"]]),
                     nStim = length(dimVals),
                     xs = dimVals)
  }
  names(out) <- groupNames
  
  layers <- list(
    geom_line(size = 1.5),
    geom_point(shape = rep(fig_shapes, each=length(dimVals)), fill = "white", size = 2),
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey"),
    scale_x_continuous(limits = c(min(dimVals), max(dimVals)), breaks = xBreaks, labels = xLab),
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)),
    scale_colour_manual(values = density_cols),
    labs(x = "stimulus", y = "mean response"),
    theme_classic(),
    guides(colour = guide_legend(override.aes = list(shape = fig_shapes)))
  )
  
  gradients <- data %>%
    arrange(subj, group, x) %>%
    mutate(group = fct_relevel(group, groupNames)) %>%
    group_by(group, x) %>%
    summarise(y = mean(y))
  
  fig <- ggplot(gradients, aes(x = x, y = y, group = group, colour = group, 
                               shape = group)) + layers
  
  return(list(out, fig))
}
#_______________________________________________________________________________

Run_Model <- function(dataList, modelName, groupName, modelFile, params) {
  
  # This function runs the model in stan
  
  stanfit <- stan(file = modelFile,
                  data = dataList, 
                  pars = params,
                  iter = n_iter, 
                  warmup = n_burnin, 
                  thin = n_thin, 
                  chains = n_chains, 
                  init = "random",
                  algorithm = "NUTS",
                  cores = parallel::detectCores())
  
  # save diagnostics
  diag <- rstan::get_sampler_params(stanfit, inc_warmup = FALSE)
  
  # save samples
  samples <- rstan::extract(stanfit)
  
  # save summary file
  summary <- rstan::summary(stanfit, probs = c(0.025, 0.50, 0.975))$summary
  write.csv(summary, file = paste0(file_name_root, modelName,"-", groupName, "-summary.csv"), row.names = TRUE)
  
  # waic and loo measures
  loglik <- extract_log_lik(stanfit)
  waic <- waic(loglik)
  loo <- loo(loglik)
  write_csv(as.data.frame(t(waic$estimates)), paste0(file_name_root, modelName, "-", groupName, "-waic.csv"))
  write_csv(as.data.frame(t(loo$estimates)), paste0(file_name_root, modelName, "-", groupName, "-loo.csv"))
  
  # output
  out <- list(stanfit, diag, samples, summary, waic, loo)
  names(out) <- c("stanfit", "diag", "samples", "summary", "waic", "loo")
  return(out)
}

#_______________________________________________________________________________

Plot_Densities <- function(samples, modelName, groupNames, paramNames) {
  
  # This function produces a multi-panelled figure with:
  # 1. Posterior distribution for the Mean parameter for each group
  # 2. Posterior distribution for the Height parameter for each group
  # 3. Posterior distribution for the Width- parameter for each group
  # 4. Posterior distribution for the Width+ parameter for each group
  # Note: length of samples and groupNames must match
  
  estimates <- data.frame()
  lengths <- rep(NA, length(groupNames))
  
  for (i in 1:length(groupNames)) {
    temp <- cbind(as.vector(samples[[i]]$M_group), 
                  as.vector(samples[[i]]$SDPlus_group), 
                  as.vector(samples[[i]]$SDMinus_group), 
                  as.vector(samples[[i]]$height_group))
    estimates <- rbind(estimates, temp)
    lengths[i] <- length(as.vector(samples[[i]]$M_group))
  }
  estimates$group <- as.factor(rep(groupNames, times = lengths))
  estimates <- estimates %>%
    mutate(group = fct_relevel(group, groupNames))
  colnames(estimates) <- c(paramNames, "group")
  
  density_layers <- list(geom_density(alpha = .5),
                         scale_fill_manual(values = density_cols),
                         theme_classic(),
                         theme(axis.title.x = element_blank()))
  
  M_fig <- ggplot(estimates, aes(M_group, fill = group)) + 
    density_layers +
    ggtitle("Mean")
  
  height_fig <- ggplot(estimates, aes(height_group, fill = group)) + 
    density_layers +
    ggtitle("Height")
  
  SDMinus_fig <- ggplot(estimates, aes(SDMinus_group, fill = group)) + 
    density_layers +
    ggtitle("Width-")
  
  SDPlus_fig <- ggplot(estimates, aes(SDPlus_group, fill = group)) + 
    density_layers +
    ggtitle("Width+")
  
  fig_panel <- M_fig + height_fig + SDMinus_fig + SDPlus_fig +
    plot_layout(nrow = 2, byrow = TRUE) + plot_layout(guides = "collect")
  return(fig_panel)
}

#_______________________________________________________________________________

Get_HDIs <- function(samples, modelName, groupName, HDIparams, ropeLow, ropeHigh, 
                     hdiLim, propPost = 1) {
  
  # This function calculates the Highest Density Intervals and 
  # p(Region of Practical Equivalence) for HDIparams
  # Note:
  # - ROPE parameters must be of length(HDIparams) and specified in the 
  # same order as HDIparams
  # - hdiLim sets the limit for the HDI
  # - propPost refers to the proportion of the posterior used to calculate p(direction) and p(ROPE)
  
  temp <- data.frame(param = HDIparams,
                     hdi_lim = hdiLim,
                     hdi_low = rep(NA, length(HDIparams)),
                     hdi_high = rep(NA, length(HDIparams)),
                     p_dir = rep(NA, length(HDIparams)),
                     rope_low = ropeLow,
                     rope_high = ropeHigh,
                     prop_rope = rep(NA, length(HDIparams)))
  for (i in 1:length(HDIparams)) {
    temp$hdi_low[i] <- hdi(as.vector(samples[[HDIparams[i]]]), ci = hdiLim)$CI_low
    temp$hdi_high[i] <- hdi(as.vector(samples[[HDIparams[i]]]), ci = hdiLim)$CI_high
    temp$p_dir[i] <- p_direction(as.vector(samples[[HDIparams[i]]]), ci = propPost)
    temp$prop_rope[i] <- rope(as.vector(samples[[HDIparams[i]]]), ci = propPost,
                              range = c(ropeLow[i], ropeHigh[i]))$ROPE_Percentage
  }
  
  write_csv(format(temp, scientific = FALSE), paste0(file_name_root, modelName, "-", groupName, "-hdis.csv"))
  
}

#_______________________________________________________________________________

Get_HDIs_diff <- function(samples_1, samples_2, modelName, comparison,
                          HDIparams, ropeLowDiffs, ropeHighDiffs, hdiLim, propPost = 1) {
  
  # This function calculates the Highest Density Intervals and 
  # p(Region of Practical Equivalence) for HDIparams for differences between 
  # samples_1 and samples_2
  # Note:
  # - ROPE parameters must be of length(HDIparams) and specified in the 
  # same order as HDIparams
  # - hdiLim sets the limit for the HDI
  # - propPost refers to the proportion of the posterior used to calculate p(direction) and p(ROPE)
  
  diff_m <- as.vector(samples_1$M_group) - c(as.vector(samples_2$M_group))
  diff_wplus <- as.vector(samples_1$SDPlus_group) - c(as.vector(samples_2$SDPlus_group))
  diff_wminus <- as.vector(samples_1$SDMinus_group) - c(as.vector(samples_2$SDMinus_group))
  diff_h <- as.vector(samples_1$height_group) - c(as.vector(samples_2$height_group))
  sample_diffs <- as.data.frame(cbind(diff_m, diff_wplus, diff_wminus, diff_h))
  colnames(sample_diffs) <- HDIparams
  temp <- data.frame(param = HDIparams,
                     hdi_lim = hdiLim,
                     hdi_low = rep(NA, length(HDIparams)),
                     hdi_high = rep(NA, length(HDIparams)),
                     p_dir = rep(NA, length(HDIparams)),
                     rope_low = ropeLowDiffs,
                     rope_high = ropeHighDiffs, 
                     prop_rope = rep(NA, length(HDIparams)))
  for (i in 1:length(HDIparams)) {
    temp$hdi_low[i] <- hdi(as.vector(sample_diffs[[HDIparams[i]]]), ci = hdiLim)$CI_low
    temp$hdi_high[i] <- hdi(as.vector(sample_diffs[[HDIparams[i]]]), ci = hdiLim)$CI_high
    temp$p_dir[i] <- p_direction(as.vector(sample_diffs[[HDIparams[i]]]), ci = propPost)
    temp$prop_rope[i] <- rope(as.vector(sample_diffs[[HDIparams[i]]]), ci = propPost,
                              range = c(temp$rope_low[i], temp$rope_high[i]))$ROPE_Percentage
  }
  write_csv(format(temp, scientific = FALSE), 
            paste0(file_name_root, modelName, "-", comparison, "-hdis.csv"))
  return(sample_diffs)
}

#_______________________________________________________________________________

Posterior_Preds <- function(samples, responses, modelName, groupName, dimVals, 
                            nSubj, subjList, nSamp, nStim, summary, 
                            nRow = nRow, figMult, labels = TRUE) {
  
  # This function plots the posterior predictives for each subject overlayed on 
  # the empirical gradients
  # - Optional labels: display the mean of the posterior for each parameter
  
  cs_loc <- (nStim-1)/2 # only if CS+ is in the middle of the dimension
  
  post_preds <- matrix(NA, nrow = nSubj*nStim*nSamp, ncol = 4)
  post_preds <- as.data.frame(post_preds)
  colnames(post_preds) <- c("subj", "dim", "samp", "pred")
  post_preds$subj <- rep(subjList, each = nSamp*nStim)
  post_preds$dim <- rep(rep(1:nStim, each = nSamp), times = nSubj)
  post_preds$samp <- rep(1:nSamp, times = nSubj*nStim)
  for (subj in 1:nSubj) {
    for (dim in 1:nStim) {
      start <- (subj-1) * nStim * nSamp+(dim-1) * nSamp + 1
      end <- (subj-1) * nStim * nSamp + (dim-1) * nSamp + nSamp
      post_preds$pred[start:end] <- sample(samples[["predR"]][, subj, dim], size = nSamp)
    }
  }
  
  # add responses
  post_preds$response <- NA
  post_preds$response <- rep(responses, each = nSamp)
  
  # add mean posterior estimates
  label <- rep(NA, nSubj)
  for (i in 1:nSubj) {
    label[i] <- paste0("M: ", round(summary$mean[i], 2),
                       " W-: ", round(summary$mean[nSubj*2 + i], 2),
                       "\n",
                       " W+: ", round(summary$mean[nSubj + i], 2),
                       " H: ", round(summary$mean[nSubj*3 + i]))
  }
  post_preds$label <- rep(label, each = nSamp * nStim)
  
  # figure layers
  grad_layers <- list(
    geom_line(stat = "identity", size = 1.25),
    labs(title = "", x = "dimension", y = "responding"),
    scale_colour_manual(values = fig_cols),
    geom_vline(xintercept = cs_loc, linetype = "dotted", colour = "black"), 
    scale_x_continuous(limits = c(1, nStim), breaks = c(1, cs_loc, nStim),
                       labels = c(min(dimVals), 0, max(dimVals))),
    scale_y_continuous(limits = c(0, 140), breaks = c(0, 50, 100)),
    theme_classic(),
    theme(panel.background = element_rect(colour = "black", size = 0.5, linetype = "solid", fill = NA),
          legend.position = "none")
  )
  
  # plot gradients alone
  grad_fig <- ggplot(post_preds, aes(y = response, x = dim)) + 
    grad_layers
  if (labels == TRUE) {
    fig <- grad_fig + 
      geom_text(data = post_preds, mapping = aes(label = label, x = cs_loc, y = 120), size = 2.5, colour = "black")
  } else {
    fig <- grad_fig
  }
  ggsave(file = paste0(file_name_root, modelName, "-", groupName, "-gradients", graph_file_type), 
         plot = fig + facet_wrap(~ subj, nrow = nRow),
         width = gg_width*figMult, height = gg_height*figMult, 
         units = "cm", dpi = dpi)
  
  pred_fig <- grad_fig + 
    geom_point(aes(y = pred, x = dim), colour = scat_col, size = scat_size, shape = scat_shape)
  if (labels == TRUE) {
    fig <- pred_fig + 
      geom_text(data = post_preds, mapping = aes(label = label, x = cs_loc, y = 120), size = 2.5, colour = "black")
  } else {
    fig <- pred_fig
  }
  ggsave(file = paste0(file_name_root, modelName, "-", groupName, "-postpreds", graph_file_type), 
         plot = fig + facet_wrap(~ subj, nrow = nRow),
         width = gg_width*figMult, height = gg_height*figMult, 
         units = "cm", dpi = dpi)
}


#_______________________________________________________________________________
Fit_Aug_Gaussian <- function(fileName, modelFile, modelName, groupNames, dimVals, params, 
                             HDIparams, ropeLow, ropeHigh, ropeLowDiffs, ropeHighDiffs,
                             hdiLim = .95, propPost = 1, 
                             nSamp, xBreaks, xLab, nRow, figMult = 2, labels = TRUE) {
  
  # Master function that reads the data file, runs the analysis, and saves output
  
  # 1. read data
  data_list <- Read_Gen_Data(
    fileName, dimVals, groupNames, xBreaks, xLab
  )
  
  mcmc_out <- vector("list", length(groupNames))
  samples <- vector("list", length(groupNames))
  
  for (i in 1:length(groupNames)) {
    
    # 2. fit models for each group
    mcmc_out[[i]] <- Run_Model(data_list[[1]][[i]], modelName, groupNames[i], modelFile, params)
    samples[[i]] <- mcmc_out[[i]]$samples
    
    # 3. HDIS and posterior summary stats for the group parameters 
    Get_HDIs(samples[[i]], modelName, groupNames[i], HDIparams, ropeLow, ropeHigh, hdi_limit)
    
    # 4. plot posterior predictives
    Posterior_Preds(samples[[i]], as.vector(t(data_list[[1]][[i]]$responses)), 
                    modelName, groupNames[i], dimVals,
                    data_list[[1]][[i]]$nSubj, unique(data_list[[1]][[i]]$subj),
                    nSamp, data_list[[1]][[i]]$nStim, as.data.frame(mcmc_out[[i]]$summary), 
                    nRow, figMult)
  }
  
  # 5. HDIs and posterior summary stats for group differences
  if (length(groupNames) == 2) {
    samples_diffs <- Get_HDIs_diff(samples[[1]], samples[[2]], modelName, paste0(groupNames[1], groupNames[2]),
                  HDIparams, ropeLowDiffs, ropeHighDiffs, hdi_limit)
  } else if (length(groupNames) == 3) {
    samples_diffs <- list()
    for (i in 1:length(groupNames)) {
      idx2 <- ifelse(i==3, 1, i+1)
      samples_diffs[[i]] <- Get_HDIs_diff(samples[[i]], samples[[idx2]], modelName, 
                                          paste0(groupNames[i], "V", groupNames[idx2]),
                                          HDIparams, ropeLowDiffs, ropeHighDiffs, hdi_limit)
    }
  } else {
    samples_diffs <- NA
    print("Call Get_HDIs_diff function for each group comparison")
  }
  
  # 6. plot gradients + posteriors
  density_fig <- Plot_Densities(samples, modelName, groupNames, HDIparams)
  fig_panel <- data_list[[2]] + density_fig + plot_layout(widths = c(1, 1.75), guides = "collect")
  ggsave(paste0(file_name_root, modelName, "-density", graph_file_type), fig_panel, 
         "jpeg", height = gg_height, width = gg_width*2.5, units = "cm", dpi = dpi)
  
  # return output list
  out <- list(data_list[[1]], data_list[[2]], mcmc_out, samples, density_fig, samples_diffs)
  names(out) <- c("data_lists", "gradients", "mcmc_out", "samples", "density_fig", "samples_diffs")
  return(out)
}

#_______________________________________________________________________________