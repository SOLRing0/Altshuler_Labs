
################################### Packages ###################################
## Specify the packages of interest
packages <-
  c("ggplot2",
    "tidyverse",
    "tidyr",
    "dplyr",
    "tibble",
    "readr",
    "readxl",
    "ggthemes",
    "magick",
    "grid",
    "cowplot",
    "metR",
    "viridis",
    "R.matlab",
    "gaussplotR")

## Load or install & load all packages
package.check <- lapply(
  packages,
  FUN = function(x)
  {
    if (!require(x, character.only = TRUE))
    {
      install.packages(x, dependencies = TRUE,
                       repos = "http://cran.us.r-project.org")
      library(x, character.only = TRUE)
    }
  }
)


#### import truncated data ####
lm_initrns_dir <-
  "./data/spatemp/LM_ER_transients_200bins_bin004_to_bin020/"

lm_initrns_files <-
  list.files(lm_initrns_dir, recursive = TRUE)

lm_initrns_filepaths <-
  paste0(
    lm_initrns_dir,
    lm_initrns_files
  )

lm_initrns_basenames <-
  basename(lm_initrns_files) %>%
  str_remove("_transients.csv")

## extract the relevant stuff
lm_initrns_dat <-
  map(lm_initrns_filepaths, read_csv, col_names = TRUE) 
names(lm_initrns_dat) <- lm_initrns_basenames


## cell1_norms
lm_initrns_cell1_norms <-
  lm_initrns_dat %>%
  map(select, logSF, logTF, ends_with("cell1_norm")) %>%
  map(drop_na)
vec <- NULL
for (i in 1:length(lm_initrns_cell1_norms)) {
  vec[i] <- ncol(lm_initrns_cell1_norms[[i]])
}
lm_initrns_cell1_norms <-lm_initrns_cell1_norms[vec > 2]
names(lm_initrns_cell1_norms) <-
  paste0(names(lm_initrns_cell1_norms),
         "-cell1")

## cell2_norms
lm_initrns_cell2_norms <-
  lm_initrns_dat %>%
  map(select, logSF, logTF, ends_with("cell2_norm")) %>%
  map(drop_na)
vec <- NULL
for (i in 1:length(lm_initrns_cell2_norms)) {
  vec[i] <- ncol(lm_initrns_cell2_norms[[i]])
}
lm_initrns_cell2_norms <-lm_initrns_cell2_norms[vec > 2]
names(lm_initrns_cell2_norms) <-
  paste0(names(lm_initrns_cell2_norms),
         "-cell2")

## cell3_norms
lm_initrns_cell3_norms <-
  lm_initrns_dat %>%
  map(select, logSF, logTF, ends_with("cell3_norm")) 
vec <- NULL
for (i in 1:length(lm_initrns_cell3_norms)) {
  vec[i] <- ncol(lm_initrns_cell3_norms[[i]])
}
lm_initrns_cell3_norms <-lm_initrns_cell3_norms[vec > 2]
names(lm_initrns_cell3_norms) <-
  paste0(names(lm_initrns_cell3_norms),
         "-cell3")

## cell4_norms
lm_initrns_cell4_norms <-
  lm_initrns_dat %>%
  map(select, logSF, logTF, ends_with("cell4_norm"))
vec <- NULL
for (i in 1:length(lm_initrns_cell4_norms)) {
  vec[i] <- ncol(lm_initrns_cell4_norms[[i]])
}
lm_initrns_cell4_norms <-lm_initrns_cell4_norms[vec > 2]
names(lm_initrns_cell4_norms) <-
  paste0(names(lm_initrns_cell4_norms),
         "-cell4")

#### get all the data ####
lm_initrns_all_cell_norms <-
  c(lm_initrns_cell1_norms,
    lm_initrns_cell2_norms,
    lm_initrns_cell3_norms,
    lm_initrns_cell4_norms)  %>%
  .[order(names(.))]


#### bring in truncation info ####
ER_initrns_truncs <-
  read_csv("./analyses_and_figures/figure 4/ER_initrns_truncs.csv")

#### truncation step ####
ER_initrns_contour_dat <- NULL
for (i in 1:length(lm_initrns_all_cell_norms)) {
  print(i)
  
  ## extract data
  dat <-
    lm_initrns_all_cell_norms[[i]] %>%
    drop_na() %>%
    as.data.frame()
  colnames(dat) <-
    c("X_values", "Y_values", "response")
  
  ## truncate if needed
  trnk <-
    ER_initrns_truncs %>%
    filter(name == names(lm_initrns_all_cell_norms)[i])
  
  datj <-
    dat %>%
    filter(X_values >= as.numeric(trnk[1, "min_sf"]) - 0.1) %>%
    filter(X_values <= as.numeric(trnk[1, "max_sf"]) + 0.1) %>%
    filter(Y_values >= as.numeric(trnk[1, "min_tf"]) - 0.1) %>%
    filter(Y_values <= as.numeric(trnk[1, "max_tf"]) + 0.1)
  
  ER_initrns_contour_dat[[i]] <- datj
  names(ER_initrns_contour_dat)[[i]] <-
    names(lm_initrns_all_cell_norms)[[i]]
  
}

#### fit unconstrained gaussians ####
ER_initrns_unconst_gauss_params <- NULL
ER_initrns_G_result_un <- NULL
ER_initrns_unconst_gauss_predicts <- NULL
for (i in 1:length(lm_initrns_all_cell_norms)) {
  print(i)
  
  tryCatch({
    
    ER_initrns_unconst_gauss_params[[i]] <-
      gaussplotR::fit_gaussian_2D(
        data = ER_initrns_contour_dat[[i]],
        method = "elliptical_log",
        constrain_amplitude = TRUE,
        constrain_orientation = "unconstrained", 
        maxiter = 10000,
        minFactor = 0.000000488281
      )
    names(ER_initrns_unconst_gauss_params)[[i]] <-
      names(lm_initrns_all_cell_norms)[[i]]
    
    ER_initrns_G_result_un[[i]] <-
      gaussplotR::predict_gaussian_2D(
        fit_object = ER_initrns_unconst_gauss_params[[i]],
        X_values = ER_initrns_contour_dat[[i]]$X_values,
        Y_values = ER_initrns_contour_dat[[i]]$Y_values
      )
    names(ER_initrns_G_result_un)[[i]] <-
      names(lm_initrns_all_cell_norms)[[i]]
    
    grid <-
      expand.grid(
        X_values = seq(
          from = round(min(ER_initrns_contour_dat[[i]]$X_values), 1),
          to = round(max(ER_initrns_contour_dat[[i]]$X_values), 1),
          by = 0.1
        ),
        Y_values = seq(
          from = round(min(ER_initrns_contour_dat[[i]]$Y_values), 1),
          to = round(max(ER_initrns_contour_dat[[i]]$Y_values), 1),
          by = 0.1
        )
      )
    
    ER_initrns_unconst_gauss_predicts[[i]] <-
      gaussplotR::predict_gaussian_2D(
        fit_object = ER_initrns_unconst_gauss_params[[i]],
        X_values = grid$X_values,
        Y_values = grid$Y_values
      )
    names(ER_initrns_unconst_gauss_predicts)[[i]] <-
      names(lm_initrns_all_cell_norms)[[i]]
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

ER_initrns_unconst_gauss_params[[147]] <- list(NULL)
ER_initrns_G_result_un[[147]] <- list(NULL)
ER_initrns_unconst_gauss_predicts[[147]] <- list(NULL)

## collect params
ER_initrns_unconst_gauss_paramdf <- NULL
for (i in 1:length(ER_initrns_unconst_gauss_params)) {
  tmp <-
    ER_initrns_unconst_gauss_params[[i]]$coefs
  tmp$name <- names(lm_initrns_all_cell_norms)[[i]]
  ER_initrns_unconst_gauss_paramdf <-
    bind_rows(
      ER_initrns_unconst_gauss_paramdf,
      tmp
    )
}

## collect summary stats
ER_initrns_unconst_gauss_fitstatsdf <- NULL
for (i in 1:length(ER_initrns_unconst_gauss_params)) {
  if (is.null(ER_initrns_unconst_gauss_params[[i]]$model_error_stats)) {
    tmp <- data.frame(
      rss = NA,
      rmse = NA,
      deviance = NA,
      AIC = NA,
      R2 = NA,
      R2_adj = NA,
      name = names(lm_initrns_all_cell_norms)[[i]]
    )
  } else{
    tmp <-
      ER_initrns_unconst_gauss_params[[i]]$model_error_stats
    tmp$name <- names(ER_initrns_unconst_gauss_params)[i]
  }
  if (is.null(ER_initrns_unconst_gauss_params[[i]]$fit_method)) {
    tmp2 <- data.frame(
      method = NA,
      amplitude = NA,
      orientation = NA,
      name = names(lm_initrns_all_cell_norms)[[i]]
    )
  } else{
    tmp2 <-
      as.data.frame(t(ER_initrns_unconst_gauss_params[[i]]$fit_method))
    tmp2$name <- names(ER_initrns_unconst_gauss_params)[i]
  }
  tmp3 <-
    left_join(tmp2, tmp, by = "name") %>%
    relocate("name")
  ER_initrns_unconst_gauss_fitstatsdf <-
    bind_rows(
      ER_initrns_unconst_gauss_fitstatsdf,
      tmp3
    )
}

## put them together
ER_initrns_unconst_gauss_summary <-
  left_join(ER_initrns_unconst_gauss_paramdf,
            ER_initrns_unconst_gauss_fitstatsdf) %>%
  relocate("name")


#### fit Q = -1 indep gaussians ####
ER_initrns_indep_gauss_params <- NULL
ER_initrns_G_result_ind <- NULL
ER_initrns_indep_gauss_predicts <- NULL
for (i in 1:length(lm_initrns_all_cell_norms)) {
  print(i)
  
  tryCatch({
    
    ER_initrns_indep_gauss_params[[i]] <-
      gaussplotR::fit_gaussian_2D(
        data = ER_initrns_contour_dat[[i]],
        method = "elliptical_log",
        constrain_amplitude = TRUE,
        constrain_orientation = -1, 
        maxiter = 10000,
        minFactor = 0.000000488281
      )
    names(ER_initrns_indep_gauss_params)[[i]] <-
      names(lm_initrns_all_cell_norms)[[i]]
    
    ER_initrns_G_result_ind[[i]] <-
      gaussplotR::predict_gaussian_2D(
        fit_object = ER_initrns_indep_gauss_params[[i]],
        X_values = ER_initrns_contour_dat[[i]]$X_values,
        Y_values = ER_initrns_contour_dat[[i]]$Y_values
      )
    names(ER_initrns_G_result_ind)[[i]] <-
      names(lm_initrns_all_cell_norms)[[i]]
    
    grid <-
      expand.grid(
        X_values = seq(
          from = round(min(ER_initrns_contour_dat[[i]]$X_values), 1),
          to = round(max(ER_initrns_contour_dat[[i]]$X_values), 1),
          by = 0.1
        ),
        Y_values = seq(
          from = round(min(ER_initrns_contour_dat[[i]]$Y_values), 1),
          to = round(max(ER_initrns_contour_dat[[i]]$Y_values), 1),
          by = 0.1
        )
      )
    
    ER_initrns_indep_gauss_predicts[[i]] <-
      gaussplotR::predict_gaussian_2D(
        fit_object = ER_initrns_indep_gauss_params[[i]],
        X_values = grid$X_values,
        Y_values = grid$Y_values
      )
    names(ER_initrns_indep_gauss_predicts)[[i]] <-
      names(lm_initrns_all_cell_norms)[[i]]
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

ER_initrns_indep_gauss_params[[147]] <- list(NULL)
ER_initrns_G_result_ind[[147]] <- list(NULL)
ER_initrns_indep_gauss_predicts[[147]] <- list(NULL)

## collect params
ER_initrns_indep_gauss_paramdf <- NULL
for (i in 1:length(ER_initrns_indep_gauss_params)) {
  tmp <-
    ER_initrns_indep_gauss_params[[i]]$coefs
  tmp$name <- names(lm_initrns_all_cell_norms)[[i]]
  ER_initrns_indep_gauss_paramdf <-
    bind_rows(
      ER_initrns_indep_gauss_paramdf,
      tmp
    )
}

## collect summary stats
ER_initrns_indep_gauss_fitstatsdf <- NULL
for (i in 1:length(ER_initrns_indep_gauss_params)) {
  if (is.null(ER_initrns_indep_gauss_params[[i]]$model_error_stats)) {
    tmp <- data.frame(
      rss = NA,
      rmse = NA,
      deviance = NA,
      AIC = NA,
      R2 = NA,
      R2_adj = NA,
      name = names(lm_initrns_all_cell_norms)[[i]]
    )
  } else{
    tmp <-
      ER_initrns_indep_gauss_params[[i]]$model_error_stats
    tmp$name <- names(ER_initrns_indep_gauss_params)[i]
  }
  if (is.null(ER_initrns_indep_gauss_params[[i]]$fit_method)) {
    tmp2 <- data.frame(
      method = NA,
      amplitude = NA,
      orientation = NA,
      name = names(lm_initrns_all_cell_norms)[[i]]
    )
  } else{
    tmp2 <-
      as.data.frame(t(ER_initrns_indep_gauss_params[[i]]$fit_method))
    tmp2$name <- names(ER_initrns_indep_gauss_params)[i]
  }
  tmp3 <-
    left_join(tmp2, tmp, by = "name") %>%
    relocate("name")
  ER_initrns_indep_gauss_fitstatsdf <-
    bind_rows(
      ER_initrns_indep_gauss_fitstatsdf,
      tmp3
    )
}

## put them together
ER_initrns_indep_gauss_summary <-
  left_join(ER_initrns_indep_gauss_paramdf, ER_initrns_indep_gauss_fitstatsdf) %>%
  relocate("name")


#### fit Q = 0 spd gaussians ####
ER_initrns_speed_gauss_params <- NULL
ER_initrns_G_result_spd <- NULL
ER_initrns_speed_gauss_predicts <- NULL
for (i in 1:length(lm_initrns_all_cell_norms)) {
  print(i)
  
  tryCatch({
    
    ER_initrns_speed_gauss_params[[i]] <-
      gaussplotR::fit_gaussian_2D(
        data = ER_initrns_contour_dat[[i]],
        method = "elliptical_log",
        constrain_amplitude = TRUE,
        constrain_orientation = 0, 
        maxiter = 10000,
        minFactor = 0.000000488281
      )
    names(ER_initrns_speed_gauss_params)[[i]] <-
      names(lm_initrns_all_cell_norms)[[i]]
    
    ER_initrns_G_result_spd[[i]] <-
      gaussplotR::predict_gaussian_2D(
        fit_object = ER_initrns_speed_gauss_params[[i]],
        X_values = ER_initrns_contour_dat[[i]]$X_values,
        Y_values = ER_initrns_contour_dat[[i]]$Y_values
      )
    names(ER_initrns_G_result_spd)[[i]] <-
      names(lm_initrns_all_cell_norms)[[i]]
    
    grid <-
      expand.grid(
        X_values = seq(
          from = round(min(ER_initrns_contour_dat[[i]]$X_values), 1),
          to = round(max(ER_initrns_contour_dat[[i]]$X_values), 1),
          by = 0.1
        ),
        Y_values = seq(
          from = round(min(ER_initrns_contour_dat[[i]]$Y_values), 1),
          to = round(max(ER_initrns_contour_dat[[i]]$Y_values), 1),
          by = 0.1
        )
      )
    
    ER_initrns_speed_gauss_predicts[[i]] <-
      gaussplotR::predict_gaussian_2D(
        fit_object = ER_initrns_speed_gauss_params[[i]],
        X_values = grid$X_values,
        Y_values = grid$Y_values
      )
    names(ER_initrns_speed_gauss_predicts)[[i]] <-
      names(lm_initrns_all_cell_norms)[[i]]
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

## collect params
ER_initrns_speed_gauss_paramdf <- NULL
for (i in 1:length(ER_initrns_speed_gauss_params)) {
  tmp <-
    ER_initrns_speed_gauss_params[[i]]$coefs
  tmp$name <- names(lm_initrns_all_cell_norms)[[i]]
  ER_initrns_speed_gauss_paramdf <-
    bind_rows(
      ER_initrns_speed_gauss_paramdf,
      tmp
    )
}

## collect summary stats
ER_initrns_speed_gauss_fitstatsdf <- NULL
for (i in 1:length(ER_initrns_speed_gauss_params)) {
  if (is.null(ER_initrns_speed_gauss_params[[i]]$model_error_stats)) {
    tmp <- data.frame(
      rss = NA,
      rmse = NA,
      deviance = NA,
      AIC = NA,
      R2 = NA,
      R2_adj = NA,
      name = names(lm_initrns_all_cell_norms)[[i]]
    )
  } else{
    tmp <-
      ER_initrns_speed_gauss_params[[i]]$model_error_stats
    tmp$name <- names(ER_initrns_speed_gauss_params)[i]
  }
  if (is.null(ER_initrns_speed_gauss_params[[i]]$fit_method)) {
    tmp2 <- data.frame(
      method = NA,
      amplitude = NA,
      orientation = NA,
      name = names(lm_initrns_all_cell_norms)[[i]]
    )
  } else{
    tmp2 <-
      as.data.frame(t(ER_initrns_speed_gauss_params[[i]]$fit_method))
    tmp2$name <- names(ER_initrns_speed_gauss_params)[i]
  }
  tmp3 <-
    left_join(tmp2, tmp, by = "name") %>%
    relocate("name")
  ER_initrns_speed_gauss_fitstatsdf <-
    bind_rows(
      ER_initrns_speed_gauss_fitstatsdf,
      tmp3
    )
}

## put them together
ER_initrns_speed_gauss_summary <-
  left_join(ER_initrns_speed_gauss_paramdf,
            ER_initrns_speed_gauss_fitstatsdf) %>%
  relocate("name")

#### characterize to get z-diffs ####
ER_initrns_Z_results <- NULL
for (i in 1:length(ER_initrns_contour_dat)) {
  print(i)
  tryCatch({
    ER_initrns_Z_results[[i]] <-
      characterize_gaussian_fits(
        data = ER_initrns_contour_dat[[i]],
        constrain_amplitude = TRUE,
        maxiter = 10000,
        minFactor = 0.000000488281
      )
  }, error = function(e) {
    cat("ERROR :", conditionMessage(e), "\n")
  })
}

ER_initrns_Z_results[[147]] <- list(NULL)

## collect params
ER_initrns_Z_paramdf <- NULL
for (i in 1:length(ER_initrns_Z_results)) {
  print (i)
  if (nrow(as_tibble(ER_initrns_Z_results[[i]][3:10])) == 0){
    tmp <- matrix(nrow = 1, ncol = 8) %>% as_tibble()
  } else{
    tmp <-
      ER_initrns_Z_results[[i]][3:10] %>% as_tibble()
  }
  
  tmp$name <- names(ER_initrns_contour_dat)[[i]]
  ER_initrns_Z_paramdf <-
    bind_rows(
      ER_initrns_Z_paramdf,
      tmp
    )[,1:9]
}

#### ALL CELLS CONTOURS ####
ER_initrns_allcells_contour_plots <- NULL
for (i in 1:length(ER_initrns_contour_dat)) {
  ER_initrns_allcells_contour_plots[[i]] <-
    ER_initrns_contour_dat[[i]] %>%
    ggplot(aes(x = X_values,
               y = Y_values,
               z = response)) +
    metR::geom_contour_fill(aes(fill = ..level..),
                            color = "black",
                            size = 0.03,
                            bins = 15) +
    xlab("log2 of SF (cpd)") +
    ylab("log2 of TF (Hz)") +
    ggtitle(
      paste0(
        names(ER_initrns_contour_dat)[i],
        " contour")
    ) +
    coord_fixed(expand = FALSE) +
    theme_classic() +
    theme(plot.title = element_text(size = 8))
}

#### Unconstrained gaussian plots ####
ER_initrns_allcells_unconst_pred_plots <- NULL
for (i in 1:length(ER_initrns_unconst_gauss_predicts)) {
  if (is.null(ER_initrns_unconst_gauss_predicts[[i]][[1]])) {
    ER_initrns_allcells_unconst_pred_plots[[i]] <- NULL
  } else {
    ER_initrns_allcells_unconst_pred_plots[[i]] <-
      ER_initrns_unconst_gauss_predicts[[i]] %>%
      ggplot(aes(x = X_values,
                 y = Y_values,
                 z = predicted_values/max(predicted_values))) +
      metR::geom_contour_fill(aes(fill = ..level..),
                              color = "black",
                              size = 0.03,
                              bins = 15) +
      xlab("log2 of SF (cpd)") +
      ylab("log2 of TF (Hz)") +
      geom_text(
        x = range(ER_initrns_unconst_gauss_predicts[[i]]$X_values)[2] - 1,
        y = range(ER_initrns_unconst_gauss_predicts[[i]]$Y_values)[1] + 0.5,
        label = paste0(
          "R2 = ", round(ER_initrns_unconst_gauss_summary[i,]$R2_adj, 2)
        ),
        color = "white"
      ) +
      ggtitle(
        paste0(
          names(ER_initrns_unconst_gauss_predicts)[i],
          " Q = ",
          round(ER_initrns_unconst_gauss_summary[i,]$Q, 2))
      ) +
      coord_fixed(expand = FALSE) +
      theme_classic() +
      theme(plot.title = element_text(size = 8))
  }
}

#### Q = -1 gaussian plots ####
ER_initrns_allcells_indep_pred_plots <- NULL
for (i in 1:length(ER_initrns_indep_gauss_predicts)) {
  if (is.null(ER_initrns_indep_gauss_predicts[[i]][[1]])) {
    ER_initrns_allcells_indep_pred_plots[[i]] <- NULL
  } else {
    ER_initrns_allcells_indep_pred_plots[[i]] <-
      ER_initrns_indep_gauss_predicts[[i]] %>%
      ggplot(aes(x = X_values,
                 y = Y_values,
                 z = predicted_values/max(predicted_values))) +
      metR::geom_contour_fill(aes(fill = ..level..),
                              color = "black",
                              size = 0.03,
                              bins = 15) +
      xlab("log2 of SF (cpd)") +
      ylab("log2 of TF (Hz)") +
      geom_text(
        x = range(ER_initrns_indep_gauss_predicts[[i]]$X_values)[2] - 1,
        y = range(ER_initrns_indep_gauss_predicts[[i]]$Y_values)[1] + 0.5,
        label = paste0(
          "R2 = ", round(ER_initrns_indep_gauss_summary[i,]$R2_adj, 2)
        ),
        color = "white"
      ) +
      ggtitle(
        paste0(
          names(ER_initrns_indep_gauss_predicts)[i],
          " Q = -1")
      ) +
      coord_fixed(expand = FALSE) +
      theme_classic() +
      theme(plot.title = element_text(size = 8))
  }
}

#### Q = 0 gaussian plots ####
ER_initrns_allcells_speed_pred_plots <- NULL
for (i in 1:length(ER_initrns_speed_gauss_predicts)) {
  if (is.null(ER_initrns_speed_gauss_predicts[[i]][[1]])) {
    ER_initrns_allcells_speed_pred_plots[[i]] <- NULL
  } else {
    ER_initrns_allcells_speed_pred_plots[[i]] <-
      ER_initrns_speed_gauss_predicts[[i]] %>%
      ggplot(aes(x = X_values,
                 y = Y_values,
                 z = predicted_values/max(predicted_values))) +
      metR::geom_contour_fill(aes(fill = ..level..),
                              color = "black",
                              size = 0.03,
                              bins = 15) +
      xlab("log2 of SF (cpd)") +
      ylab("log2 of TF (Hz)") +
      geom_text(
        x = range(ER_initrns_speed_gauss_predicts[[i]]$X_values)[2] - 1,
        y = range(ER_initrns_speed_gauss_predicts[[i]]$Y_values)[1] + 0.5,
        label = paste0(
          "R2 = ", round(ER_initrns_speed_gauss_summary[i,]$R2_adj, 2)
        ),
        color = "white"
      ) +
      ggtitle(
        paste0(
          names(ER_initrns_speed_gauss_predicts)[i],
          " Q = 0")
      ) +
      coord_fixed(expand = FALSE) +
      theme_classic() +
      theme(plot.title = element_text(size = 8))
  }
}

#### export ####

# write.csv(
#   ER_initrns_unconst_gauss_summary,
#   "./analyses_and_figures/figure 4/supplementary/ER_initrns_unconst_gauss_summary_trunc.csv"
# )
# write.csv(
#   ER_initrns_indep_gauss_summary,
#   "./analyses_and_figures/figure 4/supplementary/ER_initrns_indep_gauss_summary_trunc.csv"
# )
# write.csv(
#   ER_initrns_speed_gauss_summary,
#   "./analyses_and_figures/figure 4/supplementary/ER_initrns_speed_gauss_summary_trunc.csv"
# )
# write.csv(
#   ER_initrns_Z_paramdf,
#   "./analyses_and_figures/figure 4/supplementary/ER_initrns_Z_statistics.csv"
# )


all_initrns_cows <- NULL
for (i in 1:length(ER_initrns_contour_dat)
) {
  print(i)
  all_initrns_cows[[i]] <-
    plot_grid(
      ER_initrns_allcells_contour_plots[[i]],
      ER_initrns_allcells_unconst_pred_plots[[i]],
      ER_initrns_allcells_indep_pred_plots[[i]],
      ER_initrns_allcells_speed_pred_plots[[i]],
      ncol = 2
    )
}

# pdf("./analyses_and_figures/figure 4/LM_ER_initial_transients.pdf",
#     width = 8.5, height = 11,
#     paper = "letter", bg = "white",
#     pagecentre = TRUE, colormodel = "srgb")
# for (i in 1:length(all_initrns_cows)) {
#   plot(all_initrns_cows[[i]])
# }
# dev.off()
# beepr::beep(4)
