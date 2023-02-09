##### Options ######
## ggplot margins
## t, r, b, l
plot_margs <- unit(c(0.1, 0.1, 0.1, 0.1), "cm")

## Breaks
breaks_common = seq(0.1, 1, by = 0.1)
bins_common = 15
bins_alt = 15

## Contour line thickness
contour_wds <- 0.04

#### Functions ####
## Note from VBB: for this function to work, peak location, X- and Y- variance,
## amplitude, and Q must all be pre-determined for a 2D gaussian. This function
## then produces the predicted height at each combination of X and Y 
predict_gaussian_xl <- function(X_values,
                                Y_values,
                                A,
                                X_peak,
                                X_var,
                                Q,
                                Y_peak,
                                Y_var){
  result <- NULL
  for (i in 1:length(X_values)){
    result[i] <- A * exp(-((X_values[i]-log2(X_peak))^2)/(X_var^2)) * exp(-(Y_values[i]-log2(2^((Q+1)*(X_values[i]-log2(X_peak))+log2(Y_peak))))^2/(Y_var^2))
  }
  return(result)
}

## Note from VBB: the function that follows extracts parameters from 
## an excel file in which the solver function has been run.
extract_gaussian_parameters <- function(file, method = "unconstrained"){
  test <- read_excel(file)
  as.data.frame(test) -> test
  
  if (method == "unconstrained"){
    X_values <- as.vector(na.omit(test[,3]))
    Y_values <- as.vector(na.omit(test[,4]))
    response <- as.vector(na.omit(test[,5])) ## prob has mean as final value
    if (length(response) > length(Y_values)){
      response <- response[1:(length(response) - 1)]
    }
    norm_g_resp <- as.vector(na.omit(test[,10]))
    g_resp <- as.vector(na.omit(test[,7]))
    if (length(g_resp) > length(Y_values)){
      g_resp <- g_resp[1:(length(g_resp) - 1)]
    }
    
    A <- as.numeric(test$A[1])
    X_peak <- as.numeric(test[1, 20])
    X_var <- as.numeric(test[1, 18])
    Q <- as.numeric(test[1, 17])
    Y_peak <- as.numeric(test[1, 21])
    Y_var <- as.numeric(test[1, 19])
  }
  
  if (method == "Sp"){
    X_values <- as.vector(na.omit(test[,3]))
    Y_values <- as.vector(na.omit(test[,4]))
    response <- as.vector(na.omit(test[,5])) ## prob has mean as final value
    if (length(response) > length(Y_values)){
      response <- response[1:(length(response) - 1)]
    }
    norm_g_resp <- as.vector(na.omit(test[,11]))
    g_resp <- as.vector(na.omit(test[,8]))
    if (length(g_resp) > length(Y_values)){
      g_resp <- g_resp[1:(length(g_resp) - 1)]
    }
    A <- as.numeric(test$A[1])
    X_peak <- as.numeric(test[3, 20])
    X_var <- as.numeric(test[3, 18])
    Q <- as.numeric(test[3, 17])
    Y_peak <- as.numeric(test[3, 21])
    Y_var <- as.numeric(test[3, 19])
  }
  
  if (method == "Ind"){
    X_values <- as.vector(na.omit(test[,3]))
    Y_values <- as.vector(na.omit(test[,4]))
    response <- as.vector(na.omit(test[,5])) ## prob has mean as final value
    if (length(response) > length(Y_values)){
      response <- response[1:(length(response) - 1)]
    }
    norm_g_resp <- as.vector(na.omit(test[,12]))
    g_resp <- as.vector(na.omit(test[,9]))
    if (length(g_resp) > length(Y_values)){
      g_resp <- g_resp[1:(length(g_resp) - 1)]
    }
    A <- as.numeric(test$A[1])
    X_peak <- as.numeric(test[5, 20])
    X_var <- as.numeric(test[5, 18])
    Q <- as.numeric(test[5, 17])
    Y_peak <- as.numeric(test[5, 21])
    Y_var <- as.numeric(test[5, 19])
  }
  
  obj <- data.frame(X_values = X_values,
                    Y_values = Y_values,
                    response = response,
                    norm_g_resp = norm_g_resp,
                    g_resp = g_resp,
                    A = A,
                    X_peak = X_peak,
                    X_var = X_var,
                    Q = Q,
                    Y_peak = Y_peak,
                    Y_var = Y_var)
  return(obj)
}
##### Panel 3ABC Zebra finch #####
file_3A2E_5304 <- "./data/spatemp/LM_ER_fullsweep_20bins/ALL_LM_SPATEMP/zebrafinch_LM_spatemp_min_err_gaussian/STC02-LM-02-5304-Cell 2 redone july 2020.xls"

## Columns needed
cols_needed <- 1:12

data_3A2E_5304 <- 
  read_excel(file_3A2E_5304)  %>%
  select(all_of(cols_needed)) %>%
  drop_na()

###### _3A2E_i - contour plot #####
p_3A2E_i <-
  ggplot(data_3A2E_5304, aes(`Log2 SF`, Log2TF, z = `G-Norm-Un`)) + 
  metR::geom_contour_fill(#aes(fill = ..level..), 
                          color = "black", size = 0.03,
                          bins = 14) +
  scale_fill_gradientn(
    colors = c(
      rgb(0, 0, 0, maxColorValue = 255), #black
      rgb(255, 128, 0, maxColorValue = 255), #orange
      rgb(255, 255, 255, maxColorValue = 255)
    )) +
  xlab("spatial frequency (cpd)") +
  ylab(expression(paste("temporal \nfrequency (Hz)"))) +
  expand_limits(x = c(-6, 0),
                y = c(-5, 4)) +
  scale_y_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     labels = c(two_neg_five, two_neg_three, two_neg_one, 
                                two_one, two_three, two_four)) +
  scale_x_continuous(breaks = seq(-6, 0, by = 1), 
                     labels = c(two_neg_six, two_neg_five, two_neg_four,
                                two_neg_three,
                                two_neg_two, two_neg_one, two_zero),
                     limits = c(-6, 0)) +
  # geom_abline(slope = 1, intercept = 11,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 9,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 7,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 5,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 3,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 1,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = -1,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = -3,  colour = "grey30", size = 0.2) +
  coord_cartesian(xlim = c(-6, 0), 
                  ylim = c(-5, 4),
                  expand = FALSE) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 6),
    axis.title.y = element_text(hjust = -0.15),
    axis.title.x = element_text(hjust = 0),
    axis.title = element_text(size = 7),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    legend.position="none"
  )

params_3A2E_5304_uncon <-
  extract_gaussian_parameters(file_3A2E_5304, method = "unconstrained")

params_3A2E_5304_sp <-
  extract_gaussian_parameters(file_3A2E_5304, method = "Sp")

params_3A2E_5304_ind <-
  extract_gaussian_parameters(file_3A2E_5304, method = "Ind")

grid_3A2E_5304 <-
  expand.grid(X_values = seq(from = -5, to = 0, by = 0.1),
              Y_values = seq(from = -5, to = 4, by = 0.1))


###### _3A2E_ii - unconstrained plot #####
preds_3A2E_5304_uncon <-
  predict_gaussian_xl(
    X_values = grid_3A2E_5304$X_values,
    Y_values = grid_3A2E_5304$Y_values,
    A = params_3A2E_5304_uncon$A[1],
    X_peak = params_3A2E_5304_uncon$X_peak[1],
    Y_peak = params_3A2E_5304_uncon$Y_peak[1],
    Q = params_3A2E_5304_uncon$Q[1],
    X_var = params_3A2E_5304_uncon$X_var[1],
    Y_var = params_3A2E_5304_uncon$Y_var[1]
  )

data_3A2E_5304_uncon <- 
  data.frame(SFFIT = grid_3A2E_5304$X_values,
             TFFIT = grid_3A2E_5304$Y_values,
             FIT_NORM = preds_3A2E_5304_uncon/max(preds_3A2E_5304_uncon))

p_3A2E_ii <-
  ggplot(data_3A2E_5304_uncon, aes(SFFIT, TFFIT, z = FIT_NORM)) + 
  metR::geom_contour_fill(#aes(fill = ..level..), 
                          size = contour_wds,
                          bins = bins_common) +
  scale_fill_gradientn(
    colors = c(
      rgb(0, 0, 0, maxColorValue = 255), #black
      rgb(255, 128, 0, maxColorValue = 255), #orange
      rgb(255, 255, 255, maxColorValue = 255)
    )) + xlab(" ") +
  ylab(" ") +
  expand_limits(x = c(-6, 0),
                y = c(-5, 4)) +
  scale_y_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     labels = c(two_neg_five, two_neg_three, two_neg_one, 
                                two_one, two_three, two_four)) +
  scale_x_continuous(breaks = seq(-6, 0, by = 1), 
                     labels = c(two_neg_six, two_neg_five, two_neg_four,
                                two_neg_three,
                                two_neg_two, two_neg_one, two_zero),
                     limits = c(-6, 0)) +
  # geom_abline(slope = 1, intercept = 11,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 9,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 7,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 5,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 3,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 1,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = -1,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = -3,  colour = "grey30", size = 0.2) +
  coord_cartesian(xlim = c(-6, 0), 
                  ylim = c(-5, 4),
                  expand = FALSE) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6), #element_blank(),
    axis.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    #panel.background = element_rect(fill = col_hb, colour = col_hb),
    legend.position="none"
  )

###### _3A2E_iii - Sp plot #####
preds_3A2E_5304_sp <-
  predict_gaussian_xl(
    X_values = grid_3A2E_5304$X_values,
    Y_values = grid_3A2E_5304$Y_values,
    A = params_3A2E_5304_sp$A[1],
    X_peak = params_3A2E_5304_sp$X_peak[1],
    Y_peak = params_3A2E_5304_sp$Y_peak[1],
    Q = params_3A2E_5304_sp$Q[1],
    X_var = params_3A2E_5304_sp$X_var[1],
    Y_var = params_3A2E_5304_sp$Y_var[1]
  )

data_3A2E_5304_sp <- 
  data.frame(SFFIT = grid_3A2E_5304$X_values,
             TFFIT = grid_3A2E_5304$Y_values,
             FIT_NORM = preds_3A2E_5304_sp/max(preds_3A2E_5304_sp))

p_3A2E_iii <-
  ggplot(data_3A2E_5304_sp, aes(SFFIT, TFFIT, z = FIT_NORM)) + 
  metR::geom_contour_fill(#aes(fill = ..level..), 
                          size = contour_wds,
                          bins = bins_common) +
  scale_fill_gradientn(
    colors = c(
      rgb(0, 0, 0, maxColorValue = 255), #black
      rgb(255, 128, 0, maxColorValue = 255), #orange
      rgb(255, 255, 255, maxColorValue = 255)
    )) + xlab(" ") +
  ylab(" ") +
  expand_limits(x = c(-6, 0),
                y = c(-5, 4)) +
  scale_y_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     labels = c(two_neg_five, two_neg_three, two_neg_one, 
                                two_one, two_three, two_four)) +
  scale_x_continuous(breaks = seq(-6, 0, by = 1), 
                     labels = c(two_neg_six, two_neg_five, two_neg_four,
                                two_neg_three,
                                two_neg_two, two_neg_one, two_zero),
                     limits = c(-6, 0)) +
  coord_cartesian(xlim = c(-6, 0), 
                  ylim = c(-5, 4),
                  expand = FALSE) +
  theme_classic() +
  # geom_abline(slope = 1, intercept = 11,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 9,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 7,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 5,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 3,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 1,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = -1,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = -3,  colour = "grey30", size = 0.2) +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6), #element_blank(),
    axis.title = element_text(size = 6),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    #panel.background = element_rect(fill = col_hb, colour = col_hb),
    legend.position="none"
  )

###### _3A2E_iv - Ind plot #####
preds_3A2E_5304_ind <-
  predict_gaussian_xl(
    X_values = grid_3A2E_5304$X_values,
    Y_values = grid_3A2E_5304$Y_values,
    A = params_3A2E_5304_ind$A[1],
    X_peak = params_3A2E_5304_ind$X_peak[1],
    Y_peak = params_3A2E_5304_ind$Y_peak[1],
    Q = params_3A2E_5304_ind$Q[1],
    X_var = params_3A2E_5304_ind$X_var[1],
    Y_var = params_3A2E_5304_ind$Y_var[1]
  )

data_3A2E_5304_ind <- 
  data.frame(SFFIT = grid_3A2E_5304$X_values,
             TFFIT = grid_3A2E_5304$Y_values,
             FIT_NORM = preds_3A2E_5304_ind/max(preds_3A2E_5304_ind))

p_3A2E_iv <-
  ggplot(data_3A2E_5304_ind, aes(SFFIT, TFFIT, z = FIT_NORM)) + 
  metR::geom_contour_fill(#aes(fill = ..level..), 
                          size = contour_wds,
                          bins = 16) +
  scale_fill_gradientn(
    colors = c(
      rgb(0, 0, 0, maxColorValue = 255), #black
      rgb(255, 128, 0, maxColorValue = 255), #orange
      rgb(255, 255, 255, maxColorValue = 255)
    )) +  xlab(" ") +
  ylab(" ") +
  expand_limits(x = c(-6, 0),
                y = c(-5, 4)) +
  scale_y_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     labels = c(two_neg_five, two_neg_three, two_neg_one, 
                                two_one, two_three, two_four)) +
  scale_x_continuous(breaks = seq(-6, 0, by = 1), 
                     labels = c(two_neg_six, two_neg_five, two_neg_four,
                                two_neg_three,
                                two_neg_two, two_neg_one, two_zero),
                     limits = c(-6, 0)) +
  # geom_abline(slope = 1, intercept = 11,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 9,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 7,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 5,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 3,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 1,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = -1,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = -3,  colour = "grey30", size = 0.2) +
  coord_cartesian(xlim = c(-6, 0), 
                  ylim = c(-5, 4),
                  expand = FALSE) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6), #element_blank(),
    axis.title = element_text(size = 6),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    #panel.background = element_rect(fill = col_hb, colour = col_hb),
    legend.position="none"
  )

##### Panel 3FGH Hummingbird #####
file_3F2J <- "./data/spatemp/LM_ER_fullsweep_20bins/ALL_LM_SPATEMP/hummingbird_LM_spatemp_min_err_gaussian/CALSTC09-LM-03-3800-Cell 2 redone july 2020.xls"

## Columns needed
cols_needed <- 1:12

data_3F2J <- 
  read_excel(file_3F2J)  %>%
  select(all_of(cols_needed)) %>%
  drop_na()

###### _3F2J_i - contour plot #####
p_3F2J_i <-
  ggplot(data_3F2J, aes(`Log2 SF`, Log2TF, z = `G-Norm-Un`)) + 
  metR::geom_contour_fill(#aes(fill = ..level..), 
    color = "black", size = 0.03,
    bins = 14) +
  scale_fill_gradientn(
    colors=c(rgb(0, 0, 0, maxColorValue = 255), #black
             rgb(255, 51, 153, maxColorValue = 255), #pink
             rgb(255, 255, 255, maxColorValue = 255))
  ) +
  xlab("spatial frequency (cpd)") +
  ylab(expression(paste("temporal \nfrequency (Hz)"))) +
  expand_limits(x = c(-6, 0),
                y = c(-5, 4)) +
  scale_y_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     labels = c(two_neg_five, two_neg_three, two_neg_one, 
                                two_one, two_three, two_four)) +
  scale_x_continuous(breaks = seq(-6, 0, by = 1), 
                     labels = c(two_neg_six, two_neg_five, two_neg_four,
                                two_neg_three,
                                two_neg_two, two_neg_one, two_zero),
                     limits = c(-6.05, 0.05)) +
  # geom_abline(slope = 1, intercept = 11,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 9,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 7,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 5,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 3,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 1,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = -1,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = -3,  colour = "grey30", size = 0.2) +
  coord_cartesian(xlim = c(-6, 0), 
                  ylim = c(-5, 4),
                  expand = FALSE) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 6),
    axis.title.y = element_text(hjust = -0.15),
    axis.title.x = element_text(hjust = 0),
    axis.title = element_text(size = 7),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    legend.position="none"
  )

params_3F2J_uncon <-
  extract_gaussian_parameters(file_3F2J, method = "unconstrained")

params_3F2J_sp <-
  extract_gaussian_parameters(file_3F2J, method = "Sp")

params_3F2J_ind <-
  extract_gaussian_parameters(file_3F2J, method = "Ind")

grid_3F2J <-
  expand.grid(X_values = seq(from = -6, to = 0, by = 0.1),
              Y_values = seq(from = -5, to = 4, by = 0.1))


###### _3F2J_ii - unconstrained plot #####
preds_3F2J_uncon <-
  predict_gaussian_xl(
    X_values = grid_3F2J$X_values,
    Y_values = grid_3F2J$Y_values,
    A = params_3F2J_uncon$A[1],
    X_peak = params_3F2J_uncon$X_peak[1],
    Y_peak = params_3F2J_uncon$Y_peak[1],
    Q = params_3F2J_uncon$Q[1],
    X_var = params_3F2J_uncon$X_var[1],
    Y_var = params_3F2J_uncon$Y_var[1]
  )

data_3F2J_uncon <- 
  data.frame(SFFIT = grid_3F2J$X_values,
             TFFIT = grid_3F2J$Y_values,
             FIT_NORM = preds_3F2J_uncon/max(preds_3F2J_uncon))

p_3F2J_ii <-
  ggplot(data_3F2J_uncon, aes(SFFIT, TFFIT, z = FIT_NORM)) + 
  metR::geom_contour_fill(#aes(fill = ..level..), 
    size = contour_wds,
    bins = bins_common) +
  scale_fill_gradientn(
    colors=c(rgb(0, 0, 0, maxColorValue = 255), #black
             rgb(255, 51, 153, maxColorValue = 255), #pink
             rgb(255, 255, 255, maxColorValue = 255))
  ) + xlab(" ") +
  ylab(" ") +
  expand_limits(x = c(-6, 0),
                y = c(-5, 4)) +
  scale_y_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     labels = c(two_neg_five, two_neg_three, two_neg_one, 
                                two_one, two_three, two_four)) +
  scale_x_continuous(breaks = seq(-6, 0, by = 1), 
                     labels = c(two_neg_six, two_neg_five, two_neg_four,
                                two_neg_three,
                                two_neg_two, two_neg_one, two_zero),
                     limits = c(-6.05, 0.05)) +
  # geom_abline(slope = 1, intercept = 11,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 9,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 7,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 5,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 3,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 1,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = -1,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = -3,  colour = "grey30", size = 0.2) +
  coord_cartesian(xlim = c(-6, 0), 
                  ylim = c(-5, 4),
                  expand = FALSE) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6), #element_blank(),
    axis.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    #panel.background = element_rect(fill = col_hb, colour = col_hb),
    legend.position="none"
  )

###### _3F2J_iii - Sp plot #####
preds_3F2J_sp <-
  predict_gaussian_xl(
    X_values = grid_3F2J$X_values,
    Y_values = grid_3F2J$Y_values,
    A = params_3F2J_sp$A[1],
    X_peak = params_3F2J_sp$X_peak[1],
    Y_peak = params_3F2J_sp$Y_peak[1],
    Q = params_3F2J_sp$Q[1],
    X_var = params_3F2J_sp$X_var[1],
    Y_var = params_3F2J_sp$Y_var[1]
  )

data_3F2J_sp <- 
  data.frame(SFFIT = grid_3F2J$X_values,
             TFFIT = grid_3F2J$Y_values,
             FIT_NORM = preds_3F2J_sp/max(preds_3F2J_sp))

p_3F2J_iii <-
  ggplot(data_3F2J_sp, aes(SFFIT, TFFIT, z = FIT_NORM)) + 
  metR::geom_contour_fill(#aes(fill = ..level..), 
    size = contour_wds,
    bins = bins_common) +
  scale_fill_gradientn(
    colors=c(rgb(0, 0, 0, maxColorValue = 255), #black
             rgb(255, 51, 153, maxColorValue = 255), #pink
             rgb(255, 255, 255, maxColorValue = 255))
  ) + xlab(" ") +
  ylab(" ") +
  expand_limits(x = c(-6, 0),
                y = c(-5, 4)) +
  scale_y_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     labels = c(two_neg_five, two_neg_three, two_neg_one, 
                                two_one, two_three, two_four)) +
  scale_x_continuous(breaks = seq(-6, 0, by = 1), 
                     labels = c(two_neg_six, two_neg_five, two_neg_four,
                                two_neg_three,
                                two_neg_two, two_neg_one, two_zero),
                     limits = c(-6.05, 0.05)) +
  coord_cartesian(xlim = c(-6, 0), 
                  ylim = c(-5, 4),
                  expand = FALSE) +
  theme_classic() +
  # geom_abline(slope = 1, intercept = 11,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 9,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 7,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 5,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 3,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 1,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = -1,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = -3,  colour = "grey30", size = 0.2) +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6), #element_blank(),
    axis.title = element_text(size = 6),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    #panel.background = element_rect(fill = col_hb, colour = col_hb),
    legend.position="none"
  )

###### _3F2J_iv - Ind plot #####
preds_3F2J_ind <-
  predict_gaussian_xl(
    X_values = grid_3F2J$X_values,
    Y_values = grid_3F2J$Y_values,
    A = params_3F2J_ind$A[1],
    X_peak = params_3F2J_ind$X_peak[1],
    Y_peak = params_3F2J_ind$Y_peak[1],
    Q = params_3F2J_ind$Q[1],
    X_var = params_3F2J_ind$X_var[1],
    Y_var = params_3F2J_ind$Y_var[1]
  )

data_3F2J_ind <- 
  data.frame(SFFIT = grid_3F2J$X_values,
             TFFIT = grid_3F2J$Y_values,
             FIT_NORM = preds_3F2J_ind/max(preds_3F2J_ind))

p_3F2J_iv <-
  ggplot(data_3F2J_ind, aes(SFFIT, TFFIT, z = FIT_NORM)) + 
  metR::geom_contour_fill(#aes(fill = ..level..), 
    size = contour_wds,
    bins = 16) +
  scale_fill_gradientn(
    colors=c(rgb(0, 0, 0, maxColorValue = 255), #black
             rgb(255, 51, 153, maxColorValue = 255), #pink
             rgb(255, 255, 255, maxColorValue = 255))
  ) +  xlab(" ") +
  ylab(" ") +
  expand_limits(x = c(-6, 0),
                y = c(-5, 4)) +
  scale_y_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     labels = c(two_neg_five, two_neg_three, two_neg_one, 
                                two_one, two_three, two_four)) +
  scale_x_continuous(breaks = seq(-6, 0, by = 1), 
                     labels = c(two_neg_six, two_neg_five, two_neg_four,
                                two_neg_three,
                                two_neg_two, two_neg_one, two_zero),
                     limits = c(-6.05, 0.05)) +
  # geom_abline(slope = 1, intercept = 11,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 9,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 7,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 5,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 3,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = 1,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = -1,  colour = "grey30", size = 0.2) +
  # geom_abline(slope = 1, intercept = -3,  colour = "grey30", size = 0.2) +
  coord_cartesian(xlim = c(-6, 0), 
                  ylim = c(-5, 4),
                  expand = FALSE) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6), #element_blank(),
    axis.title = element_text(size = 6),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    #panel.background = element_rect(fill = col_hb, colour = col_hb),
    legend.position="none"
  )

#### Panel 3 DE Zebra finch #####
###### _3D - TF vs resp #####
panel_3D <-
  data_3A2E_5304 %>%
  ggplot(aes(x = Log2TF, y = Response, group = as.factor(`Log2 SF`))) +
  geom_line(aes(col = as.factor(`Log2 SF`))) +
  # geom_point(aes(fill = as.factor(`Log2 SF`)), 
  #            col = "black", pch = 21, stroke = 0.3) +
  scale_color_viridis_d(option = "C",
                        direction = -1,
                        end = 0.8,
                        name = "SF (cpd)", 
                        labels = c(two_neg_six, two_neg_five, two_neg_four,
                                   two_neg_three,two_neg_two, two_neg_one,
                                   two_zero))+
  scale_fill_viridis_d(option = "C",
                       direction = -1,
                       #end = 0.8,
                       guide = 'none') +
  xlab("temporal frequency (Hz)") +
  ylab(expression(paste("mean response \n(spikes/s)"))) +
  expand_limits(x = c(-5, 4)) +
  scale_x_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     labels = c(two_neg_five, two_neg_three, two_neg_one, 
                                two_one, two_three, two_four)) +
  theme_classic() +
  guides(shape = guide_legend(override.aes = list(size = 0.5)),
         color = guide_legend(override.aes = list(size = 0.5))) +
  theme(
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    legend.title = element_text(size = 6), 
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.3, "lines"),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(0,0,0,0)
  )


###### _3E - Vel vs resp #####
panel_3E <-
  data_3A2E_5304 %>%
  mutate(speed_raw = TF/SF) %>%
  mutate(speed_log = log10(speed_raw)) %>%
  ggplot(aes(x = speed_log, y = Response, group = as.factor(`Log2 SF`))) +
  geom_line(aes(col = as.factor(`Log2 SF`))) +
  # geom_point(aes(fill = as.factor(`Log2 SF`)), 
  #            col = "black", pch = 21, stroke = 0.3) +
  scale_color_viridis_d(option = "C",
                        direction = -1,
                        end = 0.8,
                        name = "SF (cpd)", 
                        labels = c(two_neg_six, two_neg_five, two_neg_four,
                                   two_neg_three,two_neg_two, two_neg_one,
                                   two_zero))+
  scale_fill_viridis_d(option = "C",
                       direction = -1,
                       #end = 0.8,
                       guide = 'none') +
  xlab("velocity (°/s)") +
  ylab(" ") +
  expand_limits(x = c(-2, 3)) +
  scale_x_continuous(breaks = seq(from = -2, to = 3, by = 1),
                     labels = c(ten_neg_two, ten_neg_one, ten_zero,
                                ten_one, ten_two, ten_three)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 6),
    axis.title.y = element_blank(),
    axis.title = element_text(size = 7),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    legend.position = "none"
  )


#### Panel 3 IJ Hummingbird #####
###### _3I - TF vs resp #####
panel_3I <-
  data_3F2J %>%
  ggplot(aes(x = Log2TF, y = Response, group = as.factor(`Log2 SF`))) +
  geom_line(aes(col = as.factor(`Log2 SF`))) +
  # geom_point(aes(fill = as.factor(`Log2 SF`)), 
  #            col = "black", pch = 21, stroke = 0.3) +
  scale_color_viridis_d(option = "C",
                        direction = -1,
                        end = 0.8,
                        name = "SF (cpd)", 
                        labels = c(two_neg_six, two_neg_five, two_neg_four,
                                   two_neg_three,two_neg_two, two_neg_one,
                                   two_zero))+
  scale_fill_viridis_d(option = "C",
                       direction = -1,
                       #end = 0.8,
                       guide = 'none') +
  xlab("temporal frequency (Hz)") +
  ylab(expression(paste("mean response \n(spikes/s)"))) +
  expand_limits(x = c(-5, 4)) +
  scale_x_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     labels = c(two_neg_five, two_neg_three, two_neg_one, 
                                two_one, two_three, two_four)) +
  theme_classic() +
  guides(shape = guide_legend(override.aes = list(size = 0.5)),
         color = guide_legend(override.aes = list(size = 0.5))) +
  theme(
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    legend.title = element_text(size = 6), 
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.3, "lines"),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(0,0,0,0)
  )


###### _3J - Vel vs resp #####
panel_3J <-
  data_3F2J %>%
  mutate(speed_raw = TF/SF) %>%
  mutate(speed_log = log10(speed_raw)) %>%
  ggplot(aes(x = speed_log, y = Response, group = as.factor(`Log2 SF`))) +
  geom_line(aes(col = as.factor(`Log2 SF`))) +
  # geom_point(aes(fill = as.factor(`Log2 SF`)), 
  #            col = "black", pch = 21, stroke = 0.3) +
  scale_color_viridis_d(option = "C",
                        direction = -1,
                        end = 0.8,
                        name = "SF (cpd)", 
                        labels = c(two_neg_six, two_neg_five, two_neg_four,
                                   two_neg_three,two_neg_two, two_neg_one,
                                   two_zero))+
  scale_fill_viridis_d(option = "C",
                       direction = -1,
                       #end = 0.8,
                       guide = 'none') +
  xlab("velocity (°/s)") +
  ylab(" ") +
  expand_limits(x = c(-2, 3)) +
  scale_x_continuous(breaks = seq(from = -2, to = 3, by = 1),
                     labels = c(ten_neg_two, ten_neg_one, ten_zero,
                                ten_one, ten_two, ten_three)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 6),
    axis.title.y = element_blank(),
    axis.title = element_text(size = 7),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    legend.position = "none"
  )


