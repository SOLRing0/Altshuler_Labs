####NOTE#####
## All SF and TF data in these csvs have been log2 transformed

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


## Functions
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
      g_resp <- g_resp[1:(length(g_resp) - 1)] ## prob has mean as final value
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
      g_resp <- g_resp[1:(length(g_resp) - 1)] ## prob has mean as final value
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
      g_resp <- g_resp[1:(length(g_resp) - 1)] ## prob has mean as final value
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


######## Panel 4H Hummingbird ######
## Subpanels will be labeled i through iv

#### _hb_top_left (subpanel i) ####
# data_1J_3226 <- 
#   read_csv("./Figure 4/panel_h_to_j/hummingbird/CALSTC01-LM-05-3226-Cell 1.csv") %>%
#   mutate(FIT_NORM = FIT/max(FIT))
file_1J_3226 <- "./data/spatemp/LM_ER_fullsweep_20bins/ALL_LM_SPATEMP/hummingbird_LM_spatemp_min_err_gaussian/CALSTC01-LM-05-3226-Cell 1 redone july 2020.xls"

params_1J_3226 <-
  extract_gaussian_parameters(file_1J_3226, method = "unconstrained")

grid_1J_3226 <-
  expand.grid(X_values = seq(from = -6, to = 0, by = 0.1),
              Y_values = seq(from = -5, to = 4, by = 0.1))
preds_1J_3226 <-
  predict_gaussian_xl(
    X_values = grid_1J_3226$X_values,
    Y_values = grid_1J_3226$Y_values,
    A = params_1J_3226$A[1],
    X_peak = params_1J_3226$X_peak[1],
    Y_peak = params_1J_3226$Y_peak[1],
    Q = params_1J_3226$Q[1],
    X_var = params_1J_3226$X_var[1],
    Y_var = params_1J_3226$Y_var[1]
  )

data_1J_3226 <- data.frame(SFFIT = grid_1J_3226$X_values,
                           TFFIT = grid_1J_3226$Y_values,
                           FIT_NORM = preds_1J_3226/max(preds_1J_3226))
  
p_1J_i <-
  ggplot(data_1J_3226, aes(SFFIT, TFFIT, z = FIT_NORM)) + 
  metR::geom_contour_fill(#aes(fill = ..level..), 
                          bins = bins_common) +
  scale_fill_gradientn(
    colors=c(rgb(0, 0, 0, maxColorValue = 255), #black
             rgb(255, 51, 153, maxColorValue = 255), #pink
             rgb(255, 255, 255, maxColorValue = 255))
  ) +
  xlab(" ") +
  ylab(" ") +
  expand_limits(x = c(-6, 1),
                y = c(-5, 4.58)) +
  scale_y_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     labels = c("", two_neg_three, "", 
                                two_one, "", two_four)) +
  scale_x_continuous(breaks = seq(-6, 1, by = 1), 
                     limits = c(-6, 1)) +
  coord_cartesian(xlim = c(-6, 1),
                  ylim = c(-5, 4.58),
                  expand = FALSE) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    #panel.background = element_rect(fill = col_hb, colour = col_hb),
    legend.position="none"
  )

## Replace the lowest level of the gradient with the bird color
## See here for details of steps:
## https://stackoverflow.com/questions/41940000/modifying-ggplot-objects-after-creation
# q_1J_i <- ggplot_build(p_1J_i)
# q_1J_i_colors <- q_1J_i$data[[1]]$fill 
# q_1J_i$data[[1]]$fill <- 
#   replace(q_1J_i_colors,
#           q_1J_i_colors == unique(q_1J_i$data[[1]]$fill)[1],
#           col_hb) 
# 
# panel_1J_i <- ggplot_gtable(q_1J_i)

#### _hb_top_right (subpanel ii) ####
# data_1J_3588 <- 
#   read_csv("./Figure 4/panel_h_to_j/hummingbird/CALSTC07-LM-07-3588-Cell 1.csv") %>%
#   mutate(FIT_NORM = FIT/max(FIT))
file_1J_3588 <- "./data/spatemp/LM_ER_fullsweep_20bins/ALL_LM_SPATEMP/hummingbird_LM_spatemp_min_err_gaussian/CALSTC07-LM-07-3588-Cell 1 redone july 2020.xls"

params_1J_3588 <-
  extract_gaussian_parameters(file_1J_3588, method = "unconstrained")

grid_1J_3588 <-
  expand.grid(X_values = seq(from = -6, to = 0, by = 0.1),
              Y_values = seq(from = -5, to = 4, by = 0.1))
preds_1J_3588 <-
  predict_gaussian_xl(
    X_values = grid_1J_3588$X_values,
    Y_values = grid_1J_3588$Y_values,
    A = params_1J_3588$A[1],
    X_peak = params_1J_3588$X_peak[1],
    Y_peak = params_1J_3588$Y_peak[1],
    Q = params_1J_3588$Q[1],
    X_var = params_1J_3588$X_var[1],
    Y_var = params_1J_3588$Y_var[1]
  )

data_1J_3588 <- data.frame(SFFIT = grid_1J_3588$X_values,
                           TFFIT = grid_1J_3588$Y_values,
                           FIT_NORM = preds_1J_3588/max(preds_1J_3588))

p_1J_ii <-
  ggplot(data_1J_3588, aes(SFFIT, TFFIT, z = FIT_NORM)) + 
  metR::geom_contour_fill(#aes(fill = ..level..), 
                          bins = bins_common) +
  scale_fill_gradientn(
    colors=c(rgb(0, 0, 0, maxColorValue = 255), #black
             rgb(255, 51, 153, maxColorValue = 255), #pink
             rgb(255, 255, 255, maxColorValue = 255))
  ) +
  xlab(" ") +
  ylab(" ") +
  expand_limits(x = c(-6, 1),
                y = c(-5, 4.58)) +
  scale_y_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     labels = c(two_neg_five, "" , two_neg_one, 
                                "", two_three, "")) +
  scale_x_continuous(breaks = seq(-6, 1, by = 1), 
                     limits = c(-6, 1)) +
  coord_cartesian(xlim = c(-6, 1),
                  ylim = c(-5, 4.58),
                  expand = FALSE) +
  theme_classic() +
  theme(#axis.text = element_text(size = 6, color = "white"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    legend.position="none"
  )

#### _hb_bottom_left (subpanel iii) ####
# data_1J_3070 <- 
#   read_csv("./Figure 4/panel_h_to_j/hummingbird/CALSTC08-LM-03-3070-Cell 1.csv") %>%
#   mutate(FIT_NORM = FIT/max(FIT))
file_1J_3070 <- "./data/spatemp/LM_ER_fullsweep_20bins/ALL_LM_SPATEMP/hummingbird_LM_spatemp_min_err_gaussian/CALSTC08-LM-03-3070-Cell 1 redone july 2020.xls"

params_1J_3070 <-
  extract_gaussian_parameters(file_1J_3070, method = "unconstrained")

grid_1J_3070 <-
  expand.grid(X_values = seq(from = -6, to = 0, by = 0.1),
              Y_values = seq(from = -5, to = 4, by = 0.1))
preds_1J_3070 <-
  predict_gaussian_xl(
    X_values = grid_1J_3070$X_values,
    Y_values = grid_1J_3070$Y_values,
    A = params_1J_3070$A[1],
    X_peak = params_1J_3070$X_peak[1],
    Y_peak = params_1J_3070$Y_peak[1],
    Q = params_1J_3070$Q[1],
    X_var = params_1J_3070$X_var[1],
    Y_var = params_1J_3070$Y_var[1]
  )

data_1J_3070 <- data.frame(SFFIT = grid_1J_3070$X_values,
                           TFFIT = grid_1J_3070$Y_values,
                           FIT_NORM = preds_1J_3070/max(preds_1J_3070))
  
p_1J_iii <-
  ggplot(data_1J_3070, aes(SFFIT, TFFIT, z = FIT_NORM)) + 
  metR::geom_contour_fill(#aes(fill = ..level..), 
                          bins = bins_common+2) +
  scale_fill_gradientn(
    colors=c(rgb(0, 0, 0, maxColorValue = 255), #black
             rgb(255, 51, 153, maxColorValue = 255), #pink
             rgb(255, 255, 255, maxColorValue = 255))
  ) +
  xlab(" ") +
  ylab(" ") +
  expand_limits(x = c(-6, 1),
                y = c(-5, 4.58)) +
  scale_y_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     labels = c(two_neg_five, "" , two_neg_one, 
                                "", two_three, "")) +
  scale_x_continuous(breaks = seq(-6, 1, by = 1), 
                     labels = c(two_neg_six, "", two_neg_four, "",
                                two_neg_two, "", two_zero, ""), 
                     limits = c(-6, 0)) +
  coord_cartesian(xlim = c(-6, 1),
                  ylim = c(-5, 4.58),
                  expand = FALSE) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6),
    axis.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    #panel.background = element_rect(fill = col_hb, colour = col_hb),
    legend.position="none"
  )

#### _hb_bottom_right (subpanel iv) ####
# data_1J_3551 <- 
#   read_csv("./Figure 4/panel_h_to_j/hummingbird/CALSTC08-LM-04-3551-Cell 1.csv") %>%
#   mutate(FIT_NORM = FIT/max(FIT))
file_1J_3551 <- "./data/spatemp/LM_ER_fullsweep_20bins/ALL_LM_SPATEMP/hummingbird_LM_spatemp_min_err_gaussian/CALSTC08-LM-04-3551-Cell 1 redone july 2020.xls"

params_1J_3551 <-
  extract_gaussian_parameters(file_1J_3551, method = "unconstrained")

grid_1J_3551 <-
  expand.grid(X_values = seq(from = -6, to = 0, by = 0.1),
              Y_values = seq(from = -5, to = 4, by = 0.1))
preds_1J_3551 <-
  predict_gaussian_xl(
    X_values = grid_1J_3551$X_values,
    Y_values = grid_1J_3551$Y_values,
    A = params_1J_3551$A[1],
    X_peak = params_1J_3551$X_peak[1],
    Y_peak = params_1J_3551$Y_peak[1],
    Q = params_1J_3551$Q[1],
    X_var = params_1J_3551$X_var[1],
    Y_var = params_1J_3551$Y_var[1]
  )

data_1J_3551 <- data.frame(SFFIT = grid_1J_3551$X_values,
                           TFFIT = grid_1J_3551$Y_values,
                           FIT_NORM = preds_1J_3551/max(preds_1J_3551))

p_1J_iv <-
  ggplot(data_1J_3551, aes(SFFIT, TFFIT, z = FIT_NORM)) + 
  metR::geom_contour_fill(#aes(fill = ..level..), 
                          bins = bins_common) +
  scale_fill_gradientn(
    colors=c(rgb(0, 0, 0, maxColorValue = 255), #black
             rgb(255, 51, 153, maxColorValue = 255), #pink
             rgb(255, 255, 255, maxColorValue = 255))
  ) +
  xlab(" ") +
  ylab(" ") +
  expand_limits(x = c(-6, 1),
                y = c(-5, 4.58)) +
  scale_y_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     labels = c(two_neg_five, "" , two_neg_one, 
                                "", two_three, "")) +
  scale_x_continuous(breaks = seq(-6, 1, by = 1), 
                     labels = c("", two_neg_five, "", two_neg_three,
                                "", two_neg_one, "", two_one), 
                     limits = c(-6, 0)) +
  coord_cartesian(xlim = c(-6, 1),
                  ylim = c(-5, 4.58),
                  expand = FALSE) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 6),
    axis.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    legend.position="none"
  )

######## Panel 4I Zebra finch ######
## Subpanels will be labeled i through iv

#### _zf_top_left (subpanel i) ####
# data_1K_5249 <- 
#   read_csv("./Figure 4/panel_h_to_j/zebra_finch/STC03-LM-03-5249-Cell 1.csv") %>%
#   mutate(FIT_NORM = FIT/max(FIT))
file_1K_5249 <- "./data/spatemp/LM_ER_fullsweep_20bins/ALL_LM_SPATEMP/zebrafinch_LM_spatemp_min_err_gaussian/STC03-LM-03-5249-Cell 1 redone july 2020.xls"

params_1K_5249 <-
  extract_gaussian_parameters(file_1K_5249, method = "unconstrained")

grid_1K_5249 <-
  expand.grid(X_values = seq(from = -5, to = 0, by = 0.1),
              Y_values = seq(from = -5, to = 4, by = 0.1))
preds_1K_5249 <-
  predict_gaussian_xl(
    X_values = grid_1K_5249$X_values,
    Y_values = grid_1K_5249$Y_values,
    A = params_1K_5249$A[1],
    X_peak = params_1K_5249$X_peak[1],
    Y_peak = params_1K_5249$Y_peak[1],
    Q = params_1K_5249$Q[1],
    X_var = params_1K_5249$X_var[1],
    Y_var = params_1K_5249$Y_var[1]
  )

data_1K_5249 <- data.frame(SFFIT = grid_1K_5249$X_values,
                           TFFIT = grid_1K_5249$Y_values,
                           FIT_NORM = preds_1K_5249/max(preds_1K_5249))

p_1K_i <-
  ggplot(data_1K_5249, aes(SFFIT, TFFIT, z = FIT_NORM)) + 
  metR::geom_contour_fill(#aes(fill = ..level..), 
                          bins = bins_common) +
  scale_fill_gradientn(
    colors=c(rgb(0, 0, 0, maxColorValue = 255), #black
             rgb(255, 128, 0, maxColorValue = 255), #orange
             rgb(255, 255, 255, maxColorValue = 255))
  ) +
  xlab(" ") +
  ylab(" ") +
  expand_limits(x = c(-6, 1),
                y = c(-5, 4.58)) +
  scale_y_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     labels = c("", two_neg_three, "", 
                                two_one, "", two_four)) +
  scale_x_continuous(breaks = seq(-6, 1, by = 1), 
                     limits = c(-6, 1)) +
  coord_cartesian(xlim = c(-6, 1),
                  ylim = c(-5, 4.58),
                  expand = FALSE) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    #panel.background = element_rect(fill = col_hb, colour = col_hb),
    legend.position="none"
  )

#### _zf_top_right (subpanel ii) ####
# data_1K_5483 <- 
#   read_csv("./Figure 4/panel_h_to_j/zebra_finch/STC03-LM-03-5483-Cell 1.csv") %>%
#   mutate(FIT_NORM = FIT/max(FIT))
file_1K_5483 <- "./data/spatemp/LM_ER_fullsweep_20bins/ALL_LM_SPATEMP/zebrafinch_LM_spatemp_min_err_gaussian/STC03-LM-03-5483-Cell 1 redone july 2020.xls"

params_1K_5483 <-
  extract_gaussian_parameters(file_1K_5483, method = "unconstrained")

grid_1K_5483 <-
  expand.grid(X_values = seq(from = -5, to = 0, by = 0.1),
              Y_values = seq(from = -5, to = 4, by = 0.1))
preds_1K_5483 <-
  predict_gaussian_xl(
    X_values = grid_1K_5483$X_values,
    Y_values = grid_1K_5483$Y_values,
    A = params_1K_5483$A[1],
    X_peak = params_1K_5483$X_peak[1],
    Y_peak = params_1K_5483$Y_peak[1],
    Q = params_1K_5483$Q[1],
    X_var = params_1K_5483$X_var[1],
    Y_var = params_1K_5483$Y_var[1]
  )

data_1K_5483 <- data.frame(SFFIT = grid_1K_5483$X_values,
                           TFFIT = grid_1K_5483$Y_values,
                           FIT_NORM = preds_1K_5483/max(preds_1K_5483))

p_1K_ii <-
  ggplot(data_1K_5483, aes(SFFIT, TFFIT, z = FIT_NORM)) + 
  metR::geom_contour_fill(#aes(fill = ..level..), 
                          bins = bins_common) +
  scale_fill_gradientn(
    colors=c(rgb(0, 0, 0, maxColorValue = 255), #black
             rgb(255, 128, 0, maxColorValue = 255), #orange
             rgb(255, 255, 255, maxColorValue = 255))
  ) +
  xlab(" ") +
  ylab(" ") +
  expand_limits(x = c(-6, 1),
                y = c(-5, 4.58)) +
  scale_y_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     labels = c("", two_neg_three , "", 
                                two_one, "", two_four)) +
  scale_x_continuous(breaks = seq(-6, 1, by = 1), 
                     limits = c(-6, 1)) +
  coord_cartesian(xlim = c(-6, 1),
                  ylim = c(-5, 4.58),
                  expand = FALSE) +
  theme_classic() +
  theme(#axis.text = element_text(size = 6, color = "white"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    legend.position="none"
  )

#### _zf_bottom_left (subpanel iii) ####
# data_1K_5756 <- 
#   read_csv("./Figure 4/panel_h_to_j/zebra_finch/STC07-LM-05-5756-Cell 1.csv") %>%
#   mutate(FIT_NORM = FIT/max(FIT))
file_1K_5756 <- "./data/spatemp/LM_ER_fullsweep_20bins/ALL_LM_SPATEMP/zebrafinch_LM_spatemp_min_err_gaussian/STC07-LM-05-5756-Cell 1 redone july 2020.xls"

params_1K_5756 <-
  extract_gaussian_parameters(file_1K_5756, method = "unconstrained")

grid_1K_5756 <-
  expand.grid(X_values = seq(from = -5, to = 0, by = 0.1),
              Y_values = seq(from = -5, to = 4, by = 0.1))
preds_1K_5756 <-
  predict_gaussian_xl(
    X_values = grid_1K_5756$X_values,
    Y_values = grid_1K_5756$Y_values,
    A = params_1K_5756$A[1],
    X_peak = params_1K_5756$X_peak[1],
    Y_peak = params_1K_5756$Y_peak[1],
    Q = params_1K_5756$Q[1],
    X_var = params_1K_5756$X_var[1],
    Y_var = params_1K_5756$Y_var[1]
  )

data_1K_5756 <- data.frame(SFFIT = grid_1K_5756$X_values,
                           TFFIT = grid_1K_5756$Y_values,
                           FIT_NORM = preds_1K_5756/max(preds_1K_5756))

p_1K_iii <-
  ggplot(data_1K_5756, aes(SFFIT, TFFIT, z = FIT_NORM)) + 
  metR::geom_contour_fill(#aes(fill = ..level..), 
                          bins = bins_common) +
  scale_fill_gradientn(
    colors=c(rgb(0, 0, 0, maxColorValue = 255), #black
             rgb(255, 128, 0, maxColorValue = 255), #orange
             rgb(255, 255, 255, maxColorValue = 255))
  ) +
  xlab(" ") +
  ylab(" ") +
  expand_limits(x = c(-6, 1),
                y = c(-5, 4.58)) +
  scale_y_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     labels = c(two_neg_five, "" , two_neg_one, 
                                "", two_three, "")) +
  scale_x_continuous(breaks = seq(-6, 1, by = 1), 
                     labels = c(two_neg_six, "", two_neg_four, "",
                                two_neg_two, "", two_zero, ""), 
                     limits = c(-6, 0)) +
  coord_cartesian(xlim = c(-6, 1),
                  ylim = c(-5, 4.58),
                  expand = FALSE) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6),
    axis.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    #panel.background = element_rect(fill = col_hb, colour = col_hb),
    legend.position="none"
  )


#### _zf_bottom_right (subpanel iv) ####
# data_1K_5655 <-
#   read_csv(
#     "./Figure 4/panel_h_to_j/zebra_finch/STC08-LM-05-5655-Cell 1 again.csv") %>%
#   mutate(FIT_NORM = FIT/max(FIT))

## right cell
file_1K_5655r <- "./data/spatemp/LM_ER_fullsweep_20bins/ALL_LM_SPATEMP/zebrafinch_LM_spatemp_min_err_gaussian/STC08-LM-05-5655-Cell 1-Peak 1 redone july 2020.xls"

params_1K_5655r <-
  extract_gaussian_parameters(file_1K_5655r, method = "unconstrained")

grid_1K_5655r <-
  expand.grid(X_values = seq(from = -2.5, to = 0, by = 0.1),
              Y_values = seq(from = -5, to = 4, by = 0.1))
preds_1K_5655r <-
  predict_gaussian_xl(
    X_values = grid_1K_5655r$X_values,
    Y_values = grid_1K_5655r$Y_values,
    A = params_1K_5655r$A[1],
    X_peak = params_1K_5655r$X_peak[1],
    Y_peak = params_1K_5655r$Y_peak[1],
    Q = params_1K_5655r$Q[1],
    X_var = params_1K_5655r$X_var[1],
    Y_var = params_1K_5655r$Y_var[1]
  )

## left cell
file_1K_5655l <- "./data/spatemp/LM_ER_fullsweep_20bins/ALL_LM_SPATEMP/zebrafinch_LM_spatemp_min_err_gaussian/STC08-LM-05-5655-Cell 1-Peak 2 redone july 2020.xls"

params_1K_5655l <-
  extract_gaussian_parameters(file_1K_5655l, method = "unconstrained")

grid_1K_5655l <-
  expand.grid(X_values = seq(from = -5, to = -4, by = 0.1),
              Y_values = seq(from = -5, to = 4, by = 0.1))
preds_1K_5655l <-
  predict_gaussian_xl(
    X_values = grid_1K_5655l$X_values,
    Y_values = grid_1K_5655l$Y_values,
    A = params_1K_5655l$A[1],
    X_peak = params_1K_5655l$X_peak[1],
    Y_peak = params_1K_5655l$Y_peak[1],
    Q = params_1K_5655l$Q[1],
    X_var = params_1K_5655l$X_var[1],
    Y_var = params_1K_5655l$Y_var[1]
  )

amprat_1K_5655 <- params_1K_5655l$A[1] / params_1K_5655r$A[1]

data_1K_5655r_r <- data.frame(
  SFFIT = grid_1K_5655r$X_values,
  TFFIT = grid_1K_5655r$Y_values,
  FIT_NORM = preds_1K_5655r / max(preds_1K_5655r)
)

data_1K_5655l_l <- data.frame(
  SFFIT = grid_1K_5655l$X_values,
  TFFIT = grid_1K_5655l$Y_values,
  FIT_NORM = amprat_1K_5655 * preds_1K_5655l /
    max(preds_1K_5655l)
)

data_1K_5655 <- bind_rows(data_1K_5655r_r, data_1K_5655l_l)

p_1K_iv <-
  ggplot(data_1K_5655, aes(SFFIT, TFFIT, z = FIT_NORM)) + 
  metR::geom_contour_fill(#aes(fill = ..level..), 
                          bins = bins_common) +
  scale_fill_gradientn(
    colors=c(rgb(0, 0, 0, maxColorValue = 255), #black
             rgb(255, 128, 0, maxColorValue = 255), #orange
             rgb(255, 255, 255, maxColorValue = 255))
  ) +
  xlab(" ") +
  ylab(" ") +
  expand_limits(x = c(-6, 1),
                y = c(-5, 4.58)) +
  scale_y_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     labels = c(two_neg_five, "" , two_neg_one, 
                                "", two_three, "")) +
  scale_x_continuous(breaks = seq(-6, 1, by = 1), 
                     labels = c("", two_neg_five, "", two_neg_three,
                                "", two_neg_one, "", two_one), 
                     limits = c(-6, 0)) +
  coord_cartesian(xlim = c(-6, 1),
                  ylim = c(-5, 4.58),
                  expand = FALSE) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 6),
    axis.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    legend.position="none"
  )

###### Panel 4J Pigeon ####

#### _pg_top_left (subpanel i) ####
# data_1L_nc17LM2 <- 
#   read_csv("./Figure 4/panel_h_to_j/pigeon/nc17LM2.csv") %>%
#   mutate(FIT_NORM = FIT/max(FIT))
file_1L_nc17LM2 <- "./data/spatemp/LM_ER_fullsweep_20bins/ALL_LM_SPATEMP/pigeon_LM_spatemp_min_err_gaussian/nc17LM2 redo june 2018.xls"

params_1L_nc17LM2 <-
  extract_gaussian_parameters(file_1L_nc17LM2, method = "unconstrained")

grid_1L_nc17LM2 <-
  expand.grid(X_values = seq(from = -5, to = 0, by = 0.1),
              Y_values = seq(from = -5, to = 4, by = 0.1))
preds_1L_nc17LM2 <-
  predict_gaussian_xl(
    X_values = grid_1L_nc17LM2$X_values,
    Y_values = grid_1L_nc17LM2$Y_values,
    A = params_1L_nc17LM2$A[1],
    X_peak = params_1L_nc17LM2$X_peak[1],
    Y_peak = params_1L_nc17LM2$Y_peak[1],
    Q = params_1L_nc17LM2$Q[1],
    X_var = params_1L_nc17LM2$X_var[1],
    Y_var = params_1L_nc17LM2$Y_var[1]
  )

data_1L_nc17LM2 <- data.frame(SFFIT = grid_1L_nc17LM2$X_values,
                              TFFIT = grid_1L_nc17LM2$Y_values,
                              FIT_NORM = preds_1L_nc17LM2/max(preds_1L_nc17LM2))

p_1L_i <-
  ggplot(data_1L_nc17LM2, aes(SFFIT, TFFIT, z = FIT_NORM)) +
  metR::geom_contour_fill(#aes(fill = ..level..), 
                          bins = bins_common) +
  scale_fill_gradientn(
    colors=c(rgb(0, 0, 0, maxColorValue = 255), #black
             rgb(51, 153, 255, maxColorValue = 255), #blue
             rgb(255, 255, 255, maxColorValue = 255))
  ) +
  xlab(" ") +
  #ylab(" ") +
  expand_limits(x = c(-6, 1),
                y = c(-5, 4.58)) +
  scale_y_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     labels = c("", two_neg_three, "", 
                                two_one, "", two_four)) +
  scale_x_continuous(breaks = seq(-6, 1, by = 1), 
                     limits = c(-6, 1)) +
  coord_cartesian(xlim = c(-6, 1),
                  ylim = c(-5, 4.58),
                  expand = FALSE) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    #panel.background = element_rect(fill = col_hb, colour = col_hb),
    legend.position="none"
  )

#### _pg_top_right (subpanel ii) ####
# data_1L_nc24LM2 <- 
#   read_csv("./Figure 4/panel_h_to_j/pigeon/nc24LM2.csv") %>%
#   mutate(FIT_NORM = FIT/max(FIT))
file_1L_nc24LM2 <- "./data/spatemp/LM_ER_fullsweep_20bins/ALL_LM_SPATEMP/pigeon_LM_spatemp_min_err_gaussian/nc24LM2 redo june 2018.xls"

params_1L_nc24LM2 <-
  extract_gaussian_parameters(file_1L_nc24LM2, method = "unconstrained")

grid_1L_nc24LM2 <-
  expand.grid(X_values = seq(from = -5, to = 0, by = 0.1),
              Y_values = seq(from = -5, to = 4, by = 0.1))
preds_1L_nc24LM2 <-
  predict_gaussian_xl(
    X_values = grid_1L_nc24LM2$X_values,
    Y_values = grid_1L_nc24LM2$Y_values,
    A = params_1L_nc24LM2$A[1],
    X_peak = params_1L_nc24LM2$X_peak[1],
    Y_peak = params_1L_nc24LM2$Y_peak[1],
    Q = params_1L_nc24LM2$Q[1],
    X_var = params_1L_nc24LM2$X_var[1],
    Y_var = params_1L_nc24LM2$Y_var[1]
  )

data_1L_nc24LM2 <- data.frame(SFFIT = grid_1L_nc24LM2$X_values,
                              TFFIT = grid_1L_nc24LM2$Y_values,
                              FIT_NORM = preds_1L_nc24LM2/max(preds_1L_nc24LM2))

p_1L_ii <-
  ggplot(data_1L_nc24LM2, aes(SFFIT, TFFIT, z = FIT_NORM)) +
  metR::geom_contour_fill(#aes(fill = ..level..), 
                          bins = bins_common) +
  scale_fill_gradientn(
    colors=c(rgb(0, 0, 0, maxColorValue = 255), #black
             rgb(51, 153, 255, maxColorValue = 255), #blue
             rgb(255, 255, 255, maxColorValue = 255))
  ) +
  xlab(" ") +
  #ylab(" ") +
  expand_limits(x = c(-6, 1),
                y = c(-5, 4.58)) +
  scale_y_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     labels = c("", two_neg_three , "", 
                                two_one, "", two_four)) +
  scale_x_continuous(breaks = seq(-6, 1, by = 1), 
                     limits = c(-6, 1)) +
  coord_cartesian(xlim = c(-6, 1),
                  ylim = c(-5, 4.58),
                  expand = FALSE) +
  theme_classic() +
  theme(#axis.text = element_text(size = 6, color = "white"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    #panel.background = element_rect(fill = col_pg, colour = col_pg),
    legend.position="none"
  )

#### _pg_bottom_left (subpanel iii) ####
# data_1L_nc8wu2 <- 
#   read_csv("./Figure 4/panel_h_to_j/pigeon/nc8wu2.csv") %>%
#   mutate(FIT_NORM = FIT/max(FIT))
file_1L_nc8wu2 <- "./data/spatemp/LM_ER_fullsweep_20bins/ALL_LM_SPATEMP/pigeon_LM_spatemp_min_err_gaussian/nc8wu2 redo june 2018.xls"

params_1L_nc8wu2 <-
  extract_gaussian_parameters(file_1L_nc8wu2, method = "unconstrained")

grid_1L_nc8wu2 <-
  expand.grid(X_values = seq(from = -5, to = 1, by = 0.1),
              Y_values = seq(from = -5, to = 4, by = 0.1))
preds_1L_nc8wu2 <-
  predict_gaussian_xl(
    X_values = grid_1L_nc8wu2$X_values,
    Y_values = grid_1L_nc8wu2$Y_values,
    A = params_1L_nc8wu2$A[1],
    X_peak = params_1L_nc8wu2$X_peak[1],
    Y_peak = params_1L_nc8wu2$Y_peak[1],
    Q = params_1L_nc8wu2$Q[1],
    X_var = params_1L_nc8wu2$X_var[1],
    Y_var = params_1L_nc8wu2$Y_var[1]
  )

data_1L_nc8wu2 <- data.frame(SFFIT = grid_1L_nc8wu2$X_values,
                             TFFIT = grid_1L_nc8wu2$Y_values,
                             FIT_NORM = preds_1L_nc8wu2/max(preds_1L_nc8wu2))

p_1L_iii <-
  ggplot(data_1L_nc8wu2, aes(SFFIT, TFFIT, z = FIT_NORM)) +
  metR::geom_contour_fill(#aes(fill = ..level..), 
                          bins = bins_common) +
  scale_fill_gradientn(
    colors=c(rgb(0, 0, 0, maxColorValue = 255), #black
             rgb(51, 153, 255, maxColorValue = 255), #blue
             rgb(255, 255, 255, maxColorValue = 255))
  ) +
  xlab(" ") +
  #ylab(" ") +
  expand_limits(x = c(-6, 1),
                y = c(-5, 4.58)) +
  scale_y_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     labels = c(two_neg_five, "" , two_neg_one, 
                                "", two_three, "")) +
  scale_x_continuous(breaks = seq(-6, 1, by = 1), 
                     labels = c(two_neg_six, "", two_neg_four, "",
                                two_neg_two, "", two_zero, ""), 
                     limits = c(-6, 0)) +
  coord_cartesian(xlim = c(-6, 1),
                  ylim = c(-5, 4.58),
                  expand = FALSE) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6),
    axis.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    #panel.background = element_rect(fill = col_hb, colour = col_hb),
    legend.position="none"
  )

#### _pg_bottom_right (subpanel iv) ####
# data_1L_nc8wu3 <- 
#   read_csv("./Figure 4/panel_h_to_j/pigeon/nc8wu3.csv") %>%
#   mutate(FIT_NORM = FIT/max(FIT))

## right cell
file_1L_nc7wu3t <- "./data/spatemp/LM_ER_fullsweep_20bins/ALL_LM_SPATEMP/pigeon_LM_spatemp_min_err_gaussian/nc7wu3 redo fast june 2018.xls"

params_1L_nc7wu3t <-
  extract_gaussian_parameters(file_1L_nc7wu3t, method = "unconstrained")

grid_1L_nc7wu3t <-
  expand.grid(X_values = seq(from = -5, to = 0, by = 0.1),
              Y_values = seq(from = 2.5, to = 5, by = 0.1))
preds_1L_nc7wu3t <-
  predict_gaussian_xl(
    X_values = grid_1L_nc7wu3t$X_values,
    Y_values = grid_1L_nc7wu3t$Y_values,
    A = params_1L_nc7wu3t$A[1],
    X_peak = params_1L_nc7wu3t$X_peak[1],
    Y_peak = params_1L_nc7wu3t$Y_peak[1],
    Q = params_1L_nc7wu3t$Q[1],
    X_var = params_1L_nc7wu3t$X_var[1],
    Y_var = params_1L_nc7wu3t$Y_var[1]
  )

## left cell
file_1L_nc7wu3b <- "./data/spatemp/LM_ER_fullsweep_20bins/ALL_LM_SPATEMP/pigeon_LM_spatemp_min_err_gaussian/nc7wu3 redo slow june 2018.xls"

params_1L_nc7wu3b <-
  extract_gaussian_parameters(file_1L_nc7wu3b, method = "unconstrained")

grid_1L_nc7wu3b <-
  expand.grid(X_values = seq(from = -5, to = 0, by = 0.1),
              Y_values = seq(from = -5, to = 2.7, by = 0.1))
preds_1L_nc7wu3b <-
  predict_gaussian_xl(
    X_values = grid_1L_nc7wu3b$X_values,
    Y_values = grid_1L_nc7wu3b$Y_values,
    A = params_1L_nc7wu3b$A[1],
    X_peak = params_1L_nc7wu3b$X_peak[1],
    Y_peak = params_1L_nc7wu3b$Y_peak[1],
    Q = params_1L_nc7wu3b$Q[1],
    X_var = params_1L_nc7wu3b$X_var[1],
    Y_var = params_1L_nc7wu3b$Y_var[1]
  )

amprat_1L_nc7wu3 <- params_1L_nc7wu3t$A[1] / params_1L_nc7wu3b$A[1]

data_1L_nc7wu3t <- tibble(
  SFFIT = grid_1L_nc7wu3t$X_values,
  TFFIT = grid_1L_nc7wu3t$Y_values,
  FIT_NORM = amprat_1L_nc7wu3 * preds_1L_nc7wu3t / 
    max(preds_1L_nc7wu3t)
)

data_1L_nc7wu3b <- tibble(
  SFFIT = grid_1L_nc7wu3b$X_values,
  TFFIT = grid_1L_nc7wu3b$Y_values,
  FIT_NORM = preds_1L_nc7wu3b /
    max(preds_1L_nc7wu3b)
)

## to average the intersection of cells:
top_coords <- tibble(SFFIT = data_1L_nc7wu3t$SFFIT, 
                     TFFIT = data_1L_nc7wu3t$TFFIT)
bot_coords <- tibble(SFFIT = data_1L_nc7wu3b$SFFIT, 
                     TFFIT = data_1L_nc7wu3b$TFFIT)

inter <- dplyr::inner_join(top_coords, bot_coords)

## Extract the overlap and average
top_extract <- dplyr::inner_join(data_1L_nc7wu3t, inter)
bot_extract <- dplyr::inner_join(data_1L_nc7wu3b, inter)
extract_avg <- tibble((top_extract + bot_extract)/2)

## Extract the rest of the originals
top_retain <- dplyr::anti_join(data_1L_nc7wu3t, inter,
                               by = c("SFFIT", "TFFIT"))
bot_retain <- dplyr::anti_join(data_1L_nc7wu3b, inter, 
                               by = c("SFFIT", "TFFIT"))

comb_dat <- bind_rows(top_retain, extract_avg, bot_retain)

data_1L_nc7wu3 <- 
aggregate(. ~ SFFIT + TFFIT, data = comb_dat, FUN = mean) %>% tibble()
   
p_1L_iv <-
  ggplot(data_1L_nc7wu3, aes(SFFIT, TFFIT, z = FIT_NORM)) +
  metR::geom_contour_fill(#aes(fill = ..level..), 
                          bins = bins_common) +
  scale_fill_gradientn(
    colors=c(rgb(0, 0, 0, maxColorValue = 255), #black
             rgb(51, 153, 255, maxColorValue = 255), #blue
             rgb(255, 255, 255, maxColorValue = 255))
  ) +
  xlab(" ") +
  #ylab(" ") +
  expand_limits(x = c(-6, 1),
                y = c(-5, 4.58)) +
  scale_y_continuous(breaks = c(-5, -3, -1,  1,  3,  4, 4.58),
                     labels = c("", two_neg_three, "", two_one,
                                "", two_four, "")) +
  scale_x_continuous(breaks = seq(-6, 1, by = 1), 
                     labels = c("", two_neg_five, "", two_neg_three,
                                "", two_neg_one, "", two_one), 
                     limits = c(-6, 1)) +
  coord_cartesian(xlim = c(-6, 1),
                  ylim = c(-5, 4.58),
                  expand = FALSE) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 6),
    axis.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    legend.position="none"
  )

####### ARRANGEMENT #######
panel_1J_HB_top <- plot_grid(p_1J_i, NULL, p_1J_ii,
                             ncol = 3, #align = "tblr",
                             rel_widths = c(1, 0, 0.85))
panel_1J_HB_bot <- plot_grid(p_1J_iii, NULL, p_1J_iv,
                             ncol = 3, #align = "tblr",
                             rel_widths = c(1, 0, 0.85))
panel_1J_HB <- plot_grid(panel_1J_HB_top, panel_1J_HB_bot,
                         ncol = 1,
                         rel_heights = c(0.85, 1))

panel_1K_ZF_top <- plot_grid(p_1K_i, NULL, p_1K_ii,
                             ncol = 3, #align = "tblr",
                             rel_widths = c(1, 0, 0.85))
panel_1K_ZF_bot <- plot_grid(p_1K_iii, NULL, p_1K_iv,
                             ncol = 3, #align = "tblr",
                             rel_widths = c(1, 0, 0.85))
panel_1K_ZF <- plot_grid(panel_1K_ZF_top, panel_1K_ZF_bot,
                         ncol = 1,
                         rel_heights = c(0.85, 1))

panel_1L_PG_top <- plot_grid(p_1L_i, NULL, p_1L_ii,
                             ncol = 3, #align = "tblr",
                             rel_widths = c(1, 0, 0.85))
panel_1L_PG_bot <- plot_grid(p_1L_iii, NULL, p_1L_iv,
                             ncol = 3, #align = "tblr",
                             rel_widths = c(1, 0, 0.85))
panel_1L_PG <- plot_grid(panel_1L_PG_top, panel_1L_PG_bot,
                         ncol = 1,
                         rel_heights = c(0.85, 1))

panel_1JKL <- plot_grid(NULL, panel_1J_HB, NULL, panel_1K_ZF, NULL, panel_1L_PG,
                        nrow = 1,
                        labels = c(" ", " ", " ", "K", " ", "L"),
                        label_x = c(rep(-0.08, 3), -0.05, -0.08, -0.08),
                        label_size = panel_font_size,
                        label_fontfamily = "sans",
                        rel_widths = c(0, 0.9, 0.05, 0.9, 0.05, 0.9))
