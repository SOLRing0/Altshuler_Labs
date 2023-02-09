
## Functions
## Note from VBB: for this function to work, peak location, X- and Y- variance,
## amplitude, and Q must all be pre-determined for a 2D gaussian. This function
## then produces the predicted height at each combination of X and Y This is a
## preliminary version of what ultimately became
## gaussplotR::predict_gaussian_2D()
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

## Options
## ggplot margins
## t, r, b, l
plot_margs <- unit(c(0.1, 0.1, 0.1, 0.1), "cm")

## Breaks
breaks_common = seq(0.1, 1, by = 0.1)
bins_common = 20

###### Panel G - Peristimulus time histogram #####
data_1G_matlab <-
  R.matlab::readMat("./data/spatemp/LM_ER_fullsweep_example_data/CALSTC09-LM-02-3001-ST_PSTH.mat")

data_1G_plotdat <- data_1G_matlab$ER.psth.data[,4:23]


data_1G_max_init <-
  apply(data_1G_plotdat, 2, function(x)
    max(x, na.rm = TRUE))
data_1G_max_max <- max(data_1G_max_init)
data_1G_max_corrected <-
  data_1G_max_max - data_1G_matlab$mean.baseline[1]
data_1G_basearray <- c(rep(data_1G_matlab$mean.baseline[1], 20))
data_1G_psth_letters <- letters[1:20]

data_1G_plotz <- NULL
data_1G_plt_subbase <- NULL
panel_1G_datlist <- NULL
for (i in 1:nrow(data_1G_plotdat)) {
  data_1G_plt_subbase[[i]] <- data_1G_plotdat[i, ] - data_1G_basearray
  panel_1G_datlist[[i]] <- tibble(lets = data_1G_psth_letters,
                                  psth_dat = data_1G_plt_subbase[[i]])
  data_1G_plotz[[i]] <-
    ggplot(panel_1G_datlist[[i]], aes(x = lets, y = psth_dat)) +
    # col_hb = "#ED0080"
    geom_col(fill = "#ED0080",
             color = "#ED0080",
             size = 0) +
    scale_y_continuous(breaks = seq(0, data_1G_max_corrected, by = 10)) +
    coord_cartesian(ylim = c(0, data_1G_max_corrected),
                    expand = FALSE) +
    theme_classic() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(colour = "black", size = 0.1),
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      axis.ticks = element_line(colour = "black", size = 0.1),
      axis.ticks.x = element_blank(),
      axis.ticks.length = unit(0.05, "cm")
    )
}

panel_1G <- plot_grid(plotlist = rev(data_1G_plotz),
                      nrow = 6)

#### Data for H and I #####
data_1HI <-
  read_csv(
    "./data/spatemp/LM_ER_fullsweep_example_data/CALSTC09-LM-02-3001-Cell 1 - export.csv"
  )  %>%
  mutate(FIT_NORM = G_result_Un / max(G_result_Un)) %>%
  drop_na()

###### Panel H - contour plot #####
## Options

## ggplot margins
## t, r, b, l
plot_margs <- unit(c(0.1, 0.1, 0.1, 0.1), "cm")

panel_1H <-
  ggplot(data_1HI, aes(x = Log2SF, y = Log2TF, z = R_Norm)) + 
  metR::geom_contour_fill(#aes(fill = ..level..), 
                          color = "black", size = 0.03,
                          bins = 12) +
  ## greyscale:
    # scale_fill_gradientn(
    #   colors=c(rgb(0, 0, 0, maxColorValue = 255), #black
    #            rgb(255, 255, 255, maxColorValue = 255))
    # ) +
  ## hummingbird color:
  scale_fill_gradientn(
    colors=c(rgb(0, 0, 0, maxColorValue = 255), #black
             rgb(255, 51, 153, maxColorValue = 255), #pink
             rgb(255, 255, 255, maxColorValue = 255))
  ) +
  xlab("spatial frequency (cpd)") +
  ylab("temporal frequency (Hz)") +
  expand_limits(x = c(-6, 1),
                y = c(-5, 4.58)) +
  scale_y_continuous(breaks = c(-5, -3, -1, 1, 3, 4),
                     labels = c("0.031", "0.125", "0.5", "2", "8", "16")) +
  scale_x_continuous(breaks = seq(-6, 1, by = 1), 
                     labels = c("0.0155", " ", "0.062",
                                " ", "0.25", " ",
                                "1", "")) +
  coord_cartesian(xlim = c(-6, 1), 
                  ylim = c(-5, 4.58),
                  expand = FALSE) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 6),
    axis.title.y = element_text(hjust = 0.40),
    axis.title = element_text(size = 7),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    legend.position="none"
  )

##### _4B legend #####
legend_1H <- data.frame(x = rep(1, 12),
                        y = 1:12)

panel_1H_legend <- 
  legend_1H %>%
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = y),
            color = "black") +
  scale_x_continuous(limits = c(0, 2), breaks = 1)+
  scale_fill_gradient2(low = rgb(0, 0, 0, maxColorValue = 255),
                       mid = rgb(255, 51, 153, maxColorValue = 255),
                       high = rgb(255, 255, 255, maxColorValue = 255),
                       midpoint = 6) +
  theme_void() +
  theme(legend.position = "none")

###### Panel I - gaussian plot #####
file_1I <-
  "./data/spatemp/LM_ER_fullsweep_example_data/CALSTC09-LM-02-3001-Cell 1 redone july 2020.xls"

params_1I <- extract_gaussian_parameters(file_1I, method = "unconstrained")

grid_1I <- expand.grid(X_values = seq(from = -6, to = 0, by = 0.1),
                       Y_values = seq(from = -5, to = 4, by = 0.1))
preds_1I <- predict_gaussian_xl(X_values = grid_1I$X_values,
                                Y_values = grid_1I$Y_values,
                                A = params_1I$A[1],
                                X_peak = params_1I$X_peak[1],
                                Y_peak = params_1I$Y_peak[1],
                                Q = params_1I$Q[1],
                                X_var = params_1I$X_var[1],
                                Y_var = params_1I$Y_var[1])


plotdat_1I <- data.frame(SFFIT = grid_1I$X_values,
                         TFFIT = grid_1I$Y_values,
                         FIT_NORM = preds_1I/max(preds_1I))

p_1I <-
  plotdat_1I %>%
  ggplot( aes(SFFIT, TFFIT, z = FIT_NORM)) + 
  metR::geom_contour_fill(#aes(fill = ..level..), 
                          bins = 15, na.fill = TRUE) +
  scale_fill_gradientn(
    colors=c(rgb(0, 0, 0, maxColorValue = 255), #black
             rgb(255, 51, 153, maxColorValue = 255), #pink
             rgb(255, 255, 255, maxColorValue = 255))
  ) +
  xlab("spatial frequency (cpd)") +
  ylab(" ") +
  expand_limits(x = c(-6, 1),
                y = c(-5, 4.58)) +
  scale_y_continuous(breaks = c(-5, -3, -1, 1, 3, 4),
                     labels = c(two_neg_five, two_neg_three, 
                                two_neg_one, two_one, two_three,
                                two_four)) +
  scale_x_continuous(breaks = seq(-6, 1, by = 1), 
                     labels = c(expression(paste("   ", 2^-6)), 
                                " ", two_neg_four,
                                " ", two_neg_two, " ",
                                two_zero, "")) +
  coord_cartesian(xlim = c(-6, 1), 
                  ylim = c(-5, 4.58),
                  expand = FALSE) +
  geom_abline(slope = 1, intercept = 11,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = 9,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = 7,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = 5,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = 3,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = 1,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = -1,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = -3,  colour = "grey30", size = 0.2) +
  theme_classic() +
  theme(
    # axis.text.y = element_text(size = 6),
    # axis.text.x = element_text(size = 6, margin = margin(r = 0)),
    # axis.title = element_text(size = 7),
    # plot.margin = unit(c(0.1, 0.1, 0.06, 0.1), "cm"), #trbl
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6, hjust = 0), 
    axis.title = element_text(size = 7),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs, #trbl
    #panel.background = element_rect(fill = col_hb, colour = col_hb),
    legend.position="none"
  )

# q_1I <- ggplot_build(p_1I)
# q_1I_colors <- q_1I$data[[1]]$fill 
# q_1I$data[[1]]$fill <- 
#   replace(q_1I_colors,
#           q_1I_colors == unique(q_1I$data[[1]]$fill)[1],
#           col_hb) 
# 
# panel_1I <- ggplot_gtable(q_1I)

p_1I -> panel_1I
#panel_1I


##### _1I legend #####
legend_1I <- data.frame(x = rep(1, 15),
                        y = 1:15)

panel_1I_legend <- 
  legend_1I %>%
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = y),
            color = "black") +
  scale_x_continuous(limits = c(0, 2), breaks = 1)+
  scale_fill_gradient2(low = rgb(0, 0, 0, maxColorValue = 255),
                       mid = rgb(255, 51, 153, maxColorValue = 255),
                       high = rgb(255, 255, 255, maxColorValue = 255),
                       midpoint = 8) +
  theme_void() +
  theme(legend.position = "none")
