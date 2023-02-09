bins_common = 15

p_4D2F_cell <- "STC09-LM-08-5154-cell1"

p_4D2F_fullsweep_dir <-
  "./data/spatemp/LM_ER_fullsweep_200bins_bin001_to_bin200/STC09-LM-08-5154-ST-STER_data.csv"
p_4D2F_initrns_dir <-
  "./data/spatemp/LM_ER_transients_200bins_bin004_to_bin020/STC09-LM-08-5154_transients.csv"
p_4D2F_steady_dir <-
  "./data/spatemp/LM_ER_steadystate_200bins_bin100_to_bin200/STC09-LM-08-5154_steadystate.csv"

###### import cell data ####
p_4D2F_fullsweep_dat <-
  read_csv(p_4D2F_fullsweep_dir, col_names = FALSE)
p_4D2F_fullsweep_bsln <-
  read_csv(
    "./data/spatemp/LM_ER_fullsweep_200bins_bin001_to_bin200/STC09-LM-08-5154-ST-spontaneous_rate.csv", col_names = FALSE
  )
p_4D2F_initrns_dat <- read_csv(p_4D2F_initrns_dir, col_names = TRUE)
p_4D2F_steady_dat <- read_csv(p_4D2F_steady_dir, col_names = TRUE)

##### format for gaussplotR #####
p_4D2F_fullsweep_gpr <-
  p_4D2F_fullsweep_dat %>%
  select(X3, X2, X4) %>%
  rename(raw_tf = X2, raw_sf = X3, cell1 = X4) %>%
  transmute(X_values = log2(raw_sf/0.0424),
            Y_values = log2(raw_tf),
            response = cell1/max(cell1)) %>%
  drop_na() %>%
  as.data.frame()
p_4D2F_initrns_gpr <-
  p_4D2F_initrns_dat %>%
  select(logSF, logTF, ends_with("cell1_norm")) %>%
  rename(
    X_values = logSF,
    Y_values = logTF,
    response = cell1_norm) %>%
  drop_na() %>%
  as.data.frame()
p_4D2F_steady_gpr <-
  p_4D2F_steady_dat %>%
  select(logSF, logTF, ends_with("cell1_norm")) %>%
  rename(
    X_values = logSF,
    Y_values = logTF,
    response = cell1_norm) %>%
  drop_na() %>%
  as.data.frame()

#### run fit 2D gaussians ####
p_4D2F_fullsweep_gauss <-
  gaussplotR::fit_gaussian_2D(
    p_4D2F_fullsweep_gpr, method = "elliptical_log",
    constrain_amplitude = TRUE,
    constrain_orientation = "unconstrained", 
    maxiter = 10000,
    minFactor = 0.000000488281
  )
p_4D2F_initrns_gauss <-
  gaussplotR::fit_gaussian_2D(
    p_4D2F_initrns_gpr, method = "elliptical_log",
    constrain_amplitude = TRUE,
    constrain_orientation = "unconstrained", 
    maxiter = 10000,
    minFactor = 0.000000488281
  )
p_4D2F_steady_gauss <-
  gaussplotR::fit_gaussian_2D(
    p_4D2F_steady_gpr, method = "elliptical_log",
    constrain_amplitude = TRUE,
    constrain_orientation = "unconstrained", 
    maxiter = 10000,
    minFactor = 0.000000488281
  )

## set up prediction space
grid_4D2F <-
  expand.grid(X_values = seq(from = -5, to = 0, by = 0.1),
              Y_values = seq(from = -5, to = 4, by = 0.1))

## predictions
preds_4D2F_fullsweep <-
  predict_gaussian_2D(
    p_4D2F_fullsweep_gauss,
    X_values = grid_4D2F$X_values,
    Y_values = grid_4D2F$Y_values
  )
preds_4D2F_initrns <-
  predict_gaussian_2D(
    p_4D2F_initrns_gauss,
    X_values = grid_4D2F$X_values,
    Y_values = grid_4D2F$Y_values
  )
preds_4D2F_steady <-
  predict_gaussian_2D(
    p_4D2F_steady_gauss,
    X_values = grid_4D2F$X_values,
    Y_values = grid_4D2F$Y_values
  )

##### plots #####
panel_4D <-
  ggplot(preds_4D2F_fullsweep, aes(X_values, Y_values, z = predicted_values)) + 
  metR::geom_contour_fill(#aes(fill = ..level..), 
    size = contour_wds,
    bins = bins_common) +
  scale_fill_gradientn(
    colors = c(
      rgb(0, 0, 0, maxColorValue = 255), #black
      rgb(255, 128, 0, maxColorValue = 255), #orange
      rgb(255, 255, 255, maxColorValue = 255)
    )) + 
  xlab("spatial frequency (cpd)") +
  ylab(expression(paste("temporal frequency (Hz)"))) +
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
  geom_abline(slope = 1, intercept = 11,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = 9,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = 7,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = 5,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = 3,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = 1,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = -1,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = -3,  colour = "grey30", size = 0.2) +
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

panel_4E <-
  ggplot(preds_4D2F_initrns, aes(X_values, Y_values, z = predicted_values)) + 
  metR::geom_contour_fill(#aes(fill = ..level..), 
    size = contour_wds,
    bins = bins_common) +
  scale_fill_gradientn(
    colors = c(
      rgb(0, 0, 0, maxColorValue = 255), #black
      rgb(255, 128, 0, maxColorValue = 255), #orange
      rgb(255, 255, 255, maxColorValue = 255)
    )) + 
  xlab(" ") +
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
  geom_abline(slope = 1, intercept = 11,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = 9,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = 7,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = 5,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = 3,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = 1,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = -1,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = -3,  colour = "grey30", size = 0.2) +
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

panel_4F <-
  ggplot(preds_4D2F_steady, aes(X_values, Y_values, z = predicted_values)) + 
  metR::geom_contour_fill(#aes(fill = ..level..), 
    size = contour_wds,
    bins = bins_common) +
  scale_fill_gradientn(
    colors = c(
      rgb(0, 0, 0, maxColorValue = 255), #black
      rgb(255, 128, 0, maxColorValue = 255), #orange
      rgb(255, 255, 255, maxColorValue = 255)
    )) + 
  xlab(" ") +
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
  geom_abline(slope = 1, intercept = 11,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = 9,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = 7,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = 5,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = 3,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = 1,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = -1,  colour = "grey30", size = 0.2) +
  geom_abline(slope = 1, intercept = -3,  colour = "grey30", size = 0.2) +
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

