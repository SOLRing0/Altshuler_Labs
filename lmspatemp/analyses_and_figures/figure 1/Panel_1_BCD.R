## polar histograms of LM preferred directions

#### Functions ####
sem <- function(x){
  sd(x)/sqrt(length(x))
}  

get_prefdir <- function(file_name, species = "unknown"){
  
  ## get single value for mean spontaneous rate
  baseline <- 
    read_csv(file_name, col_names = FALSE, 
             col_types = cols(.default = col_double())) %>% 
    filter(X1 == 'NaN') %>% ## X1 is direction X2 is speed
    select(-X1, -X2) %>% 
    pivot_longer(everything()) %>% 
    summarise(mean_baseline = mean(value)) %>% 
    pull()  

  ## Get mean firing rate for each direction trial. 
  ## Stim_per indicates the order in which it was tested during the experiment
  firing_frame <- 
    read_csv(file_name, col_names = FALSE, 
             col_types = cols(.default = col_double())) %>% 
    rename(direction = X1, speed = X2) %>% 
    filter(direction != 'NaN') %>% 
    mutate(stim_per = row_number()) %>% 
    pivot_longer(contains("X")) %>%  ## this selects all the 250ms columns
    group_by(direction, speed, stim_per) %>% 
    summarise(firing = mean(value - baseline)) %>% ## subtract baseline 
    ungroup() %>% 
    # need equal number of replicates for each direction. The modulus gives us
    # the number of rows (x) more than a complete set of directions. 
    # Exclude the last x trials, so that each direction is equally weighted
    filter(
      stim_per <= max(stim_per) - (max(stim_per) %% n_distinct(direction))) %>%
    group_by(direction, speed) %>%
    summarise(firing = mean(firing)) %>%
    ungroup() 
  
  ## Exclude file if all the firing rates are too close to baseline. 
  ## This looks to see if all firing rates occur outside 0.5 * baseline and 1.5 
  ## * baseline -- if they do, then the cell is included
  firing_frame$firing_plus_bsl <- firing_frame$firing + baseline
  
  if(
    firing_frame %>%
    pull(firing_plus_bsl) %>%
    {any(
      . < 0.8 * baseline,
      . > 1.2 * baseline
    )}
  ) {
    print(str_replace(file_name, ".*//(.*)", "\\1"))
  } else{
    #print("this cell is a dad")
    return(NULL)
  }
  
  firing_frame_shift <- firing_frame %>%
    mutate(firing = 
      case_when(min(firing) < 0 ~ firing + abs(min(firing)), TRUE ~ firing))
  
  circ_data <- firing_frame_shift %>%
    pmap(
      function(direction, firing, ...) { 
        # ... captures unused columns from firing_frame
        rep(direction, round(firing * 10)) 
        # multiply by 10 to increase resolution
      }
    ) %>%
    unlist() %>%
    circular(units = "degrees", rotation = "clock")
  
  if(rayleigh.test(circ_data)$p.value > 0.05){
    warning("pvalue greater than 0.05")
    return(NULL)
  }
  
  vector_sum <- firing_frame_shift %>% 
    group_by(direction) %>%
    summarize(firing = mean(firing)) %>%
    transmute(
      x = cos(direction * pi/180) * firing,
      y = sin(direction * pi/180) * firing
    ) %>% 
    summarise_all(mean) %>% 
    transmute(prefdir_sum = (atan2(y,x) * 180/pi) %% 360)
  
  si_value <- firing_frame_shift %>% 
    group_by(direction) %>% 
    summarize(firing = mean(firing)) %>% 
    transmute(
      a = (sin(direction * pi/180) * firing),
      b = (cos(direction * pi/180) * firing),
      c = mean(firing)
    ) %>% 
    summarise_all(mean) %>% 
    transmute(si = sqrt(a^2 + b^2)/c)
  
  vm <- circ_data %>%
    mle.vonmises %>%
    unlist %>%
    `[`(c('mu','kappa')) %>%
    as_tibble %>%
    mutate(mu = mu %% 360)
  
  firing_by_angle <- firing_frame_shift %>%
    group_by(direction) %>%
    summarize(value = mean(firing)) %>%
    pivot_wider(everything(), 
                names_from = direction, 
                names_prefix = "firing_at_")
  
  return(
    tibble(
      filename = str_replace(file_name, ".*//(.*)", "\\1"),
      bird_num = str_replace(file_name, ".*//[A-Za-z]+([0-9]+).*", "\\1"),
      track = str_replace(file_name, ".*(trk|nbor-)([0-9]+).*", "\\2"),
      site = str_replace(file_name, ".*(trk|nbor-)[0-9]+-([0-9]+).*", "\\2"),
      cell = str_replace(file_name, ".*cell([0-9]+).*", "\\1"),
      species = species,
      baseline = baseline
    ) %>%
      bind_cols(
        vector_sum,
        si_value,
        vm,
        firing_by_angle
      )
  )
}

## GS's files have names formatted differently, hence a separate get_prefdir()
## function. Key difference is at the end, where regex patterns are revised
## to deal with GS's file naming system.
gs_get_prefdir <- function(file_name, species = "unknown"){
  ## get single value for mean spontaneous rate
  baseline <- read_csv(file_name, col_names = FALSE, 
                       col_types = cols(.default = col_double())) %>% 
    filter(X1 == 'NaN') %>% ## X1 is direction X2 is speed
    select(-X1, -X2) %>% 
    pivot_longer(everything()) %>% 
    summarise(mean_baseline = mean(value)) %>% 
    pull()  

  firing_frame <- 
    read_csv(file_name, col_names = FALSE, 
             col_types = cols(.default = col_double())) %>% 
    rename(direction = X1, speed = X2) %>% 
    filter(direction != 'NaN') %>% 
    mutate(stim_per = row_number()) %>% 
    pivot_longer(contains("X")) %>%  ## this selects all the 250ms columns
    group_by(direction, speed, stim_per) %>% 
    summarise(firing = mean(value - baseline)) %>% ## subtract baseline 
    ungroup() %>% 
    # need equal number of replicates for each direction. The modulus gives us
    # the number of rows (x) more than a complete set of directions. 
    # Exclude the last x trials, so that each direction is equally weighted
    filter(
      stim_per <= max(stim_per) - (max(stim_per) %% n_distinct(direction))) %>%
    group_by(direction, speed) %>%
    summarise(firing = mean(firing)) %>%
    ungroup() 
  
  ## Exclude file if all the firing rates are too close to baseline. 
  ## This looks to see if all firing rates occur outside 0.5 * baseline and 1.5 
  ## * baseline -- if they do, then the cell is included
  
  firing_frame$firing_plus_bsl <- firing_frame$firing + baseline
  
  if(
    firing_frame %>%
    pull(firing_plus_bsl) %>%
    {any(
      . < 0.8 * baseline,
      . > 1.2 * baseline
    )}
  ) {
    print(str_replace(file_name, ".*//(.*)", "\\1"))
  } else{
    print("this cell is a dad")
    return(NULL)
  }
  
  firing_frame_shift <- firing_frame %>%
    mutate(firing = case_when(
      min(firing) < 0 ~ firing + abs(min(firing)), TRUE ~ firing))
  
  circ_data <- firing_frame_shift %>%
    pmap(
      function(direction, firing, ...) {
        # ... captures unused columns from firing_frame
        rep(direction, round(firing * 10)) 
        # multiply by 10 to increase resolution
      }
    ) %>%
    unlist() %>%
    circular(units = "degrees", rotation = "clock")
  
  if(rayleigh.test(circ_data)$p.value > 0.05){
    warning("pvalue greater than 0.05")
    return(NULL)
  }
  
  vector_sum <- firing_frame_shift %>% 
    group_by(direction) %>%
    summarize(firing = mean(firing)) %>%
    transmute(
      x = cos(direction * pi/180) * firing,
      y = sin(direction * pi/180) * firing
    ) %>% 
    summarise_all(mean) %>% 
    transmute(prefdir_sum = (atan2(y,x) * 180/pi) %% 360)
  
  vm <- circ_data %>%
    mle.vonmises %>%
    unlist %>%
    `[`(c('mu','kappa')) %>%
    as_tibble %>%
    mutate(mu = mu %% 360)
  
  firing_by_angle <- firing_frame_shift %>%
    group_by(direction) %>%
    summarize(value = mean(firing)) %>%
    pivot_wider(everything(), names_from = direction,
                names_prefix = "firing_at_")
  
  global_max_firing <- max(firing_frame$firing)
  
  return(
    tibble(
      filename = str_replace(file_name, ".*//(.*)", "\\1"),
      bird_num = str_split(filename, "-")[[1]][1],
      track = str_split(filename, "-")[[1]][3],
      site = str_split(filename, "-")[[1]][4],
      cell = str_replace(filename, ".*cell-([0-9]+).*", "\\1"),
      species = species,
      baseline = baseline,
      global_max_firing = global_max_firing
    ) %>%
      bind_cols(
        vector_sum,
        vm,
        firing_by_angle
      )
  )
}

get_vm_curve <- function(filename, species, kappa, 
                         steps = seq(-180, 180, by = 5), ...){
  tibble(
    filename = filename,
    species = species,
    y_temp = dvonmises(x = steps * pi / 180, mu = 0, kappa = kappa), 
    #suppressWarnings(
    y = y_temp + (1 - max(y_temp)),
    x = steps,
    kappa = kappa
  ) %>%
    select(-y_temp)
}

#### Theme ####
polar_cell_theme_grp <-   theme_minimal() +
  theme(axis.title = element_blank(),
        #panel.border = element_blank(),
        legend.key = element_blank(),
        axis.ticks.x = element_blank(),
        # axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(), 
        #element_line(color = "gray90", size = 0.3),
        plot.title = element_blank())

#### Colors ####
lt_pink <- "#fccae4"
pink <-  "#ffbddf"
lt_magenta  <-  "#fc4ea9"
lt_yellow <- "#fff799"
yellow <- "#ffe352"
lt_orange <- "#FFE4B5"

###### Import LM data for hummingbirds and zebra finches ####
## Analyses will be done on all LM cells from hummingbirds and zebra finches
## and then later separated out to create species-specific plots

#### __AG data first ####
hb_file_list_LM_ag <- 
  list.files(path = "./data/direction/2015-CALAN-PSTH-direction/",
             full.names = TRUE) 

zf_file_list_LM_ag <- 
  list.files(path = "./data/direction/2015-ZF-PSTH-direction/",
             full.names = TRUE) 

hb_prefdir_ag <- map_dfr(hb_file_list_LM_ag, get_prefdir, species = "hb")
zf_prefdir_ag <- map_dfr(zf_file_list_LM_ag, get_prefdir, species = "zf")

prefdir_df_lm_ag <- bind_rows(zf_prefdir_ag, hb_prefdir_ag) 

## add file path to data frame
filepath_hb_lm_ag <- 
  list.files(path = "./data/direction/2015-CALAN-PSTH-direction/",
             full.names = F) %>% 
  cbind(file_path = list.files(
    path = "./data/direction/2015-CALAN-PSTH-direction/",
    full.names = T))

filepath_zf_lm_ag <- 
  list.files(path = "./data/direction/2015-ZF-PSTH-direction/",
             full.names = F) %>% 
  cbind(file_path = list.files(
    path = "./data/direction/2015-ZF-PSTH-direction/",
    full.names = T))

filepath_all_lm_ag <- rbind(filepath_zf_lm_ag, filepath_hb_lm_ag)
filepath_all_lm_ag <- 
  as_tibble(filepath_all_lm_ag) %>% 
  rename(filename = ".")

calc_tw_lm_ag <- left_join(prefdir_df_lm_ag, 
                           filepath_all_lm_ag, by = "filename")

## calculate widths of direction tuning curves
dir_width_lm_ag <- NULL
dir_spline_fits_lm_ag <- NULL
source("./data/direction/loglik_smooth_spline.R")
for (i in 1:dim(calc_tw_lm_ag)[1]) {
  one_cell <-  calc_tw_lm_ag$file_path[i] 
  prefdir_sum <- calc_tw_lm_ag$prefdir_sum[i]
  
  baseline <- 
    suppressMessages(read_csv(one_cell, col_names = FALSE)) %>% 
    filter(X1 == 'NaN') %>% 
    select(-X1, -X2) %>% 
    pivot_longer(everything()) %>% 
    summarise(mean_baseline = mean(value)) %>% 
    pull()
  
  sem_baseline <- 
    suppressMessages(read_csv(one_cell, col_names = FALSE)) %>% 
    filter(X1 == 'NaN') %>% 
    select(-X1, -X2) %>% 
    pivot_longer(everything()) %>% 
    summarise(sem_baseline = sem(value)) %>% 
    pull()
  
  test_cell_mean <- 
    suppressMessages(read_csv(one_cell, 
                              col_names = FALSE, 
                              col_types = cols(.default = col_double()))) %>% 
    rename(direction = X1, speed = X2) %>% 
    filter(direction != 'NaN') %>% 
    mutate(stim_per = row_number()) %>% 
    pivot_longer(contains("X")) %>%  ## this selects all the 250ms columns
    group_by(direction, speed, stim_per) %>% 
    summarise(firing = mean(value - baseline)) %>% ## subtract baseline 
    ungroup() %>% 
    # need equal number of replicates for each direction. The modulus gives us
    # the number of rows (x) more than a complete set of directions. 
    # Exclude the last x trials, so that each direction is equally weighted
    filter(
      stim_per <= max(stim_per) - (max(stim_per) %% n_distinct(direction))) %>%
    group_by(direction, speed) %>%
    summarise(firing = mean(firing)) %>%
    ungroup() 
  
  test_cell_raw <- 
    suppressMessages(read_csv(one_cell, 
                              col_names = FALSE, 
                              col_types = cols(.default = col_double()))) %>% 
    rename(direction = X1, speed = X2) %>% 
    filter(direction != 'NaN') %>% 
    mutate(stim_per = row_number()) %>% 
    pivot_longer(contains("X")) %>%  ## this selects all the 250ms columns
    group_by(direction, speed, stim_per) %>% 
    summarise(firing = mean(value - baseline)) %>% ## subtract baseline 
    ungroup() %>% 
    # need equal number of replicates for each direction. The modulus gives us
    # the number of rows (x) more than a complete set of directions. 
    # Exclude the last x trials, so that each direction is equally weighted
    filter(stim_per <= max(stim_per) - (max(stim_per) %% n_distinct(direction)))
  
  ## Adding two additional loops to ensure wrapping between 0 and 360
  x_raw <- c(test_cell_raw$direction-360, 
             test_cell_raw$direction, test_cell_raw$direction + 360)
  y_raw <- c(test_cell_raw$firing, test_cell_raw$firing, test_cell_raw$firing)
  
  fit5 <- smooth.spline(x = x_raw, y = y_raw, df = 5)
  fit6 <- smooth.spline(x = x_raw, y = y_raw, df = 6)
  fit7 <- smooth.spline(x = x_raw, y = y_raw, df = 7)
  fit8 <- smooth.spline(x = x_raw, y = y_raw, df = 8)
  fit9 <- smooth.spline(x = x_raw, y = y_raw, df = 9)
  fit10 <- smooth.spline(x = x_raw, y = y_raw, df = 10)
  fit11 <- smooth.spline(x = x_raw, y = y_raw, df = 11)
  fit12 <- smooth.spline(x = x_raw, y = y_raw, df = 12)
  fit13 <- smooth.spline(x = x_raw, y = y_raw, df = 13)
  fit14 <- smooth.spline(x = x_raw, y = y_raw, df = 14)
  fit15 <- smooth.spline(x = x_raw, y = y_raw, df = 15)
  fit16 <- smooth.spline(x = x_raw, y = y_raw, df = 16)
  fit17 <- smooth.spline(x = x_raw, y = y_raw, df = 17)
  fit18 <- smooth.spline(x = x_raw, y = y_raw, df = 18)
  fit19 <- smooth.spline(x = x_raw, y = y_raw, df = 19)
  fit20 <- smooth.spline(x = x_raw, y = y_raw, df = 20)
  
  fit_list <- list("5" = fit5, "6" = fit6, "7" = fit7, 
                   "8" = fit8, "9" = fit9, "10" = fit10, 
                   "11" = fit11, "12" = fit12, "13" = fit13, 
                   "14" = fit14, "15" = fit15, "16" = fit16, 
                   "17" = fit17, "18" = fit18, "19" = fit19, 
                   "20" = fit20)
  
  best_fit <- map_dfr(fit_list, AICc_smooth_spline) %>% 
    pivot_longer(cols = everything(), 
                 names_to = "fit", 
                 values_to = "AICc") %>% 
    #mutate(fit = as.numeric(fit)) %>% 
    filter(AICc == min(AICc))
  
  best_fit_contents <- fit_list[best_fit$fit][[1]]
  
  y_vals <- stats:::predict.smooth.spline(best_fit_contents, 
                                          x = seq(0, 359, by = 1))
  predicted_values <- tibble(x = y_vals$x, y = y_vals$y)
  
  predicted_raw_vals <- stats:::predict.smooth.spline(best_fit_contents, 
                                                      x = seq(-360, 720, 1))
  predicted_raw_df <- tibble(x = predicted_raw_vals$x, 
                             y = predicted_raw_vals$y)
  
  max_y <- predicted_raw_df %>% 
    filter(x >= 0 & x < 360) %>% 
    filter(y == max(y))
  
  min_y <- predicted_raw_df %>% 
    filter(x >= 0 & x < 360) %>% 
    filter(y == min(y))
  
  y_range <- max_y$y - min_y$y
  
  lower_curve_of_interest <- predicted_raw_df %>% 
    filter(x >= max_y$x - 180 & x<= max_y$x)
  
  upper_curve_of_interest <- predicted_raw_df %>% 
    filter(x >= max_y$x & x<= max_y$x + 180)
  
  predicted_curve_center <- predicted_raw_df %>% 
    filter(x >= 0 & x < 360) %>% 
    rename(spline_x = x, spline_y = y) %>% 
    mutate(
      shifted_spline_y = spline_y- min(spline_y),
      standardized_x = spline_x/360,
    standardized_y = (spline_y - min(spline_y))/(max(spline_y) - min(spline_y)))
  
  translated_x <- predicted_raw_df %>% 
    filter(x >= max_y$x - 180 & x < max_y$x + 180) %>% 
    rename(trans_spline_x = x, trans_spline_y = y) %>% 
    mutate(
      trans_std_x = trans_spline_x - max_y$x, 
      trans_std_y = (trans_spline_y - min(trans_spline_y))/(max(trans_spline_y) - min(trans_spline_y)))
  
  firing_ratio_y <- predicted_raw_df  %>% 
    filter(x >= max_y$x - 180 & x < max_y$x + 180) %>% 
    mutate(bsln_ratio_y = (y + baseline)/baseline) %>% 
    select(bsln_ratio_y)
  
  ### If using prefdir_sum instead of spline peak (optional)
  prefdir_lower_curve_of_interest <- predicted_raw_df %>% 
    filter(x >= prefdir_sum - 180 & x<= prefdir_sum)
  
  prefdir_upper_curve_of_interest <- predicted_raw_df %>% 
    filter(x >= prefdir_sum & x<= prefdir_sum + 180)
  
  # predicted_curve_center_prefdir <- predicted_raw_df %>% 
  #   filter(x >= 0 & x <= 360) %>% 
  #   rename(spline_x = x, spline_y = y) %>% 
  #   mutate(shifted_spline_y = spline_y- min(spline_y), standardized_x = spline_x/360, standardized_y = (spline_y - min(spline_y))/(max(spline_y) - min(spline_y)))
  
  prefdir_translated_x <- predicted_raw_df %>% 
    filter(x >= prefdir_sum - 180 & x <= prefdir_sum + 180) %>% 
    rename(prefdir_trans_spline_x = x, prefdir_trans_spline_y = y) %>% 
    mutate(
      prefdir_trans_std_x = prefdir_trans_spline_x - prefdir_sum, 
      prefdir_trans_std_y = (prefdir_trans_spline_y - min(prefdir_trans_spline_y))/(max(prefdir_trans_spline_y) - min(prefdir_trans_spline_y)))
  
  # firing_ratio_y <- predicted_raw_df  %>% 
  #   filter(x >= max_y$x - 180 & x <= max_y$x + 180) %>% 
  #   mutate(bsln_ratio_y = (y + baseline)/baseline) %>% 
  #   select(bsln_ratio_y)
  
  si_value <- predicted_curve_center %>% 
    transmute(
      a = (sin(spline_x * pi/180) * shifted_spline_y), 
      b = (cos(spline_x * pi/180) * shifted_spline_y),
      c = shifted_spline_y
    ) %>% 
    summarise_all(sum) %>% 
    transmute(si = sqrt(a^2 + b^2)/c)
  
  half_range_y <- y_range/2
  
  y50 <- max_y$y - half_range_y
  
  desc_lower_curve_of_interest <- lower_curve_of_interest %>% 
    arrange(desc(x))
  
  desc_low_curve_y <- desc_lower_curve_of_interest$y - y50
  lower_x_index <- desc_low_curve_y %>% 
    detect_index(function (x) x < 0, .dir = "forward")
  
  upper_curve_y <- upper_curve_of_interest$y - y50
  upper_x_index <- upper_curve_y %>% 
    detect_index(function(x) x < 0, .dir = "forward")
  
  lower_x <- desc_lower_curve_of_interest[lower_x_index,]
  upper_x <- upper_curve_of_interest[upper_x_index,]
  
  if (dim(lower_x)[1] == 0
  ) {
    lower_x <- tibble(x = NA, y = NA)
  }
  
  if (dim(upper_x)[1] == 0
  ) {
    upper_x <- tibble(x = NA, y = NA)
  }
  
  ### preferred direction as the center of the spline
  prefdir_desc_lower_curve_of_interest <- 
    prefdir_lower_curve_of_interest %>% 
    arrange(desc(x))
  
  prefdir_desc_low_curve_y <- 
    prefdir_desc_lower_curve_of_interest$y - y50
  prefdir_lower_x_index <- 
    prefdir_desc_low_curve_y %>% 
    detect_index(function (x) x < 0, .dir = "forward")
  
  prefdir_upper_curve_y <- 
    prefdir_upper_curve_of_interest$y - y50
  prefdir_upper_x_index <- 
    prefdir_upper_curve_y %>% 
    detect_index(function(x) x < 0, .dir = "forward")
  
  prefdir_lower_x <- 
    prefdir_desc_lower_curve_of_interest[prefdir_lower_x_index,]
  prefdir_upper_x <- 
    prefdir_upper_curve_of_interest[prefdir_upper_x_index,]
  
  if (dim(prefdir_lower_x)[1] == 0
  ) {
    prefdir_lower_x <- tibble(x = NA, y = NA)
  }
  
  if (dim(prefdir_upper_x)[1] == 0
  ) {
    prefdir_upper_x <- tibble(x = NA, y = NA)
  }
  
  ## shift firing frame so only positive values
  firing_frame_shift <- test_cell_raw %>%
    mutate(
      firing = case_when(min(firing) < 0 ~ firing + abs(min(firing)), TRUE ~ firing))
  
  circ_data <- firing_frame_shift %>%
    pmap(
      function(direction, firing, ...) { 
        # ... captures unused columns from firing_frame
        rep(direction, round(firing * 10)) 
        # multiply by 10 to increase resolution
      }
    ) %>%
    unlist() %>%
    circular(units = "degrees", rotation = "clock")
  
  vm <- circ_data %>%
    mle.vonmises %>%
    unlist %>%
    `[`(c('mu','kappa')) %>%
    as_tibble %>%
    mutate(mu = mu %% 360)
  
  ## rows are mismatched because we translated the x-axis, 
  ## binding them together for plotting convenience
  dir_spline_fits_lm_ag[[i]] <- 
    bind_cols(predicted_curve_center, translated_x, 
              firing_ratio_y, prefdir_translated_x)  %>% 
    mutate(species = calc_tw_lm_ag$species[i], 
           file_path = calc_tw_lm_ag$file_path[i],
           filename = calc_tw_lm_ag$filename[i],
           prefdir_sum = calc_tw_lm_ag$prefdir_sum[i],
           mu = vm$mu
    )
  
  dir_width <- bind_cols(lower_x, upper_x, 
                         min_y, max_y, 
                         si_value, 
                         prefdir_lower_x, prefdir_upper_x) %>% 
    #select(x, x1) %>% 
    rename(half_pwr_lower_x = x...1, 
           half_pwr_lower_y = y...2, 
           half_pwr_upper_x = x...3, 
           half_pwr_upper_y = y...4, 
           curve_trough_x = x...5, 
           curve_trough_y = y...6, 
           curve_peak_x = x...7, 
           curve_peak_y = y...8,
           prefdir_half_pwr_lower_x = x...10,
           prefdir_half_pwr_lower_y = y...11,
           prefdir_half_pwr_upper_x = x...12,
           prefdir_half_pwr_upper_y = y...13
    ) %>% 
    mutate(
      width = half_pwr_upper_x - half_pwr_lower_x,
      prefdir_width = prefdir_half_pwr_upper_x - prefdir_half_pwr_lower_x) %>% 
    bind_cols(vm)
  
  dir_width_lm_ag <- bind_rows(dir_width_lm_ag, dir_width)
}

## combine with metadata
all_lm_dir_ag <- 
  bind_cols(calc_tw_lm_ag, dir_width_lm_ag) %>% 
  rename(si = si...9, 
         mu = mu...10, 
         kappa = kappa...11,
         si_full = si...29,
         mu_full = mu...36,
         kappa_full = kappa...37) %>%
  rowid_to_column(var = "id") %>% 
  mutate(investigator = "ag")

#### __AG data curation ####
## LM cells to be excluded
## These cells were determined to have at least one of the following issues:
## a) firing rates that were too low, b) poor fits and/or high noise:signal
exclude_all_ag <- c(182, 25, 128, 144, 155, 193, 196, 191, 6, 38, 29, 28, 27, 
                    26, 24, 54, 51, 42, 68, 89, 85, 84, 83,82, 81, 106,  
                    179, 176, 173, 172, 171, 170, 199, 197, 194, 189, 187, 
                    185, 181, 213, 209, 208, 204)

## exclude for tuning width (see Fig 3)  
exclude_width_only_ag <- c(147, 50, 168, 183, 215, 11, 34, 46, 87, 107, 104, 
                           135, 169, 210)
## combine with previous list
exclude_width_ag <- c(exclude_all_ag, exclude_width_only_ag)

## 2020-07-02 - what if we exclude all prolematic cells, period? 
## Curated data frame of cells for determining preferred direction
all_lm_dir_ag_curated <- all_lm_dir_ag %>% 
  filter(! id %in% exclude_width_ag)

## ORIGINAL BLOCK, IN CASE OF NEED TO REVERT:
# ## Curated data frame of cells for determining preferred direction
# all_lm_dir_ag_curated <- all_lm_dir_ag %>% 
#   filter(! id %in% exclude_all_ag)

## Curated data frame of cells for determining direction tuning width
all_lm_dir_width_ag_curated <- all_lm_dir_ag %>% 
  filter(! id %in% exclude_width_ag)

## And corresponding spline fits
spline_fits_lm_dir_width_ag_curated <- 
  dir_spline_fits_lm_ag[-exclude_width_ag]


## Code to plot any given cell:
plot_LM_dir_tuning <- function(i) { # i is the desired case number
  one_cell <-  calc_tw_lm_ag$file_path[i] 
  prefdir_sum <- calc_tw_lm_ag$prefdir_sum[i]
  
  baseline <- 
    suppressMessages(read_csv(one_cell, col_names = FALSE)) %>% 
    filter(X1 == 'NaN') %>% 
    select(-X1, -X2) %>% 
    pivot_longer(everything()) %>% 
    summarise(mean_baseline = mean(value)) %>% 
    pull()
  
  sem_baseline <- 
    suppressMessages(read_csv(one_cell, col_names = FALSE)) %>% 
    filter(X1 == 'NaN') %>% 
    select(-X1, -X2) %>% 
    pivot_longer(everything()) %>% 
    summarise(sem_baseline = sem(value)) %>% 
    pull()
  
  test_cell_mean <- 
    suppressMessages(read_csv(one_cell, 
                              col_names = FALSE, 
                              col_types = cols(.default = col_double()))) %>% 
    rename(direction = X1, speed = X2) %>% 
    filter(direction != 'NaN') %>% 
    mutate(stim_per = row_number()) %>% 
    pivot_longer(contains("X")) %>%  ## this selects all the 250ms columns
    group_by(direction, speed, stim_per) %>% 
    summarise(firing = mean(value - baseline)) %>% ## subtract baseline 
    ungroup() %>% 
    # need equal number of replicates for each direction. The modulus gives us
    # the number of rows (x) more than a complete set of directions. 
    # Exclude the last x trials, so that each direction is equally weighted
    filter(
      stim_per <= max(stim_per) - (max(stim_per) %% n_distinct(direction))) %>%
    group_by(direction, speed) %>%
    summarise(firing = mean(firing)) %>%
    ungroup() 
  
  test_cell_raw <- 
    suppressMessages(read_csv(one_cell, 
                              col_names = FALSE, 
                              col_types = cols(.default = col_double()))) %>% 
    rename(direction = X1, speed = X2) %>% 
    filter(direction != 'NaN') %>% 
    mutate(stim_per = row_number()) %>% 
    pivot_longer(contains("X")) %>%  ## this selects all the 250ms columns
    group_by(direction, speed, stim_per) %>% 
    summarise(firing = mean(value - baseline)) %>% ## subtract baseline 
    ungroup() %>% 
    # need equal number of replicates for each direction. The modulus gives us
    # the number of rows (x) more than a complete set of directions. 
    # Exclude the last x trials, so that each direction is equally weighted
    filter(stim_per <= max(stim_per) - (max(stim_per) %% n_distinct(direction)))
  
  ## Adding two additional loops to ensure wrapping between 0 and 360
  x_raw <- c(test_cell_raw$direction-360, 
             test_cell_raw$direction, test_cell_raw$direction + 360)
  y_raw <- c(test_cell_raw$firing, test_cell_raw$firing, test_cell_raw$firing)
  
  plot(x_raw, y_raw)
  abline(v = 0) ## lower bound of tested domain
  abline(v = 315) ## upper bound of tested domain
  points(dir_spline_fits_lm_ag[[i]]$spline_x,
         dir_spline_fits_lm_ag[[i]]$spline_y,
         col = "purple", cex = 0.5, pch = 19)
}  

## Examples:
# plot_LM_dir_tuning(1)   # good cell
# plot_LM_dir_tuning(50)


## Quick plot of AG direction histogram
# all_lm_dir_ag_curated %>%
#   filter(species == 'hb') %>%
#   ggplot(aes(x = prefdir_sum)) +
#   geom_histogram(binwidth = 15, boundary = 0,
#                  fill = col_hb,
#                  color = "black",
#                  size = 0.15) +
#   coord_polar(direction = 1, start = pi/2) +
#   scale_y_continuous(trans = "sqrt") +
#   scale_x_continuous(breaks = c(0, 90, 180, 270),
#                      expand = c(0,0),
#                      limits = c(0,360)) +
#   theme_bw() +
#   theme(
#     #legend.key = element_blank(),
#     axis.text.x = element_blank(),
#     axis.text.y = element_blank(),
#     axis.title = element_blank(),
#     panel.border = element_blank(),
#     axis.ticks.y = element_blank()
#   ) +
#   geom_text(x = 0, y = sqrt(2), label = "2", size = 3) +
#   geom_text(x = 0, y = sqrt(4), label = "4", size = 3) +
#   geom_text(x = 0, y = sqrt(12), label = "F") +
#   geom_text(x = 270, y = sqrt(12), label = "U") +
#   geom_text(x = 180, y = sqrt(12), label = "B") +
#   geom_text(x = 90, y = sqrt(12), label = "D")


#### __GS data ####
hb_file_list_LM_gs <- 
  list.files(path = "./data/direction/gs_LM_dir_hb/",
             full.names = TRUE) 

zf_file_list_LM_gs <- 
  list.files(path = "./data/direction/gs_LM_dir_zf/",
             full.names = TRUE) 

hb_prefdir_gs <- map_dfr(hb_file_list_LM_gs, gs_get_prefdir, species = "hb")
zf_prefdir_gs <- map_dfr(zf_file_list_LM_gs, gs_get_prefdir, species = "zf")

prefdir_df_lm_gs <- bind_rows(zf_prefdir_gs, hb_prefdir_gs) 

## add file path to data frame
filepath_hb_lm_gs <- 
  list.files(path = "./data/direction/gs_LM_dir_hb/",
             full.names = F) %>% 
  cbind(file_path = list.files(
    path = "./data/direction/gs_LM_dir_hb/",
    full.names = T))

filepath_zf_lm_gs <- 
  list.files(path = "./data/direction/gs_LM_dir_zf/",
             full.names = F) %>% 
  cbind(file_path = list.files(
    path = "./data/direction/gs_LM_dir_zf/",
    full.names = T))

filepath_all_lm_gs <- rbind(filepath_zf_lm_gs, filepath_hb_lm_gs)
filepath_all_lm_gs <- 
  as_tibble(filepath_all_lm_gs) %>% 
  rename(filename = ".")

calc_tw_lm_gs <- left_join(prefdir_df_lm_gs, 
                           filepath_all_lm_gs, by = "filename")

## calculate widths of direction tuning curves
GS_dir_width_lm <- NULL
GS_dir_spline_fits_lm <- NULL
for (i in 1:dim(calc_tw_lm_gs)[1]) {
  one_cell <-  calc_tw_lm_gs$file_path[i] 
  prefdir_sum <- calc_tw_lm_gs$prefdir_sum[i]
  
  baseline <- 
    suppressMessages(read_csv(one_cell, col_names = FALSE)) %>% 
    filter(X1 == 'NaN') %>% 
    select(-X1, -X2) %>% 
    pivot_longer(everything()) %>% 
    summarise(mean_baseline = mean(value)) %>% 
    pull()
  
  sem_baseline <- 
    suppressMessages(read_csv(one_cell, col_names = FALSE)) %>% 
    filter(X1 == 'NaN') %>% 
    select(-X1, -X2) %>% 
    pivot_longer(everything()) %>% 
    summarise(sem_baseline = sem(value)) %>% 
    pull()
  
  test_cell_mean <- 
    suppressMessages(read_csv(one_cell, 
                              col_names = FALSE, 
                              col_types = cols(.default = col_double()))) %>% 
    rename(direction = X1, speed = X2) %>% 
    filter(direction != 'NaN') %>% 
    mutate(stim_per = row_number()) %>% 
    pivot_longer(contains("X")) %>%  ## this selects all the 250ms columns
    group_by(direction, speed, stim_per) %>% 
    summarise(firing = mean(value - baseline)) %>% ## subtract baseline 
    ungroup() %>% 
    # need equal number of replicates for each direction. The modulus gives us
    # the number of rows (x) more than a complete set of directions. 
    # Exclude the last x trials, so that each direction is equally weighted
    filter(
      stim_per <= max(stim_per) - (max(stim_per) %% n_distinct(direction))) %>%
    group_by(direction, speed) %>%
    summarise(firing = mean(firing)) %>%
    ungroup() 
  
  test_cell_raw <- 
    suppressMessages(read_csv(one_cell, 
                              col_names = FALSE, 
                              col_types = cols(.default = col_double()))) %>% 
    rename(direction = X1, speed = X2) %>% 
    filter(direction != 'NaN') %>% 
    mutate(stim_per = row_number()) %>% 
    pivot_longer(contains("X")) %>%  ## this selects all the 250ms columns
    group_by(direction, speed, stim_per) %>% 
    summarise(firing = mean(value - baseline)) %>% ## subtract baseline 
    ungroup() %>% 
    # need equal number of replicates for each direction. The modulus gives us
    # the number of rows (x) more than a complete set of directions. 
    # Exclude the last x trials, so that each direction is equally weighted
    filter(
      stim_per <= max(stim_per) - (max(stim_per) %% n_distinct(direction)))
  
  x_raw <- c(test_cell_raw$direction - 360,
             test_cell_raw$direction, 
             test_cell_raw$direction + 360)
  y_raw <- c(test_cell_raw$firing, test_cell_raw$firing, test_cell_raw$firing)
  
  fit5 <- smooth.spline(x = x_raw, y = y_raw, df = 5)
  fit6 <- smooth.spline(x = x_raw, y = y_raw, df = 6)
  fit7 <- smooth.spline(x = x_raw, y = y_raw, df = 7)
  fit8 <- smooth.spline(x = x_raw, y = y_raw, df = 8)
  fit9 <- smooth.spline(x = x_raw, y = y_raw, df = 9)
  fit10 <- smooth.spline(x = x_raw, y = y_raw, df = 10)
  fit11 <- smooth.spline(x = x_raw, y = y_raw, df = 11)
  fit12 <- smooth.spline(x = x_raw, y = y_raw, df = 12)
  fit13 <- smooth.spline(x = x_raw, y = y_raw, df = 13)
  fit14 <- smooth.spline(x = x_raw, y = y_raw, df = 14)
  fit15 <- smooth.spline(x = x_raw, y = y_raw, df = 15)
  fit16 <- smooth.spline(x = x_raw, y = y_raw, df = 16)
  fit17 <- smooth.spline(x = x_raw, y = y_raw, df = 17)
  fit18 <- smooth.spline(x = x_raw, y = y_raw, df = 18)
  fit19 <- smooth.spline(x = x_raw, y = y_raw, df = 19)
  fit20 <- smooth.spline(x = x_raw, y = y_raw, df = 20)
  
  fit_list <- list("5" = fit5, "6" = fit6, "7" = fit7, "8" = fit8, 
                   "9" = fit9, "10" = fit10, "11" = fit11, "12" = fit12,
                   "13" = fit13, "14" = fit14, "15" = fit15, "16" = fit16,
                   "17" = fit17, "18" = fit18, "19" = fit19, "20" = fit20)
  
  best_fit <- map_dfr(fit_list, AICc_smooth_spline) %>% 
    pivot_longer(cols = everything(), 
                 names_to = "fit", 
                 values_to = "AICc") %>% 
    #mutate(fit = as.numeric(fit)) %>% 
    filter(AICc == min(AICc))
  
  best_fit_contents <- fit_list[best_fit$fit][[1]]
  
  y_vals <- stats:::predict.smooth.spline(best_fit_contents,
                                          x = seq(0, 359, by = 1))
  predicted_values <- tibble(x = y_vals$x, y = y_vals$y)
  
  predicted_raw_vals <- stats:::predict.smooth.spline(best_fit_contents, 
                                                      x = seq(-360, 720, 1))
  predicted_raw_df <- tibble(x = predicted_raw_vals$x, 
                             y = predicted_raw_vals$y)
  
  max_y <- predicted_raw_df %>% 
    filter(x >= 0 & x < 360) %>% 
    filter(y == max(y))
  
  min_y <- predicted_raw_df %>% 
    filter(x >= 0 & x < 360) %>% 
    filter(y == min(y))
  
  y_range <- max_y$y - min_y$y
  
  lower_curve_of_interest <- predicted_raw_df %>% 
    filter(x >= max_y$x - 180 & x<= max_y$x)
  
  upper_curve_of_interest <- predicted_raw_df %>% 
    filter(x >= max_y$x & x<= max_y$x + 180)
  
  predicted_curve_center <- predicted_raw_df %>% 
    filter(x >= 0 & x < 360) %>% 
    rename(spline_x = x, spline_y = y) %>% 
    mutate(
      shifted_spline_y = spline_y- min(spline_y), 
      standardized_x = spline_x/360, 
  standardized_y = (spline_y - min(spline_y))/(max(spline_y) - min(spline_y)))
  
  translated_x <- predicted_raw_df %>% 
    filter(x >= max_y$x - 180 & x < max_y$x + 180) %>% 
    rename(trans_spline_x = x, trans_spline_y = y) %>% 
    mutate(
      trans_std_x = trans_spline_x - max_y$x, 
      trans_std_y = (trans_spline_y - min(trans_spline_y))/(max(trans_spline_y) - min(trans_spline_y)))
  
  firing_ratio_y <- predicted_raw_df  %>% 
    filter(x >= max_y$x - 180 & x < max_y$x + 180) %>% 
    mutate(bsln_ratio_y = (y + baseline)/baseline) %>% 
    select(bsln_ratio_y)
  
  ### using prefdir_sum instead of spline peak
  prefdir_lower_curve_of_interest <- predicted_raw_df %>% 
    filter(x >= prefdir_sum - 180 & x<= prefdir_sum)
  
  prefdir_upper_curve_of_interest <- predicted_raw_df %>% 
    filter(x >= prefdir_sum & x<= prefdir_sum + 180)
  
  prefdir_translated_x <- predicted_raw_df %>% 
    filter(x >= prefdir_sum - 180 & x <= prefdir_sum + 180) %>% 
    rename(prefdir_trans_spline_x = x, prefdir_trans_spline_y = y) %>% 
    mutate(
      prefdir_trans_std_x = prefdir_trans_spline_x - prefdir_sum, 
      prefdir_trans_std_y = (prefdir_trans_spline_y - min(prefdir_trans_spline_y))/(max(prefdir_trans_spline_y) - min(prefdir_trans_spline_y)))
  
  si_value <- predicted_curve_center %>% 
    transmute(
      a = (sin(spline_x * pi/180) * shifted_spline_y), 
      b = (cos(spline_x * pi/180) * shifted_spline_y),
      c = shifted_spline_y
    ) %>% 
    summarise_all(sum) %>% 
    transmute(si = sqrt(a^2 + b^2)/c)
  
  half_range_y <- y_range/2
  
  y50 <- max_y$y - half_range_y
  
  desc_lower_curve_of_interest <- lower_curve_of_interest %>% 
    arrange(desc(x))
  
  desc_low_curve_y <- desc_lower_curve_of_interest$y - y50
  lower_x_index <- desc_low_curve_y %>% 
    detect_index(function (x) x < 0, .dir = "forward")
  
  upper_curve_y <- upper_curve_of_interest$y - y50
  upper_x_index <- upper_curve_y %>% 
    detect_index(function(x) x < 0, .dir = "forward")
  
  lower_x <- desc_lower_curve_of_interest[lower_x_index,]
  upper_x <- upper_curve_of_interest[upper_x_index,]
  
  if (dim(lower_x)[1] == 0
  ) {
    lower_x <- tibble(x = NA, y = NA)
  }
  
  if (dim(upper_x)[1] == 0
  ) {
    upper_x <- tibble(x = NA, y = NA)
  }
  
  ### preferred direction as the center of the spline
  prefdir_desc_lower_curve_of_interest <- 
    prefdir_lower_curve_of_interest %>% 
    arrange(desc(x))
  
  prefdir_desc_low_curve_y <- 
    prefdir_desc_lower_curve_of_interest$y - y50
  prefdir_lower_x_index <- 
    prefdir_desc_low_curve_y %>% 
    detect_index(function (x) x < 0, .dir = "forward")
  
  prefdir_upper_curve_y <- 
    prefdir_upper_curve_of_interest$y - y50
  prefdir_upper_x_index <- 
    prefdir_upper_curve_y %>% 
    detect_index(function(x) x < 0, .dir = "forward")
  
  prefdir_lower_x <- 
    prefdir_desc_lower_curve_of_interest[prefdir_lower_x_index,]
  prefdir_upper_x <- 
    prefdir_upper_curve_of_interest[prefdir_upper_x_index,]
  
  if (dim(prefdir_lower_x)[1] == 0
  ) {
    prefdir_lower_x <- tibble(x = NA, y = NA)
  }
  
  if (dim(prefdir_upper_x)[1] == 0
  ) {
    prefdir_upper_x <- tibble(x = NA, y = NA)
  }
  
  ## shift firing frame so only positive values
  firing_frame_shift <- test_cell_raw %>%
    mutate(
      firing = 
        case_when(min(firing) < 0 ~ firing + abs(min(firing)), TRUE ~ firing))
  
  circ_data <- firing_frame_shift %>%
    pmap(
      function(direction, firing, ...) { 
        # ... captures unused columns from firing_frame
        rep(direction, round(firing * 10)) 
        # multiply by 10 to increase resolution
      }
    ) %>%
    unlist() %>%
    circular(units = "degrees", rotation = "clock")
  
  vm <- circ_data %>%
    mle.vonmises %>%
    unlist %>%
    `[`(c('mu','kappa')) %>%
    as_tibble %>%
    mutate(mu = mu %% 360)
  
  ## rows are mismatched because we translated the x-axis, 
  ## binding them together for plotting convenience
  GS_dir_spline_fits_lm[[i]] <- 
    bind_cols(predicted_curve_center, 
              translated_x, 
              firing_ratio_y, 
              prefdir_translated_x)  %>% 
    mutate(species = calc_tw_lm_gs$species[i], 
           file_path = calc_tw_lm_gs$file_path[i],
           filename = calc_tw_lm_gs$filename[i],
           prefdir_sum = calc_tw_lm_gs$prefdir_sum[i],
           mu = vm$mu
    )
  
  dir_width <- bind_cols(lower_x, upper_x, 
                         min_y, max_y, 
                         si_value, 
                         prefdir_lower_x, prefdir_upper_x) %>% 
    #select(x, x1) %>% 
    rename(half_pwr_lower_x = x...1, 
           half_pwr_lower_y = y...2, 
           half_pwr_upper_x = x...3, 
           half_pwr_upper_y = y...4, 
           curve_trough_x = x...5, 
           curve_trough_y = y...6, 
           curve_peak_x = x...7, 
           curve_peak_y = y...8,
           prefdir_half_pwr_lower_x = x...10,
           prefdir_half_pwr_lower_y = y...11,
           prefdir_half_pwr_upper_x = x...12,
           prefdir_half_pwr_upper_y = y...13
    ) %>% 
    mutate(
      width = half_pwr_upper_x - half_pwr_lower_x,
      prefdir_width = prefdir_half_pwr_upper_x - prefdir_half_pwr_lower_x) %>% 
    bind_cols(vm)
  
  GS_dir_width_lm <- bind_rows(GS_dir_width_lm, dir_width)
  
}

## combine with metadata
all_lm_dir_gs <- 
  bind_cols(calc_tw_lm_gs, GS_dir_width_lm) %>%
  rename(mu = mu...10, 
         kappa = kappa...11,
         mu_full = mu...36,
         kappa_full = kappa...37) %>%
  rowid_to_column(var = "id") %>% 
  mutate(investigator = "gs")

#### __GS data curation ####
## LM cells to be excluded
GS_excluded <- c(35, 31, 24, 59, 49, 43, 79, 97, 96, 95, 87, 119, 108, 103, 157,
                 156, 148, 141, 174, 172, 171, 166, 160)

## Curated data frame of cells for determining preferred direction
all_lm_dir_gs_curated <- all_lm_dir_gs %>% 
  filter(! id %in% GS_excluded)

## Quick plot of GS direction histogram
# all_lm_dir_gs_curated %>%
#   filter(species == 'hb') %>% 
#   ggplot(aes(x = prefdir_sum)) +
#   geom_histogram(binwidth = 15, boundary = 0, 
#                  fill = col_hb, 
#                  color = "black",
#                  size = 0.15) +
#   coord_polar(direction = 1, start = pi/2) + 
#   scale_y_continuous(trans = "sqrt") +
#   scale_x_continuous(breaks = c(0, 90,180,270), 
#                      expand = c(0,0), limits = c(0,360)) +
#   theme_bw() +
#   theme(
#     #legend.key = element_blank(),
#     axis.text.x = element_blank(),
#     axis.text.y = element_blank(),
#     axis.title = element_blank(),
#     panel.border = element_blank(),
#     axis.ticks.y = element_blank()
#   ) +
#   geom_text(x = 0, y = sqrt(2), label = "2", size = 3) +
#   geom_text(x = 0, y = sqrt(4), label = "4", size = 3) +
#   geom_text(x = 0, y = sqrt(12), label = "F") +
#   geom_text(x = 270, y = sqrt(12), label = "U") +
#   geom_text(x = 180, y = sqrt(12), label = "B") +
#   geom_text(x = 90, y = sqrt(12), label = "D")

## Curated data frame of cells for determining direction tuning width
all_lm_dir_width_gs_curated <- all_lm_dir_gs_curated  

## And corresponding spline fits
spline_fits_lm_dir_width_gs_curated <- GS_dir_spline_fits_lm[-GS_excluded]

#### Combine AG and GS #####
all_lm_dir_width_tmp <- 
  bind_rows(all_lm_dir_ag_curated,
            all_lm_dir_gs_curated) %>%
  rename(investigator_id = id) %>% 
  ## add ID number now that it's all combined
  rowid_to_column(var = "combined_id") %>%
  ## add column of differentials between width computation methods
  mutate(width_differential = abs(width - prefdir_width)) 

## Remove problematic cells from the polar histogram data frame
all_lm_dir_width <- 
  all_lm_dir_width_tmp %>%
  filter(!width_differential > 0)

all_lm_tuning_width_tmp <-
  bind_rows(all_lm_dir_width_ag_curated,
            all_lm_dir_width_gs_curated) %>%
  ## add ID number from _ag or _gs data frames
  rename(investigator_id = id) %>% 
  ## add ID number now that it's all combined
  rowid_to_column(var = "combined_id") %>%
  ## add column of differentials between width computation methods
  mutate(width_differential = abs(width - prefdir_width)) 

## Identify the problematic cells
mismatch_cells <- all_lm_tuning_width_tmp %>%
  filter(width_differential > 1 | is.na(width_differential))

## Remove the problematic cells from the tuning width data frame
all_lm_tuning_width <- 
  all_lm_tuning_width_tmp %>%
  filter(!width_differential > 0)

all_lm_dir_spline_fits_list_tmp <-
  c(spline_fits_lm_dir_width_ag_curated,
    spline_fits_lm_dir_width_gs_curated)

all_lm_dir_spline_fits_list <-
  all_lm_dir_spline_fits_list_tmp[-mismatch_cells$combined_id]

all_lm_dir_spline_fits_df <-
  bind_rows(all_lm_dir_spline_fits_list, .id = "combined_id")


#### Pigeon LM direction tuning data ####
pg_lm_dir_raw <- 
  "./data/direction/Pigeon_LM_data_2020/recovering_pg_LM_data_april_2020_AHG.xlsx"

## Change column names, put the preferred direction into AG's reference frame,
## and ensure it is on a 0-360 scale
pg_lm_dir_metadata <- 
  read_excel(pg_lm_dir_raw, sheet = 1) %>% 
  select(cell, type, file_name, screen_corr, 
         old_cell_id, prefdir_sum_AHG_frame) %>%
  mutate(species = "pg")

## Read in the pigeon cells: each sheet is an individual cell
sheets <- excel_sheets(pg_lm_dir_raw)
pg_dir_df_lm <- NULL

for (i in 2:length(sheets)) {
  cell <- read_excel(pg_lm_dir_raw, sheet = sheets[i])
  pg_dir_df_lm <- bind_rows(pg_dir_df_lm, cell)
}

## read in the csv files and get mean response at each direction
pg_dir_csv_df_lm <- NULL
for (i in 1:dim(pg_lm_dir_metadata %>% 
                filter(!is.na(file_name)))[1]) {
  one_cell <- 
    paste0("./data/direction/Pigeon_LM_data_2020/pigeon_psth_confirmed_DRW_2020/", pg_lm_dir_metadata$file_name[i])
  
  baseline <- 
    suppressMessages(read_csv(one_cell, 
                              col_names = FALSE)) %>% 
    filter(X1 == 'NaN') %>% 
    select(-X1, -X2) %>% 
    pivot_longer(everything()) %>% 
    summarise(mean_baseline = mean(value)) %>% 
    pull()
  
  sem_baseline <- 
    suppressMessages(read_csv(one_cell, 
                              col_names = FALSE)) %>% 
    filter(X1 == 'NaN') %>% 
    select(-X1, -X2) %>% 
    pivot_longer(everything()) %>% 
    summarise(sem_baseline = sem(value)) %>% 
    pull()
  
  test_cell_mean <- 
    suppressMessages(read_csv(one_cell, 
                              col_names = FALSE, 
                              col_types = cols(.default = col_double()))) %>% 
    rename(direction_DRW_frame = X1, speed = X2) %>% 
    filter(direction_DRW_frame != 'NaN') %>% 
    mutate(stim_per = row_number()) %>% 
    pivot_longer(contains("X")) %>%  ## this selects all the 250ms columns
    group_by(direction_DRW_frame, speed, stim_per) %>% 
    summarise(firing = mean(value - baseline)) %>% ## subtract baseline 
    ungroup() %>% 
    # need equal number of replicates for each direction. The modulus gives us
    # the number of rows (x) more than a complete set of directions. 
    # Exclude the last x trials, so that each direction is equally weighted
    filter(
      stim_per <= max(stim_per) - (max(stim_per) %% n_distinct(direction_DRW_frame))) %>%
    group_by(direction_DRW_frame, speed) %>%
    summarise(firing = mean(firing)) %>%
    ungroup() %>% 
    mutate(
      direction_corr = direction_DRW_frame - pg_lm_dir_metadata$screen_corr[i],
      direction_AHG_frame = (direction_corr-90)%%360,
      cell = pg_lm_dir_metadata$cell[i],
      spike_less_bsln = firing,
      spontaneous_rate = baseline
    )
  
  pg_dir_csv_df_lm[[i]] <- test_cell_mean
  
}

pg_dir_csv_lm_collapsed <- bind_rows(pg_dir_csv_df_lm)

pg_key_csv_info_lm <- 
  pg_dir_csv_lm_collapsed %>% 
  select(cell, direction_AHG_frame, spike_less_bsln, spontaneous_rate) %>% 
  left_join(pg_lm_dir_metadata, by = "cell")

pg_dir_df_lm_final <- 
  pg_dir_df_lm %>% 
  mutate(spike_less_bsln = mean_spike_rate - spontaneous_rate) %>% 
  left_join(pg_lm_dir_metadata, by = "cell") %>% 
  mutate(direction_corr = (direction_DRW_frame - screen_corr)%%360,
         direction_AHG_frame = (direction_corr-90)%%360)

pg_lm_dir_all_data <- bind_rows(pg_key_csv_info_lm, pg_dir_df_lm_final)
pg_lm_dir_cell_names <- unique(pg_lm_dir_all_data$cell)

pg_dir_width_lm <- NULL
pg_spline_fits_lm_dir <- NULL
for (i in 1:length(pg_lm_dir_cell_names)) {
  one_cell <- pg_lm_dir_all_data %>% 
    filter(cell == pg_lm_dir_cell_names[i])
  
  baseline <- mean(one_cell$spontaneous_rate)
  prefdir_sum <- mean(one_cell$prefdir_sum_AHG_frame)
  
  x_mean <- c(one_cell$direction_AHG_frame-360, 
              one_cell$direction_AHG_frame, 
              one_cell$direction_AHG_frame + 360)
  
  y_mean <- c(one_cell$spike_less_bsln, 
              one_cell$spike_less_bsln,
              one_cell$spike_less_bsln)
  
  fit5 <- smooth.spline(x = x_mean, y = y_mean, df = 5)
  fit6 <- smooth.spline(x = x_mean, y = y_mean, df = 6)
  fit7 <- smooth.spline(x = x_mean, y = y_mean, df = 7)
  fit8 <- smooth.spline(x = x_mean, y = y_mean, df = 8)
  fit9 <- smooth.spline(x = x_mean, y = y_mean, df = 9)
  fit10 <- smooth.spline(x = x_mean, y = y_mean, df = 10)
  fit11 <- smooth.spline(x = x_mean, y = y_mean, df = 11)
  fit12 <- smooth.spline(x = x_mean, y = y_mean, df = 12)
  fit13 <- smooth.spline(x = x_mean, y = y_mean, df = 13)
  fit14 <- smooth.spline(x = x_mean, y = y_mean, df = 14)
  fit15 <- smooth.spline(x = x_mean, y = y_mean, df = 15)
  fit16 <- smooth.spline(x = x_mean, y = y_mean, df = 16)
  fit17 <- smooth.spline(x = x_mean, y = y_mean, df = 17)
  fit18 <- smooth.spline(x = x_mean, y = y_mean, df = 18)
  fit19 <- smooth.spline(x = x_mean, y = y_mean, df = 19)
  fit20 <- smooth.spline(x = x_mean, y = y_mean, df = 20)
  
  fit_list <- list("5" = fit5, "6" = fit6, "7" = fit7, 
                   "8" = fit8, "9" = fit9, "10" = fit10, 
                   "11" = fit11, "12" = fit12, "13" = fit13, 
                   "14" = fit14, "15" = fit15, "16" = fit16, 
                   "17" = fit17, "18" = fit18, "19" = fit19, 
                   "20" = fit20)
  
  best_fit <- map_dfr(fit_list, AICc_smooth_spline) %>% 
    pivot_longer(cols = everything(), 
                 names_to = "fit", 
                 values_to = "AICc") %>% 
    filter(AICc == min(AICc))
  
  best_fit_contents <- fit_list[best_fit$fit][[1]]
  
  y_vals <- stats:::predict.smooth.spline(best_fit_contents, 
                                          x = seq(0, 359, by = 1))
  predicted_values <- tibble(x = y_vals$x, y = y_vals$y)
  
  predicted_raw_vals <- stats:::predict.smooth.spline(best_fit_contents, 
                                                      x = seq(-360, 720, 1))
  predicted_raw_df <- tibble(x = predicted_raw_vals$x, 
                             y = predicted_raw_vals$y)
  
  max_y <- predicted_raw_df %>% 
    filter(x >= 0 & x< 360) %>% 
    filter(y == max(y))
  
  min_y <- predicted_raw_df %>% 
    filter(x >= 0 & x< 360) %>% 
    filter(y == min(y))
  
  y_range <- max_y$y - min_y$y
  
  lower_curve_of_interest <- predicted_raw_df %>% 
    filter(x >= max_y$x - 180 & x<= max_y$x)
  
  upper_curve_of_interest <- predicted_raw_df %>% 
    filter(x >= max_y$x & x<= max_y$x + 180)
  
  predicted_curve_center <- predicted_raw_df %>% 
    filter(x >= 0 & x < 360) %>% 
    rename(spline_x = x, spline_y = y) %>% 
    mutate(shifted_spline_y = spline_y- min(spline_y), 
           standardized_x = spline_x/360, 
           standardized_y = 
             (spline_y - min(spline_y))/(max(spline_y) - min(spline_y)))
  
  translated_x <- predicted_raw_df %>% 
    filter(x >= max_y$x - 180 & x < max_y$x + 180) %>% 
    rename(trans_spline_x = x, trans_spline_y = y) %>% 
    mutate(trans_std_x = trans_spline_x - max_y$x, 
           trans_std_y = 
             (trans_spline_y - min(trans_spline_y))/(max(trans_spline_y) - min(trans_spline_y)))
  
  firing_ratio_y <- predicted_raw_df  %>% 
    filter(x >= max_y$x - 180 & x < max_y$x + 180) %>% 
    mutate(bsln_ratio_y = (y + baseline)/baseline) %>% 
    select(bsln_ratio_y)
  
  ### If using prefdir_sum instead of spline peak
  prefdir_lower_curve_of_interest <- predicted_raw_df %>% 
    filter(x >= prefdir_sum - 180 & x<= prefdir_sum)
  
  prefdir_upper_curve_of_interest <- predicted_raw_df %>% 
    filter(x >= prefdir_sum & x<= prefdir_sum + 180)

  prefdir_translated_x <- predicted_raw_df %>% 
    filter(x >= prefdir_sum - 180 & x <= prefdir_sum + 180) %>% 
    rename(prefdir_trans_spline_x = x, prefdir_trans_spline_y = y) %>% 
    mutate(prefdir_trans_std_x = prefdir_trans_spline_x - prefdir_sum, 
           prefdir_trans_std_y = 
             (prefdir_trans_spline_y - min(prefdir_trans_spline_y))/(max(prefdir_trans_spline_y) - min(prefdir_trans_spline_y)))
  
  si_value <- predicted_curve_center %>% 
    transmute(
      a = (sin(spline_x * pi/180) * shifted_spline_y), 
      b = (cos(spline_x * pi/180) * shifted_spline_y),
      c = shifted_spline_y
    ) %>% 
    summarise_all(sum) %>% 
    transmute(si = sqrt(a^2 + b^2)/c)
  
  half_range_y <- y_range/2
  
  y50 <- max_y$y - half_range_y
  
  desc_lower_curve_of_interest <- lower_curve_of_interest %>% 
    arrange(desc(x))
  
  desc_low_curve_y <- desc_lower_curve_of_interest$y - y50
  lower_x_index <- desc_low_curve_y %>% 
    detect_index(function (x) x < 0, .dir = "forward")
  
  upper_curve_y <- upper_curve_of_interest$y - y50
  upper_x_index <- upper_curve_y %>% 
    detect_index(function(x) x < 0, .dir = "forward")
  
  lower_x <- desc_lower_curve_of_interest[lower_x_index,]
  upper_x <- upper_curve_of_interest[upper_x_index,]

  if (dim(lower_x)[1] == 0
  ) {
    lower_x <- tibble(x = NA, y = NA)
  }
  
  if (dim(upper_x)[1] == 0
  ) {
    upper_x <- tibble(x = NA, y = NA)
  }

  ### preferred direction as the center of the spline
  prefdir_desc_lower_curve_of_interest <- 
    prefdir_lower_curve_of_interest %>% 
    arrange(desc(x))
  
  prefdir_desc_low_curve_y <- 
    prefdir_desc_lower_curve_of_interest$y - y50
  prefdir_lower_x_index <- prefdir_desc_low_curve_y %>% 
    detect_index(function (x) x < 0, .dir = "forward")
  
  prefdir_upper_curve_y <- prefdir_upper_curve_of_interest$y - y50
  prefdir_upper_x_index <- prefdir_upper_curve_y %>% 
    detect_index(function(x) x < 0, .dir = "forward")
  
  prefdir_lower_x <- 
    prefdir_desc_lower_curve_of_interest[prefdir_lower_x_index,]
  prefdir_upper_x <- 
    prefdir_upper_curve_of_interest[prefdir_upper_x_index,]
  
  if (dim(prefdir_lower_x)[1] == 0
  ) {
    prefdir_lower_x <- tibble(x = NA, y = NA)
  }
  
  if (dim(prefdir_upper_x)[1] == 0
  ) {
    prefdir_upper_x <- tibble(x = NA, y = NA)
  }
  
  ## shift firing frame so only positive values
  firing_frame_shift <- one_cell %>%
    mutate(
      firing = case_when(min(spike_less_bsln) < 0 ~ spike_less_bsln + abs(min(spike_less_bsln)), TRUE ~ spike_less_bsln))
  
  circ_data <- firing_frame_shift %>%
    pmap(
      function(direction_AHG_frame, firing, ...) { 
        # ... captures unused columns from firing_frame
        rep(direction_AHG_frame, round(firing * 10)) 
        # multiply by 10 to increase resolution
      }
    ) %>%
    unlist() %>%
    circular(units = "degrees", rotation = "clock")
  
  vm <- circ_data %>%
    mle.vonmises %>%
    unlist %>%
    `[`(c('mu','kappa')) %>%
    as_tibble %>%
    mutate(mu = mu %% 360)
  
  ## rows are mismatched because we translated the x-axis, 
  ## binding them together for plotting convenience
  pg_spline_fits_lm_dir[[i]] <- 
    bind_cols(predicted_curve_center, 
              translated_x, 
              firing_ratio_y, 
              prefdir_translated_x)  %>% 
    mutate(species = "pg", 
           file_path = pg_lm_dir_cell_names[i],
           filename = pg_lm_dir_cell_names[i],
           prefdir_sum = prefdir_sum,
           mu = vm$mu
    )
  
  dir_width <- bind_cols(cell_id = pg_lm_dir_cell_names[i],
                         species = "pg",
                         prefdir_sum = prefdir_sum, 
                         lower_x, upper_x, 
                         min_y, max_y, 
                         si_value, 
                         prefdir_lower_x, prefdir_upper_x) %>% 
    rename(half_pwr_lower_x = x...4, 
           half_pwr_lower_y = y...5, 
           half_pwr_upper_x = x...6, 
           half_pwr_upper_y = y...7, 
           curve_trough_x = x...8, 
           curve_trough_y = y...9, 
           curve_peak_x = x...10, 
           curve_peak_y = y...11,
           prefdir_half_pwr_lower_x = x...13,
           prefdir_half_pwr_lower_y = y...14,
           prefdir_half_pwr_upper_x = x...15,
           prefdir_half_pwr_upper_y = y...16
    ) %>% 
    mutate(
      width = half_pwr_upper_x - half_pwr_lower_x,
      prefdir_width = prefdir_half_pwr_upper_x - prefdir_half_pwr_lower_x) %>% 
    bind_cols(vm)
  
  pg_dir_width_lm <- bind_rows(pg_dir_width_lm, dir_width)
}

## Bind the rows of the list of the plotting objects
pg_spline_plotting_df_lm <- bind_rows(pg_spline_fits_lm_dir, .id = "cell")

#### Combined objects for all species #####
## Combine pg dir width metrics with metadata
pg_dir_width_lm_renamed <-
  pg_dir_width_lm %>%
  rename(cell = cell_id)

pg_LM_tuning_width_summary <- left_join(pg_lm_dir_metadata,
                                        pg_dir_width_lm_renamed,
                                        by = c("cell", "species"))

## Update splines lists and dfs
all_lm_dir_spline_fits_list_pg <-
  c(all_lm_dir_spline_fits_list,
    pg_spline_fits_lm_dir)

all_lm_dir_spline_fits_df <-
  bind_rows(all_lm_dir_spline_fits_list_pg, .id = "combined_id")


###### Separate LM data by species ########
hb_prefdir_lm_all <- 
  all_lm_dir_width %>%
  filter(species == "hb")

zf_prefdir_lm_all <- 
  all_lm_dir_width %>%
  filter(species == "zf")

pg_prefdir_lm <-
  pg_LM_tuning_width_summary %>%
  filter(species == "pg")

###### Panel B Hummingbird LM ########
panel_1B <- 
  all_lm_dir_width %>%
  filter(species == "hb") %>%
  #filter(investigator == "ag") %>%
  ggplot(aes(x = prefdir_sum)) +
  geom_hline(yintercept = c(0.10, 0.25), 
             colour = "grey90", size = 0.2) +
  geom_vline(xintercept = seq(0, 360-1, by = 90), 
             colour = "grey90", size = 0.2) +
  geom_histogram(aes(y = stat(width*density)), 
                 binwidth = 15, 
                 boundary = 0, 
                 fill = col_hb, 
                 color = "black",
                 size = 0.15) +
  coord_polar(direction = 1, start = pi/2) + 
  scale_y_continuous(trans = "sqrt", 
                     breaks = c(0.10, 0.25),
                     expand = c(0,0), limits = c(0, 0.42)) +
  scale_x_continuous(breaks = c(0, 90,180,270), expand = c(0,0), 
                     limits = c(0,360)) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "cm"),
    panel.grid  = element_blank()
  )


###### Panel C Zeebie LM ########
panel_1C <- 
  all_lm_dir_width %>%
  filter(species == "zf") %>%
  #filter(investigator == "ag") %>%
  ggplot(aes(x = prefdir_sum)) +
  geom_hline(yintercept = c(0.10, 0.25), 
             colour = "grey90", size = 0.2) +
  geom_vline(xintercept = seq(0, 360-1, by = 90), 
             colour = "grey90", size = 0.2) +
  geom_histogram(aes(y = stat(width*density)), 
                 binwidth = 15, 
                 boundary = 0, 
                 fill = col_zb, 
                 color = "black",
                 size = 0.15) +
  coord_polar(direction = 1, start = pi/2) + 
  scale_y_continuous(trans = "sqrt", 
                     breaks = c(0.10, 0.25),
                     expand = c(0,0), limits = c(0, 0.42)) +
  scale_x_continuous(breaks = c(0, 90,180,270), expand = c(0,0), 
                     limits = c(0,360)) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "cm"),
    panel.grid  = element_blank()
  )

###### Panel D Pigeon LM ########
panel_1D <- 
  pg_LM_tuning_width_summary %>% 
  ggplot(aes(x = prefdir_sum)) +
  geom_hline(yintercept = c(0.10, 0.25), 
             colour = "grey90", size = 0.2) +
  geom_vline(xintercept = seq(0, 360-1, by = 90), 
             colour = "grey90", size = 0.2) +
  geom_histogram(aes(y = stat(width*density)), 
                 binwidth = 15, 
                 boundary = 0, 
                 fill = col_pg, 
                 color = "black",
                 size = 0.15) +
  coord_polar(direction = 1, start = pi/2) + 
  scale_y_continuous(trans = "sqrt", 
                     breaks = c(0.10, 0.25),
                     expand = c(0,0), limits = c(0, 0.42)) +
  scale_x_continuous(breaks = c(0, 90,180,270), expand = c(0,0), 
                     limits = c(0,360)) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "cm"),
    panel.grid  = element_blank()
  )

