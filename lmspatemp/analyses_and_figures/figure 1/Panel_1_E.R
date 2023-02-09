#### IMPORTANT NOTE #####
## All code in this script depends on objects that are created for various
## panels in Figure 1. For things to work here, the following line MUST
## be run, which will then create all the necessary objects in the environment.

#source("./analyses_and_figures/figure 1/Panel_1_BCD.R")

################ Panel 1E Polar unaligned tuning width curves ##################
## LM direction tuning width curves (spline fits)
## Firing rates are normalized so that min = 0 and max = 1
## Now in polar format
panel_1E <- 
  all_lm_dir_spline_fits_df %>%
  ggplot(
         aes(x = spline_x, 
             y = standardized_y, 
             color = species, 
             group = filename)) +
  geom_hline(yintercept = c(0, 1), 
             colour = "black", size = 0.2) +
  geom_vline(xintercept = seq(0, 360-1, by = 90), 
             colour = "grey90", size = 0.2) +
  geom_line(aes(alpha = species)) +
  scale_alpha_manual(values = c(0.09, 0.12, 0.09)) +
  scale_color_manual(values = c(col_hb, col_pg, col_zf)) +
  geom_quantile(quantiles = 0.5,
                aes(group = species),
                method = 'rqss',
                lambda = 50,
                size = 1.3, color = "black") +
  geom_quantile(quantiles = 0.5,
                aes(group = species),
                method = 'rqss',
                lambda = 50,
                size = 1) +
  #geom_smooth(aes(group = species)) +
  coord_polar(direction = 1, start = pi/2) + 
  scale_x_continuous(breaks = c(0, 90, 180, 270), 
                     expand = c(0, 0), 
                     #labels = c("F", "D", "B", "U"),
                     limits = c(0, 360)) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "pt"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid  = element_blank()
  )
