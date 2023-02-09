#### IMPORTANT NOTE #####
## All code in this script depends on objects that are created for various
## panels in Figure 1. For things to work here, create Fig 1 first.

#### __2O: LM ####
## LM direction tuning width curves
## "Peak direction" = max firing rate (Y-axis) of fitted spline
panel_S1A <- 
  ggplot(all_lm_dir_spline_fits_df,
         aes(x = trans_std_x, 
             y = trans_std_y, 
             color = species, 
             group = filename)) +
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
  scale_x_continuous(breaks = seq(-180, 180, by = 90),
                     limits = c(-180, 180),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25),
                     limits = c(0, 1.01),
                     expand = c(0, 0)) +
  geom_hline(aes(yintercept = -Inf)) +
  geom_vline(aes(xintercept = -Inf)) +
  labs(x = "° from peak direction", 
       y = "spikes/s (normalized)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.ticks = element_line(size = 0.3),
    axis.ticks.length = unit(0.1, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    axis.line = element_line(size = 0.3)
  )


#### __2P: LM SI boxplot ####
hb_zf_LM_dir_SI_vals <- 
  all_lm_tuning_width %>%
  select(species, si)

pg_LM_dir_SI_vals <- 
  pg_LM_tuning_width_summary %>%
  select(species, si)

LM_dir_SI_df <- 
  bind_rows(hb_zf_LM_dir_SI_vals, pg_LM_dir_SI_vals) %>%
  mutate(species = factor(species, levels = c("hb", "zf", "pg")))

panel_S1B <-   
  LM_dir_SI_df %>%
  ggplot(aes(x = species, 
             y = si, 
             fill = factor(species))) +  
  geom_boxplot(colour="black", size = 0.1, alpha = 0.3,
               outlier.size = 0, outlier.alpha = 0) +
  stat_boxplot(geom ='errorbar', size = 0.1, width = 0.5) + 
  geom_jitter(aes(fill = species, alpha = 0.3), 
              shape = 21, stroke = 0.1,
              width = 0.1, size = 1) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 1, by = 0.25), 
                     limits = c(0, 1)) +
  scale_color_manual(values = c(col_hb, col_zf, col_pg)) +
  scale_fill_manual(values = c(col_hb, col_zf, col_pg)) +
  xlab("species") +
  ylab("sensitivity index") +
  #coord_flip() +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.3)
  )

