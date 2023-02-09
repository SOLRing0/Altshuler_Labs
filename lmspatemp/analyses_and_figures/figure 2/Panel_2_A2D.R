################################### Panel 2A ###################################
#geom_text_size = 1.5

## Import data
LM_gaussian_peak_data <- 
  read_csv(
    "./data/spatemp/LM_ER_fullsweep_20bins/PZH - LM Combined Population spatemp Data - july 2021.csv")
## Log transform
LM_gaussian_peak_data$logVel <- log2(LM_gaussian_peak_data$Velocity)
LM_gaussian_peak_data$logSF <- log2(LM_gaussian_peak_data$SF)
LM_gaussian_peak_data$logTF <- log2(LM_gaussian_peak_data$TF)

## Set up grouping
Panel_2A_data <-
  LM_gaussian_peak_data %>%
  mutate(edge_by_sp = paste(Edge_case, Bird, sep = "_"))

Panel_2A_color_map <- c(col_hb, col_pg, col_zb,
                        "black", "black", "black")
names(Panel_2A_color_map) <- c("Yes_HB", "Yes_PG", "Yes_ZF",
                               "No_HB", "No_PG", "No_ZF")

Panel_2A_fill_map <- c("white", "white", "white",
                       col_hb, col_pg, col_zb)
names(Panel_2A_fill_map) <- c("Yes_HB", "Yes_PG", "Yes_ZF",
                              "No_HB", "No_PG", "No_ZF")

Panel_2A_main <- 
  Panel_2A_data %>%
  ggplot(aes(x = logSF, y = logTF)) +
  geom_point(size = 1.65, aes(fill = factor(edge_by_sp), 
                              color = factor(edge_by_sp)), pch = 21) +
  scale_fill_manual(name = "edge_by_sp", values = Panel_2A_fill_map) +
  scale_color_manual(name = "edge_by_sp", values = Panel_2A_color_map) +
  expand_limits(y = c(log2(0.031), log2(16.5))) +
  scale_x_continuous(breaks = c(log2(0.0155), log2(0.062), log2(0.25),
                                log2(0.5) , log2(1)),
                     labels = c("0.0155", "0.062", "0.25", 
                                "0.5", "1"),
                     limits = c(log2(0.0150), log2(1.07))) + 
  scale_y_continuous(breaks = c(log2(0.031), log2(0.125), log2(0.5),
                                log2(2) , log2(8), log2(16)),
                     labels = c("0.031", "0.125", "0.5", 
                                "2", "8", "16"),
                     limits = c(log2(0.030), log2(16.5))) +
  coord_cartesian(xlim = c(log2(0.0155), log2(1.07)), 
                  ylim = c(log2(0.031), log2(16.5))) +
  xlab("spatial frequency (cpd)") +
  ylab("temporal frequency (Hz)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text=element_text(size=6),
        axis.title=element_text(size=7),
        axis.ticks = element_line(size = 0.3),
        axis.text.y = element_text(hjust = 0.5, angle = 90)
  )

## Add kernel densities along margins
Panel_5A_fill_map_edges <-
  c("#ED0080","#12bcff","#F48D00","#ED0080","#12bcff","#F48D00")
names(Panel_5A_fill_map_edges) <-
  c("Yes_HB", "Yes_PG", "Yes_ZF", "No_HB", "No_PG", "No_ZF")

Panel_5A_fill_map_edges_birds <-
  c("#ED0080","#12bcff","#F48D00")
names(Panel_5A_fill_map_edges_birds) <-
  c("HB", "PG", "ZF")


## Set up data sets
LM_spatemp_data <-
  Panel_2A_data %>%
  dplyr::select(Bird, Cell, logVel, logSF, logTF, 
                Edge_case, Rsp, Rind, zdiff, edge_by_sp) %>%
  mutate(nucleus = "LM") 

LM_spatemp_dat_no_edge <-
  Panel_2A_data %>%
  dplyr::select(Bird, Cell, logVel, logSF, logTF, 
                Edge_case, Rsp, Rind, zdiff, edge_by_sp) %>%
  mutate(nucleus = "LM") %>%
  filter(Edge_case == "No")

## Set up a relatively uninformative prior
biv_prior <- list(R = list(V = diag(2) * 0.8,
                           nu = 2),
                  G = list(G1 = list(V = diag(2) * 0.2,
                                     nu = 2))
)
univ_prior <- list(G = list(G1 = list(V = 1,
                                      nu = 0.2)),
                   R = list(V = 1, nu = 0.2))


# Marginal densities along x axis
xdens <- axis_canvas(Panel_2A_main, axis = "x") +
  geom_density(
    data = LM_spatemp_data,
    aes(x = logSF, fill = Bird),
    alpha = 0.7,
    size = 0.2
  ) +
  scale_fill_manual(name = "Bird", values = Panel_5A_fill_map_edges_birds)
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(Panel_2A_main, axis = "y", coord_flip = TRUE) +
  geom_density(
    data = LM_spatemp_data,
    aes(x = logTF, fill = Bird),
    alpha = 0.7,
    size = 0.2
  ) +
  scale_fill_manual(name = "Bird", values = Panel_5A_fill_map_edges_birds) +
  coord_flip()
p1 <-
  insert_xaxis_grob(Panel_2A_main, xdens, 
                    grid::unit(.2, "null"), position = "top")
p2 <-
  insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
Panel_2A_new_all <- ggdraw(p2)

#### __2A LM velocity density ####
Panel_2A_velocity <-
  ggplot(data = Panel_2A_data
         ) +
  geom_density(aes(x = logVel, fill = Bird),
               alpha = 0.7,
               size = 0.2) +
  scale_x_continuous(
    breaks = c(-2, 1, 4, 7),
    labels = c("0.25", "2", "16", "128"),
    expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + 
  xlab("velocity (°/s)") +
  ylab("density") +
  scale_fill_manual(name = "Bird", values = Panel_5A_fill_map_edges_birds) +
  theme_classic() +
  coord_flip() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 6),
    axis.title.x = element_text(size = 7),
    axis.ticks = element_line(size = 0.3),
    axis.title.y = element_text(size = 7, angle = 270),
    axis.text.y = element_text(
      size = 6,
      hjust = 0.5,
      angle = 270
    ),
    panel.background = element_rect(fill = "transparent"),
    # bg of the panel
    plot.background = element_rect(fill = "transparent",
                                   color = NA),
    # bg of the plot
    panel.grid.major = element_blank(),
    # get rid of major grid
    panel.grid.minor = element_blank(),
    # get rid of minor grid
    legend.background = element_rect(fill = "transparent"),
    # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") 
    # get rid of legend panel bg
  )

# ggsave("./analyses_and_figures/figure 2/Panel_2A_velocity.png", 
#        Panel_2A_velocity, width = 1, height = 1, units = "in",
#        dpi = "retina" #, bg = "transparent
# )
# 
# pdf(file = "./analyses_and_figures/figure 2/Panel_2A_velocity.pdf", 
#     width = 1, height = 1,
#     title = "Figure 5", paper = "letter", bg = "white",
#     pagecentre = TRUE, colormodel = "srgb")
# Panel_2A_velocity
# dev.off()

#### __Panel_2A_complete ####
p_2A_vel <-
  magick::image_read("./analyses_and_figures/figure 2/Panel_2A_velocity.png") %>%
  #magick::image_crop("310x297") %>%
  magick::image_rotate(-45) %>%
  image_trim()
Panel_2A_complete <- 
  plot_grid(Panel_2A_new_all, NULL,
            nrow = 1, rel_widths = c(1, 0.2)) + 
  draw_image(p_2A_vel,
             x = 0.37,
             y = 0.42,
             scale = 0.3)

#### Hypothesis testing for spatemp ####
#### __model velocity LM ####
p2_vel_mod_lm <-
  MCMCglmm::MCMCglmm(
    logVel ~ Bird - 1,
    #random =  ~ bird_id,
    data = LM_spatemp_data,
    nitt = 13000, thin = 10, burnin = 3000, 
    verbose = FALSE
  )
#summary(p2_vel_mod_lm)

p2_vel_mod_lm_dat <- 
  p2_vel_mod_lm$Sol %>%
  as_tibble() %>%
  rename(hb = BirdHB, pg = BirdPG, zf = BirdZF) %>%
  gather() %>%
  mutate(key = as_factor(key),
         key = fct_relevel(key, "hb", "zf", "pg"))


##### __model LM sf-tf ####
## Step 1: reshape data
LM_spatemp_reshaped <-
  reshape2::melt(
    LM_spatemp_data,
    measure.vars = c("logSF", "logTF")
  )

## Step 2: fit models
p2A_mod_sp <- 
  MCMCglmm::MCMCglmm(
    value ~ variable:Bird  - 1, #random = ~bird_id,
    data = LM_spatemp_reshaped, scale = TRUE,
    nitt = 130000, thin = 100, burnin = 30000, 
    verbose = FALSE
  )
#summary(p2A_mod_sp)

## Step 3: plot effects from best model
p2A_model_dat <- 
  p2A_mod_sp$Sol %>%
  as_tibble() %>%
  rename(`SF:HB` = `variablelogSF:BirdHB`, 
         `SF:PG` = `variablelogSF:BirdPG`, 
         `SF:ZF` = `variablelogSF:BirdZF`,
         `TF:HB` = `variablelogTF:BirdHB`, 
         `TF:PG` = `variablelogTF:BirdPG`, 
         `TF:ZF` = `variablelogTF:BirdZF`) %>%
  gather() 

## Now split that tibble based on SF vs. TF effects
p2A_model_dat_SF <-
  p2A_model_dat %>%
  filter(!grepl("TF", key)) %>%
  mutate(across("key", str_replace, "SF:HB", "hb")) %>%
  mutate(across("key", str_replace, "SF:PG", "pg")) %>%
  mutate(across("key", str_replace, "SF:ZF", "zf")) %>%
  mutate(key = as_factor(key),
         key = fct_relevel(key, "hb", "zf", "pg"))

p2A_model_dat_TF <-
  p2A_model_dat %>%
  filter(!grepl("SF", key)) %>%
  mutate(across("key", str_replace, "TF:HB", "hb")) %>%
  mutate(across("key", str_replace, "TF:PG", "pg")) %>%
  mutate(across("key", str_replace, "TF:ZF", "zf")) %>%
  mutate(key = as_factor(key),
         key = fct_relevel(key, "hb", "zf", "pg"))

#### set panel ranges for velocity ####
## min vel based on nbor slow
models_vel_min <- min(p2_vel_mod_lm_dat$value) - 0.01
## max vel based on nbor fast
models_vel_max <- max(p2_vel_mod_lm_dat$value) + 0.01

#### set panel ranges for SF TF ####
## min sf based on nbor fast
models_sf_min <- min(p2A_model_dat_SF$value) - 0.01
## max sf based on nbor slow
models_sf_max <- max(p2A_model_dat_SF$value) + 0.01
## min tf based on nbor slow
models_tf_min <- min(p2A_model_dat_SF$value) - 0.01
## max tf based on nbor fast
models_tf_max <- max(p2A_model_dat_SF$value) + 0.01


#### __Panel_2B_LM_vel ####
Panel_2B_LM_vel <-
  ggplot(p2_vel_mod_lm_dat, aes(x = value, y = key, fill = key)) +
  stat_halfeye(.width = 0.95, slab_colour = "black", point_size = 0.15,
               interval_size = 0.5, slab_size = 0.25) +
  scale_fill_manual(values = c(col_hb, col_zf, col_pg)) +
  scale_x_continuous(breaks = c(1, 4, 7),
                     labels = c("2", "16", "128"),
                     limits = c(1, 7)) +
  scale_y_discrete(expand = c(0.02, 0))+
  xlab("velocity (°/s)") +
  ylab("species effect") +
  theme_minimal() +
  theme(
    legend.position = 'none',
    axis.ticks.x = element_line(size = 0.3),
    axis.ticks.y = element_blank(),
    axis.ticks.length = unit(0.1, "cm"),
    axis.text = element_text(size = 6),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(size = 7)
  )


#### __2C LM SF ####
Panel_2C_LM_SF <-
  ggplot(p2A_model_dat_SF, aes(x = value, y = key, fill = key)) +
  stat_halfeye(.width = 0.95, slab_colour = "black", point_size = 0.15,
               interval_size = 0.5, slab_size = 0.25) +
  scale_fill_manual(values = c(col_hb, col_zf, col_pg,
                               col_hb, col_zf, col_pg)) +
  scale_x_continuous(breaks = c(-6, -4, -2),
                     labels = c("0.016", "0.062", "0.25"),
                     limits = c(-6.1, -1.4)) +
  scale_y_discrete(expand = c(0.02, 0))+
  xlab("spatial frequency (cpd)") +
  ylab(" ") +
  theme_minimal() +
  theme(
    legend.position = 'none',
    axis.ticks.x = element_line(size = 0.3),
    axis.ticks.y = element_blank(),
    axis.ticks.length = unit(0.1, "cm"),
    axis.text = element_text(size = 6),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(size = 7)
  )

#### __2D LM TF ####
Panel_2D_LM_TF <-
  ggplot(p2A_model_dat_TF, aes(x = value, y = key, fill = key)) +
  stat_halfeye(.width = 0.95, slab_colour = "black", point_size = 0.15,
               interval_size = 0.5, slab_size = 0.25) +
  scale_fill_manual(values = c(col_hb, col_zf, col_pg,
                               col_hb, col_zf, col_pg)) +
  scale_x_continuous(breaks = c(-1, 1, 3),
                     labels = c("0.5", "2", "8"),
                     limits = c(-1.1, 3.1)) +
  scale_y_discrete(expand = c(0.02, 0))+
  xlab("temporal frequency (Hz)") +
  ylab(" ") +
  theme_minimal() +
  theme(
    legend.position = 'none',
    axis.ticks.x = element_line(size = 0.3),
    axis.ticks.y = element_blank(),
    axis.ticks.length = unit(0.1, "cm"),
    axis.text = element_text(size = 6),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(size = 7)
  )

