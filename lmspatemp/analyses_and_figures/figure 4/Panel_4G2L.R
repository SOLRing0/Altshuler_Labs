## geom_point sizing
p4g2l_dotsize <- 0.8
p4g2l_stroke <- 0.5

##### Import data ####
LM_initial_trans <-
  read_xlsx("./data/spatemp/LM_ER_phasic_summaries/ER_initial_transients_summary.xlsx") %>%
  mutate(edge_by_sp = paste(Edge_case, Bird, sep = "_")) %>%
  filter(R2_un > 0.5) %>% # keep those that have good gaussian fits
  #filter(Edge_case == "No") %>%
  drop_na(Q_un) %>% # just in case
  add_column(source = "IT")
LM_steady_state <-
  read_xlsx("./data/spatemp/LM_ER_phasic_summaries/ER_steady_state_summary.xlsx") %>%
  mutate(edge_by_sp = paste(Edge_case, Bird, sep = "_")) %>%
  filter(R2_un > 0.5) %>% # keep those that have good gaussian fits
  #filter(Edge_case == "No") %>%
  drop_na(Q_un) %>% # just in case
  add_column(source = "SS")

LM_phasic_data_wide <-
  left_join(LM_initial_trans,
            LM_steady_state,
            by = c("name", "Bird"),
            suffix = c(".IT", ".SS")) %>%
  ## describe full condition
  mutate(phase_edges = 
           paste(Edge_case.IT, source.IT, 
                 Edge_case.SS, source.SS, 
                 sep = "_"))
LM_phase_edges <- 
  LM_phasic_data_wide %>%
  select(name, Bird, phase_edges)
## No_IT_No_SS = Full Gaussian for both transient and steady state
## Yes_IT_No_SS = Edge case for transient only
## No_IT_Yes_SS = Edge case for steady only
## Yes_IT_Yes_SS = Edge case for both phases

## Which cells are edge cases in either IT or SS?
# summary(as.factor(LM_phase_edges$phase_edges))

## which cells appear in both?
ssit_intersect <- intersect(LM_initial_trans$name, LM_steady_state$name)
LM_phasic_data_tall <-
  bind_rows(LM_initial_trans[LM_initial_trans$name %in% ssit_intersect,], 
            LM_steady_state[LM_steady_state$name %in% ssit_intersect,]) %>%
  left_join(LM_phase_edges, by = c("name", "Bird"))

## sort by phasic characteristic
## this will make plotting better
LM_phasic_data_tall <-
  arrange(LM_phasic_data_tall, desc(phase_edges))
LM_phasic_data_wide <-
  arrange(LM_phasic_data_wide, desc(phase_edges))


##### Sort out colors #####
## Black lines and dot encircling indicate full gaussians for both IT and SS
## Grey fill = IT; species-color fill = SS
  
p_4G_color_map <- c("black", ## no edge cases
                    alpha("grey80", 0), 
                    alpha("grey80", 0), 
                    alpha("grey80", 0)
                    )
names(p_4G_color_map) <- c("No_IT_No_SS", ## no edge cases
                           "Yes_IT_No_SS", 
                           "No_IT_Yes_SS", 
                           "Yes_IT_Yes_SS")

p_4G_hb_fill_map <- c("grey80", col_hb)
names(p_4G_hb_fill_map) <- c("IT", "SS")
p_4G_zf_fill_map <- c("grey80", col_zf)
names(p_4G_zf_fill_map) <- c("IT", "SS")

##### Panel 4G paired peaks plots ######
###### __p_4G_hb #####
p_4G_hb <-
  LM_phasic_data_tall %>%
  filter(Bird == "HB") %>%
  ggplot(aes(x = SF_peak_un, y = TF_peak_un, group = name)) +
  geom_line(aes(linetype = phase_edges, color = phase_edges), 
            alpha = 1) +
  geom_point(aes(fill = source, color = phase_edges), 
             #size = p4g2l_dotsize, 
             pch = 21) +
  scale_linetype_manual(values = c("solid", "dotted", "dotted")) +
  scale_fill_manual(name = "Bird", values = p_4G_hb_fill_map) +
  scale_color_manual(name = "source", values = p_4G_color_map) +
  xlab("spatial frequency (cpd)") +
  ylab(expression(paste("temporal frequency (Hz)"))) +
  scale_y_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     labels = c(two_neg_five, two_neg_three, two_neg_one, 
                                two_one, two_three, two_four),
                     limits = c(-5, 4)) +
  scale_x_continuous(breaks = seq(-6, 0, by = 1), 
                     labels = c(two_neg_six, two_neg_five, two_neg_four,
                                two_neg_three,
                                two_neg_two, two_neg_one, two_zero),
                     limits = c(-6, 0)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    legend.position = "none"
  )

# ## how many No_IT_No_SS cells total here (black lines)?
# LM_phasic_data_tall %>%
#   filter(Bird == "HB") %>%
#   select(phase_edges) %>%
#   as_vector() %>%
#   as.factor() %>%
#   summary() / 2


###### __p_4G_zf #####
p_4G_zf <-
  LM_phasic_data_tall %>%
  filter(Bird == "ZF") %>%
  ggplot(aes(x = SF_peak_un, y = TF_peak_un, group = name)) +
  geom_line(aes(linetype = phase_edges, color = phase_edges), 
            alpha = 1) +
  geom_point(aes(fill = source, 
                 color = phase_edges), 
             #size = p4g2l_dotsize, 
             pch = 21) +
  scale_linetype_manual(values = c("solid", "dotted", "dotted", "dotted")) +
  scale_fill_manual(name = "Bird", values = p_4G_zf_fill_map) +
  scale_color_manual(name = "source", values = p_4G_color_map) +
  xlab(" ") +
  ylab(expression(paste(" "))) +
  # expand_limits(x = c(-6, 0),
  #               y = c(-5, 4)) +
  scale_y_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     # labels = c(two_neg_five, two_neg_three, two_neg_one, 
                     #            two_one, two_three, two_four),
                     labels = c(rep("", 6)),
                     limits = c(-5, 4)) +
  scale_x_continuous(breaks = seq(-6, 0, by = 1), 
                     labels = c(two_neg_six, two_neg_five, two_neg_four,
                                two_neg_three,
                                two_neg_two, two_neg_one, two_zero),
                     limits = c(-6, 0)) +
  # coord_cartesian(xlim = c(-6, 0), 
  #                 ylim = c(-5, 4),
  #                 expand = FALSE) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 6),
    axis.title = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    legend.position="none"
  )

# ## how many No_IT_No_SS cells total here (black lines)?
# LM_phasic_data_tall %>%
#   filter(Bird == "ZF") %>%
#   select(phase_edges) %>%
#   as_vector() %>%
#   as.factor() %>%
#   summary() / 2

##### Paired peak hypothesis testing ######
phase_vel_test <- LM_phasic_data_tall %>%
  group_by(Bird) %>%
  pairwise_t_test(
    LnVel_un ~ source,
    paired = TRUE,
    p.adjust.method = "bonferroni",
    detailed = TRUE
  )
phase_vel_test

phase_SF_test <- LM_phasic_data_tall %>%
  group_by(Bird) %>%
  pairwise_t_test(
    SF_peak_un ~ source,
    paired = TRUE,
    p.adjust.method = "bonferroni",
    detailed = TRUE
  )
phase_SF_test

phase_TF_test <- LM_phasic_data_tall %>%
  group_by(Bird) %>%
  pairwise_t_test(
    TF_peak_un ~ source,
    paired = TRUE,
    p.adjust.method = "bonferroni",
    detailed = TRUE
  )
phase_TF_test

phase_peaks_pvals <-
  c(phase_vel_test$p, phase_SF_test$p, phase_TF_test$p)
format(phase_peaks_pvals, scientific = FALSE)

######## Panel 4H connected dot plots ########
## get y-axis minima and maxima, then add small adjustment to expand bounds
xpans <- 0.1
p_4H_vel_min <- min(LM_phasic_data_tall$LnVel_un) - xpans
p_4H_vel_max <- max(LM_phasic_data_tall$LnVel_un) + xpans
p_4H_SF_min <- min(LM_phasic_data_tall$SF_peak_un) - xpans
p_4H_SF_max <- max(LM_phasic_data_tall$SF_peak_un) + xpans
p_4H_TF_min <- min(LM_phasic_data_tall$TF_peak_un) - xpans
p_4H_TF_max <- max(LM_phasic_data_tall$TF_peak_un) + xpans

#### __p_4H_vel_hb ####
p_4H_vel_hb <- 
  LM_phasic_data_tall %>%
  filter(Bird == "HB") %>%
  ggplot(aes(x = source, y = LnVel_un)) +
  geom_boxplot(width = 0.75) +
  geom_line(aes(linetype = phase_edges, 
                color = phase_edges,
                group = name), 
            alpha = 1, size = 0.1) +
  geom_point(aes(fill = source,
                 color = phase_edges), 
             ##size = p4g2l_dotsize,  
             pch = 21) +
  scale_linetype_manual(values = c("solid", "dotted", "dotted", "dotted")) +
  scale_fill_manual(name = "Bird", values = p_4G_hb_fill_map) +
  scale_color_manual(name = "source", values = p_4G_color_map) +
  expand_limits(y = c(p_4H_vel_min, p_4H_vel_max)) +
  scale_y_continuous(breaks = c(0, 3, 6), 
                     labels = c(two_zero, two_three, two_six)) +
  ylab("velocity (°/s)") +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 6), #element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    #panel.background = element_rect(fill = col_hb, colour = col_hb),
    legend.position="none"
  )
#### __p_4H_vel_zf ####
p_4H_vel_zf <- 
  LM_phasic_data_tall %>%
  filter(Bird == "ZF") %>%
  ggplot(aes(x = source, y = LnVel_un)) +
  geom_boxplot(width = 0.75) +
  geom_line(aes(linetype = phase_edges, 
                color = phase_edges,
                group = name), 
            alpha = 1, size = 0.1) +
  geom_point(aes(fill = source,
                 color = phase_edges), 
             #size = p4g2l_dotsize, 
             pch = 21) +
  scale_linetype_manual(values = c("solid", "dotted", "dotted", "dotted")) +
  scale_fill_manual(name = "Bird", values = p_4G_zf_fill_map) +
  scale_color_manual(name = "source", values = p_4G_color_map) +
  scale_y_continuous(breaks = c(0, 3, 6), 
                     labels = c(two_zero, two_three, two_six)) +
  ylab(" ") +
  expand_limits(y = c(p_4H_vel_min, p_4H_vel_max)) +   
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 6),
    axis.title = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    #panel.background = element_rect(fill = col_hb, colour = col_hb),
    legend.position="none"
  )

#### __p_4H_SF_hb ####
p_4H_SF_hb <- 
  LM_phasic_data_tall %>%
  filter(Bird == "HB") %>%
  ggplot(aes(x = source, y = SF_peak_un)) +
  geom_boxplot(width = 0.75) +
  geom_line(aes(linetype = phase_edges, 
                color = phase_edges,
                group = name), 
            alpha = 1, size = 0.1) +
  geom_point(aes(fill = source, 
                 color = phase_edges), 
             #size = p4g2l_dotsize, 
             pch = 21) +
  scale_linetype_manual(values = c("solid", "dotted", "dotted", "dotted")) +
  scale_fill_manual(name = "Bird", values = p_4G_hb_fill_map) +
  scale_color_manual(name = "source", values = p_4G_color_map) +
  scale_y_continuous(breaks = seq(-6, 0, by = 1), 
                     labels = c(two_neg_six, two_neg_five, two_neg_four,
                                two_neg_three,
                                two_neg_two, two_neg_one, two_zero),
                     limits = c(-6, 0)) +
  expand_limits(y = c(p_4H_SF_min, p_4H_SF_max)) +   
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 6), #element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    #panel.background = element_rect(fill = col_hb, colour = col_hb),
    legend.position="none"
  )
#### __p_4H_SF_zf ####
p_4H_SF_zf <- 
  LM_phasic_data_tall %>%
  filter(Bird == "ZF") %>%
  ggplot(aes(x = source, y = SF_peak_un)) +
  geom_boxplot(width = 0.75) +
  geom_line(aes(linetype = phase_edges, 
                color = phase_edges,
                group = name), 
            alpha = 1, size = 0.1) +
  geom_point(aes(fill = source, 
                 color = phase_edges), 
             #size = p4g2l_dotsize, 
             pch = 21) +
  scale_linetype_manual(values = c("solid", "dotted", "dotted", "dotted")) +
  scale_fill_manual(name = "Bird", values = p_4G_zf_fill_map) +
  scale_color_manual(name = "source", values = p_4G_color_map) +
  scale_y_continuous(breaks = seq(-6, 0, by = 1), 
                     labels = c(two_neg_six, two_neg_five, two_neg_four,
                                two_neg_three,
                                two_neg_two, two_neg_one, two_zero),
                     limits = c(-6, 0)) +
  expand_limits(y = c(p_4H_SF_min, p_4H_SF_max)) +   
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 6), #element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    #panel.background = element_rect(fill = col_hb, colour = col_hb),
    legend.position="none"
  )

#### __p_4H_TF_hb ####
p_4H_TF_hb <- 
  LM_phasic_data_tall %>%
  filter(Bird == "HB") %>%
  ggplot(aes(x = source, y = TF_peak_un)) +
  geom_boxplot(width = 0.75) +
  geom_line(aes(linetype = phase_edges, 
                color = phase_edges,
                group = name), 
            alpha = 1, size = 0.1) +
  geom_point(aes(fill = source, 
                 color = phase_edges), 
             #size = p4g2l_dotsize, 
             pch = 21) +
  scale_linetype_manual(values = c("solid", "dotted", "dotted", "dotted")) +
  scale_fill_manual(name = "Bird", values = p_4G_hb_fill_map) +
  scale_color_manual(name = "source", values = p_4G_color_map) +
  scale_y_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     # labels = c(two_neg_five, two_neg_three, two_neg_one, 
                     #            two_one, two_three, two_four),
                     labels = c(rep("", 6)),
                     limits = c(-5, 4)) +
  expand_limits(y = c(p_4H_TF_min, p_4H_TF_max)) +  
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 6), #element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    #panel.background = element_rect(fill = col_hb, colour = col_hb),
    legend.position="none"
  )
#### __p_4H_TF_zf ####
p_4H_TF_zf <- 
  LM_phasic_data_tall %>%
  filter(Bird == "ZF") %>%
  ggplot(aes(x = source, y = TF_peak_un)) +
  geom_boxplot(width = 0.75) +
  geom_line(aes(linetype = phase_edges, 
                color = phase_edges,
                group = name), 
            alpha = 1, size = 0.1) +
  geom_point(aes(fill = source, 
                 color = phase_edges), 
             #size = p4g2l_dotsize, 
             pch = 21) +
  scale_linetype_manual(values = c("solid", "dotted", "dotted", "dotted")) +
  scale_fill_manual(name = "Bird", values = p_4G_zf_fill_map) +
  scale_color_manual(name = "source", values = p_4G_color_map) +
  scale_y_continuous(breaks = c(-5, -3, -1,  1,  3,  4),
                     # labels = c(two_neg_five, two_neg_three, two_neg_one, 
                     #            two_one, two_three, two_four),
                     labels = c(rep("", 6)),
                     limits = c(-5, 4)) +
  expand_limits(y = c(p_4H_TF_min, p_4H_TF_max)) +  
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 6), #element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    #panel.background = element_rect(fill = col_hb, colour = col_hb),
    legend.position="none"
  )

#### all plots, quick view ###
# plot_grid(p_4H_vel_hb,
#           p_4H_SF_hb,
#           p_4H_TF_hb,
#           p_4H_vel_zf,
#           p_4H_SF_zf,
#           p_4H_TF_zf)

##### Orientation data ######
LM_phasic_q_data <-
  LM_phasic_data_tall %>%
  filter(phase_edges == "No_IT_No_SS")
  
##### Panel 4I paired orientation plots ######
###### __p_4I_hb #####
p_4I_hb <-
  LM_phasic_q_data %>%
  filter(Bird == "HB") %>%
  drop_na(Z_diff) %>%
  ggplot(aes(x = R_indp, y = R_diag, group = source)) +
  geom_line(aes(group = name, alpha = 0.1)) +
  geom_point(aes(fill = source,
                 color = source), 
             #size = p4g2l_dotsize, 
             pch = 21) +
  scale_linetype_manual(values = c("solid")) +
  scale_fill_manual(name = "Bird", values = p_4G_hb_fill_map) +
  scale_color_manual(name = "source", values = c("black", "black")) +
  expand_limits(x = c(-1, 1),
                y = c(-1, 1)) +
  xlab("independent correlation (Rind)") +
  ylab("velocity correlation (Rvel)") +
  coord_fixed() +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    legend.position = "none"
  )

###### __p_4I_zf #####
p_4I_zf <-
  LM_phasic_q_data %>%
  filter(Bird == "ZF") %>%
  drop_na(Z_diff) %>%
  ggplot(aes(x = R_indp, y = R_diag, group = source)) +
  geom_line(aes(group = name, alpha = 0.1)) +
  geom_point(aes(fill = source,
                 color = source), 
             #size = p4g2l_dotsize, 
             pch = 21) +
  scale_linetype_manual(values = c("solid")) +
  scale_fill_manual(name = "Bird", values = p_4G_zf_fill_map) +
  scale_color_manual(name = "source", values = c("black", "black")) +
  expand_limits(x = c(-1, 1),
                y = c(-1, 1)) +
  xlab(" ") +
  ylab(" ") +
  coord_fixed() +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    legend.position = "none"
  )

##### Paired orientation hypothesis testing ######
phase_Rind_test <- 
  LM_phasic_q_data %>%
  group_by(Bird) %>%
  pairwise_t_test(
    R_indp ~ source,
    paired = TRUE,
    p.adjust.method = "bonferroni",
    detailed = TRUE
  )
phase_Rind_test

phase_Rvel_test <- 
  LM_phasic_q_data %>%
  group_by(Bird) %>%
  pairwise_t_test(
    R_diag ~ source,
    paired = TRUE,
    p.adjust.method = "bonferroni",
    detailed = TRUE
  )
phase_Rvel_test

phase_Qs_pvals <-
  c(phase_Rind_test$p, phase_Rvel_test$p)
format(phase_Qs_pvals, scientific = FALSE)

  
######## Panel 4J connected dot plots ########
xpans <- 0.1
p_4J_Rvel_min <- min(LM_phasic_q_data$R_diag) - xpans
p_4J_Rind_min <- min(LM_phasic_q_data$R_indp) - xpans
p_4J_Yval_min <- min(p_4J_Rvel_min, p_4J_Rind_min)
p_4J_Qs_max <- 1

#### __p_4J_Rvel_hb ####
p_4J_Rvel_hb <- 
  LM_phasic_q_data %>%
  filter(Bird == "HB") %>%
  ggplot(aes(x = source, y = R_diag)) +
  geom_boxplot(width = 0.75) +
  geom_line(aes(linetype = phase_edges, 
                color = phase_edges,
                group = name), 
            alpha = 1, size = 0.1) +
  geom_point(aes(fill = source, 
                 color = phase_edges), 
             #size = p4g2l_dotsize, 
             pch = 21) +
  scale_linetype_manual(values = c("solid", "dotted", "dotted", "dotted")) +
  scale_fill_manual(name = "Bird", values = p_4G_hb_fill_map) +
  scale_color_manual(name = "source", values = p_4G_color_map) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  expand_limits(y = c(p_4J_Yval_min, p_4J_Qs_max)) +   
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_line(size = 0.3)
  )
#### __p_4J_Rvel_zf ####
p_4J_Rvel_zf <- 
  LM_phasic_q_data %>%
  filter(Bird == "ZF") %>%
  ggplot(aes(x = source, y = R_diag)) +
  geom_boxplot(width = 0.75) +
  geom_line(aes(linetype = phase_edges, 
                color = phase_edges,
                group = name), 
            alpha = 1, size = 0.1) +
  geom_point(aes(fill = source, 
                 color = phase_edges), 
             #size = p4g2l_dotsize, 
             pch = 21) +
  scale_linetype_manual(values = c("solid", "dotted", "dotted", "dotted")) +
  scale_fill_manual(name = "Bird", values = p_4G_zf_fill_map) +
  scale_color_manual(name = "source", values = p_4G_color_map) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  expand_limits(y = c(p_4J_Yval_min, p_4J_Qs_max)) +   
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_line(size = 0.3)
  )
#### __p_4J_Rind_hb ####
p_4J_Rind_hb <- 
  LM_phasic_q_data %>%
  filter(Bird == "HB") %>%
  ggplot(aes(x = source, y = R_indp)) +
  geom_boxplot(width = 0.75) +
  geom_line(aes(linetype = phase_edges, 
                color = phase_edges,
                group = name), 
            alpha = 1, size = 0.1) +
  geom_point(aes(fill = source,
                 color = phase_edges),
             #size = p4g2l_dotsize, 
             pch = 21) +
  scale_linetype_manual(values = c("solid", "dotted", "dotted", "dotted")) +
  scale_fill_manual(name = "Bird", values = p_4G_hb_fill_map) +
  scale_color_manual(name = "source", values = p_4G_color_map) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  expand_limits(y = c(p_4J_Yval_min, p_4J_Qs_max)) +   
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_line(size = 0.3)
  )
#### __p_4J_Rind_zf ####
p_4J_Rind_zf <- 
  LM_phasic_q_data %>%
  filter(Bird == "ZF") %>%
  ggplot(aes(x = source, y = R_indp)) +
  geom_boxplot(width = 0.75) +
  geom_line(aes(linetype = phase_edges, 
                color = phase_edges,
                group = name), 
            alpha = 1, size = 0.1) +
  geom_point(aes(fill = source,
                 color = phase_edges),
             #size = p4g2l_dotsize, 
             pch = 21) +
  scale_linetype_manual(values = c("solid", "dotted", "dotted", "dotted")) +
  scale_fill_manual(name = "Bird", values = p_4G_zf_fill_map) +
  scale_color_manual(name = "source", values = p_4G_color_map) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  expand_limits(y = c(p_4J_Yval_min, p_4J_Qs_max)) +   
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_line(size = 0.3)
  )

##### Category count munging ######
IT_counts <-
  LM_phasic_q_data %>%
  filter(source == "IT") %>%
  mutate(q_cat = ifelse(Z_diff > 1.65, "C", #Rind = C, #Rvel = A, Unk = B
                        ifelse(Z_diff < -1.65, "A", "B"))) %>%
  filter(R_indp > -1) %>%
  mutate(volume_n = 2 * pi * sqrt(abs(SF_sig_un)) * sqrt(abs(TF_sig_un))) %>%
  drop_na(Z_diff) %>%
  group_by(Bird, q_cat) %>%
  summarize(n = n())
SS_counts <-
  LM_phasic_q_data %>%
  filter(source == "SS") %>%
  mutate(q_cat = ifelse(Z_diff > 1.65, "C", #Rind = C, #Rvel = A, Unk = B
                        ifelse(Z_diff < -1.65, "A", "B"))) %>%
  filter(R_indp > -1) %>%
  mutate(volume_n = 2 * pi * sqrt(abs(SF_sig_un)) * sqrt(abs(TF_sig_un))) %>%
  drop_na(Z_diff) %>%
  group_by(Bird, q_cat) %>%
  summarize(n = n())

IT_counts_w <-
  IT_counts %>%
  pivot_wider(names_from = q_cat, values_from = n) %>%
  column_to_rownames(var = "Bird") %>%
  as.data.frame()
SS_counts_w <-
  SS_counts %>%
  pivot_wider(names_from = q_cat, values_from = n) %>%
  column_to_rownames(var = "Bird") %>%
  as.data.frame()

phasic_counts <-
  bind_rows(IT_counts_w, SS_counts_w)
chisq_test(phasic_counts)
chisq.test(IT_counts_w)
chisq.test(SS_counts_w)

IT_counts$percent <-
  round(100*IT_counts$n / c(22, 22, 22, 45, 45, 45), 0) 
SS_counts$percent <-
  round(100*SS_counts$n / c(22, 22, 22, 45, 45, 45), 0) 

IT_counts$Bird <-
  factor(IT_counts$Bird, levels = c("HB", "ZF"))
SS_counts$Bird <-
  factor(SS_counts$Bird, levels = c("HB", "ZF"))

#### Panel 4K ####
Panel_4K <-
  ggplot(IT_counts,aes(y = n,
                         x = factor(Bird),
                         fill = factor(q_cat)))+
  geom_bar(position="fill", 
           stat = "identity")+
  geom_text(data = IT_counts,
            aes(y = n, label = percent),
            position = position_fill(vjust = 0.5),
            size = geom_text_size) +
  scale_y_continuous(labels = scales::percent,
                     breaks = c(0, 0.5, 1),
                     expand = c(0, 0)) +
  xlab('') + ylab('percent of cells') +
  scale_fill_manual(values = 
                      #CBB8D7 = purple
                      #E6E6E6 = grey
                      #B8D6BE = green
                      c("#B8D6BE", "#E6E6E6", "#CBB8D7")) +
  #coord_flip() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 6),
    axis.title = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )

#### Panel 4L ####
Panel_4L <-
  ggplot(SS_counts,aes(y = n,
                       x = factor(Bird),
                       fill = factor(q_cat)))+
  geom_bar(position="fill", 
           stat = "identity")+
  geom_text(data = SS_counts,
            aes(y = n, label = percent),
            position = position_fill(vjust = 0.5),
            size = geom_text_size) +
  scale_y_continuous(labels = scales::percent,
                     breaks = c(0, 0.5, 1),
                     expand = c(0, 0)) +
  xlab('') + ylab('percent of cells') +
  scale_fill_manual(values = 
                      #CBB8D7 = purple
                      #E6E6E6 = grey
                      #B8D6BE = green
                      c("#B8D6BE", "#E6E6E6", "#CBB8D7")) +
  #coord_flip() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 6),
    axis.title = element_blank(),
    axis.ticks = element_line(size = 0.3),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )
