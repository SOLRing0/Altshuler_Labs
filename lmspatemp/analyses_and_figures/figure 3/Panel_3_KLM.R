#### __3L counts panel ####
data_3K <- 
  read_csv("./data/spatemp/LM_ER_fullsweep_20bins/PZH - LM Combined Population spatemp Data - july 2021.csv") %>%
  mutate(log2SF = log2(SF),
         log2TF = log2(TF),
         q_cat = ifelse(zdiff > 1.65, "C", #Rind = C, #Rvel = A, Unk = B
                        ifelse(zdiff < -1.65, "A", "B")) 
  ) %>% 
  filter(Rind > -1) %>% 
  mutate(volume_n = 2*pi*sqrt(abs(Sfvar))*sqrt(abs(Tfvar))) %>%
  drop_na(zdiff) 

LM_q_cat_counts <-
  data_3K %>%
  group_by(Bird, q_cat) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = q_cat, values_from = n) %>%
  column_to_rownames(var = "Bird") %>%
  as.data.frame()

#chisq.test(LM_q_cat_counts)

p_3L_counts <-
  data_3K %>%
  group_by(Bird, q_cat) %>%
  summarize(n = n()) %>%
  mutate(Bird = as.factor(Bird))

p_3L_counts$percent <-
  round(100*p_3L_counts$n / c(61, 61, 61, 46, 46, 46, 61, 61, 61), 0) 

p_3L_counts$Bird <-
  factor(p_3L_counts$Bird, levels = c("HB", "ZF", "PG"))

Panel_3L <-
  ggplot(p_3L_counts,aes(y = n,
                         x = factor(Bird),
                         fill = factor(q_cat)))+
  geom_bar(position="fill", 
           stat = "identity")+
  geom_text(data = p_3L_counts,
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
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6), #element_blank(),
    axis.title = element_text(size = 7),
    axis.ticks = element_line(size = 0.3),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )

###### __3K - LM plot #####
p_5xE_ii <-
  ggplot(data_3K, aes(x = Rind, y = Rsp, group = Bird)) +
  geom_point(aes(fill = Bird), pch = 21) +
  scale_fill_manual(values = c(col_hb, col_pg, col_zf)) +
  #geom_density_2d(aes(color = Bird))+
  #scale_color_manual(values = c(col_hb, col_pg, col_zf)) +
  expand_limits(x = c(-1, 1),
                y = c(-1, 1)) +
  geom_hline(yintercept = 0.36) +
  geom_vline(xintercept = 0.36) +
  xlab("independent correlation (Rind)") +
  ylab("velocity correlation (Rvel)") +
  coord_fixed() +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6), #element_blank(),
    axis.title = element_text(size = 7),
    axis.ticks = element_line(size = 0.3),
    plot.margin = plot_margs,
    #panel.background = element_rect(fill = col_hb, colour = col_hb),
    legend.position="none"
  )

# Marginal densities along x axis
p_5xE_ii_xdens <- axis_canvas(p_5xE_ii, axis = "x") +
  geom_density(
    data = data_3K,
    aes(x = Rind, fill = Bird),
    alpha = 0.7,
    size = 0.2
  ) +
  scale_fill_manual(name = "Bird", values = c(col_hb, col_pg, col_zf))
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
p_5xE_ii_ydens <- axis_canvas(p_5xE_ii, axis = "y", coord_flip = TRUE) +
  geom_density(
    data = data_3K,
    aes(x = Rsp, fill = Bird),
    alpha = 0.7,
    size = 0.2
  ) +
  scale_fill_manual(name = "Bird", values = c(col_hb, col_pg, col_zf)) +
  coord_flip()
p_5xE_ii_p1 <-
  insert_xaxis_grob(p_5xE_ii, p_5xE_ii_xdens, 
                    grid::unit(.2, "null"), position = "top")
p_5xE_ii_p2 <-
  insert_yaxis_grob(p_5xE_ii_p1, p_5xE_ii_ydens, 
                    grid::unit(.2, "null"), position = "right")
Panel_3K <- ggdraw(p_5xE_ii_p2)

##### Hypothesis testing #####
## Step 1: reshape data
LM_spatempQ_reshaped <-
  reshape2::melt(
    data_3K,
    measure.vars = c("Rsp", "Rind")
  )

## Step 2: fit models
p3M_mod_sp <- 
  MCMCglmm::MCMCglmm(
    value ~ variable:Bird  - 1, #random = ~bird_id,
    data = LM_spatempQ_reshaped, scale = TRUE,
    nitt = 130000, thin = 100, burnin = 30000, 
    verbose = FALSE
  )

p3M_mod_sp_dat <- 
  p3M_mod_sp$Sol %>%
  as_tibble() %>%
  rename(`Rvel:hb`  = `variableRsp:BirdHB`, 
         `Rvel:pg`  = `variableRsp:BirdPG`, 
         `Rvel:zf`  = `variableRsp:BirdZF`,
         `Rind:hb` = `variableRind:BirdHB`, 
         `Rind:pg` = `variableRind:BirdPG`, 
         `Rind:zf` = `variableRind:BirdZF`) %>%
  gather()

#### __3M effects panel ####
Panel_3M <-
  ggplot(p3M_mod_sp_dat, aes(x = value, y = key, fill = key)) +
  stat_halfeye(.width = 0.95, slab_colour = "black", point_size = 0.15,
               interval_size = 0.5, slab_size = 0.25) +
  scale_fill_manual(values = c(col_hb, col_pg, col_zf, 
                               col_hb, col_pg, col_zf)) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_discrete(expand = c(0.02, 0))+
  xlab("correlation") +
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
