cols <- c("#ED0080" = "hb", "#F48D00" = "zb", "#12bcff" = "pg")

###### Panel 1M #######
panel_1M_LM_data <- 
  read_csv("./data/spatemp/LM_ER_fullsweep_20bins/PZH - LM Combined Population spatemp Data - july 2021.csv")

panel_1M_LM_filtered <-
  panel_1M_LM_data %>%
  filter(Edge_case == "No") %>% ## Remove all cells at the boundaries of SF, TF
  #filter(`number of peaks` == 1) %>% ## Only primary peaks allowed
  filter(Fit_procedure == "unconstrained") %>% ## Remove any constrained fits
  mutate(volume_n = 2*pi*sqrt(abs(Sfvar))*sqrt(abs(Tfvar)))

panel_1M <-
  panel_1M_LM_filtered %>%
  ggplot(aes(x = Bird, y = volume_n,
             fill = factor(Bird),
             group = factor(Bird))) + 
  geom_boxplot(colour="black", size = 0.1, alpha = 0.3,
               outlier.size = 0, outlier.alpha = 0) +
  stat_boxplot(geom ='errorbar', size = 0.1, width = 0.5) + 
  # geom_jitter(aes(fill = Bird, alpha = 0.3), 
  #             shape = 21, stroke = 0.1,
  #             width = 0.1, size = 1) +
  geom_dotplot(
    binaxis = 'y',
    stackdir = 'center',
    method = 'histodot',
    dotsize = 0.5,
    stroke = 0
  ) +
  expand_limits(y = c(0, 7.5))  +
  scale_y_continuous(limits = c(0, 40),
    breaks = c(0, 20, 40)) +
  scale_fill_manual(values = c(col_hb,col_pg,col_zb))+
  scale_color_manual(values = c(col_hb,col_pg,col_zb))+
  scale_x_discrete(limits=c("HB", "ZF", "PG")) +
  xlab("species") +
  ylab("volume under gaussian fit") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 6),
        axis.title.x = element_text(vjust = -1),
        axis.ticks = element_line(size = 0.3),
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.3)
  )

## Species-specific data sets
HB_LM_gaussian_summaries <- 
  panel_1M_LM_data %>%
  filter(Edge_case == "No") %>% ## Remove all cells at the boundaries of SF, TF
  #filter(`number of peaks` == 1) %>% ## Only primary peaks allowed
  filter(Fit_procedure == "unconstrained") %>% ## Remove any constrained fits
  filter(Bird == "HB")
#length(unique(HB_LM_gaussian_summaries$Ind)) ## number of individuals
#nrow(HB_LM_gaussian_summaries) ## Total number of cells

ZF_LM_gaussian_summaries <- 
  panel_1M_LM_data %>%
  filter(Edge_case == "No") %>% ## Remove all cells at the boundaries of SF, TF
  #filter(`number of peaks` == 1) %>% ## Only primary peaks allowed
  filter(Fit_procedure == "unconstrained") %>% ## Remove any constrained fits
  filter(Bird == "ZF")
#length(unique(ZF_LM_gaussian_summaries$Ind)) ## number of individuals
#nrow(ZF_LM_gaussian_summaries) ## Total number of cells

PG_LM_gaussian_summaries <- 
  panel_1M_LM_data %>%
  filter(Edge_case == "No") %>% ## Remove all cells at the boundaries of SF, TF
  #filter(`number of peaks` == 1) %>% ## Only primary peaks allowed
  filter(Fit_procedure == "unconstrained") %>% ## Remove any constrained fits
  filter(Bird == "PG")
#length(unique(PG_LM_gaussian_summaries$Ind)) ## number of individuals
#nrow(PG_LM_gaussian_summaries) ## Total number of cells