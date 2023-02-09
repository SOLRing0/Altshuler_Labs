##### file lists #####
all_spatemp <- list.files(
  "./data/spatemp/LM_ER_fullsweep_200bins_bin001_to_bin200", 
  recursive = TRUE, full.names = TRUE
)

all_psth <-
  all_spatemp[str_which(all_spatemp, "-ALL_CELLS-SUBSET_PSTH_MEAN.csv")]
all_spont <-
  all_spatemp[str_which(all_spatemp, "-spontaneous_rate.csv")]

all_psth_hb <-
  all_psth[str_which(all_psth, "CALSTC0")]
all_spont_hb <-
  all_spont[1:length(all_psth_hb)]

all_psth_zf <-
  all_psth[(length(all_psth_hb)+1):length(all_psth)]
all_spont_zf <-
  all_spont[(length(all_psth_hb)+1):length(all_psth)]


##### panel 4Ci and 4Cii HB #####
case_num <- 1
data_4Ai_csv <-
  read_csv(all_psth_hb[case_num], col_names = FALSE
  ) %>% as.data.frame()
data_4Ai_spont <-
  read_csv(all_spont_hb[case_num], col_names = FALSE
  ) %>% as.data.frame()

cell_num <- 1 # cell number
nbins <- 200 # total number of bins

data_4Ai_plotdat <- data_4Ai_csv[,4:ncol(data_4Ai_csv)]
data_4Ai_max_init <-
  apply(data_4Ai_plotdat, 2, function(x)
    max(x, na.rm = TRUE))
data_4Ai_max_max <- max(data_4Ai_max_init)
data_4Ai_max_corrected <-
  as.numeric(data_4Ai_max_max - data_4Ai_spont[cell_num,])
data_4Ai_basearray <- c(rep(data_4Ai_spont[cell_num,], nbins))
data_4Ai_psth_letters <- seq(1:nbins)
data_4Ai_spontplot <- as.numeric(data_4Ai_spont[cell_num,])

data_4Ai_plotz <- NULL
data_4Ai_plt_subbase <- NULL
panel_4Ci_datlist <- NULL
for (i in 1:nrow(data_4Ai_plotdat)) {
  data_4Ai_plt_subbase[[i]] <-
    data_4Ai_plotdat[i, 1:nbins] #- data_4Ai_basearray
  tmp <- as.vector(t(data_4Ai_plt_subbase[[i]]))
  names(tmp) <- NULL
  panel_4Ci_datlist[[i]] <- tibble(lets = data_4Ai_psth_letters,
                                  psth_dat = tmp)
  data_4Ai_plotz[[i]] <-
    ggplot(panel_4Ci_datlist[[i]], aes(x = lets, y = psth_dat)) +
    geom_col(fill = col_hb,
             color = col_hb,
             size = 0) +
    #scale_y_continuous(breaks = seq(0, data_4Ai_max_corrected, by = 10)) +
    scale_x_continuous(breaks = c(0, 100, 200), 
                       labels = c("0", "1", "2")) +
    coord_cartesian(ylim =
                      c(0,#-1*data_4Ai_spontplot,
                        200 #data_4Ai_max_corrected
                        ),
                        expand = FALSE) +
    geom_hline(yintercept = data_4Ai_spontplot,
               linetype = "dotted",
               color = "grey20",
    ) +
    ylab("spikes/s") +
    xlab("time (s)") +
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
}

panel_4Ci_ii_psth <- plot_grid(plotlist = rev(data_4Ai_plotz),
                               nrow = 6)
panel_4Ci <- data_4Ai_plotz[[30]] ## strong initial transient
panel_4Cii <- data_4Ai_plotz[[4]] ## weak IT; consistent throughout
# case  1: 30 and 4
# case 14: 33 and 17

##### panel 4Ciii and 4Civ ZF #####
case_num <- 56
data_4Aiii_csv <-
  read_csv(all_psth_zf[case_num], col_names = FALSE
  ) %>% as.data.frame()
data_4Aiii_spont <-
  read_csv(all_spont_zf[case_num], col_names = FALSE
  ) %>% as.data.frame()

cell_num <- 1 # cell number
nbins <- 200 # total number of bins

data_4Aiii_plotdat <- data_4Aiii_csv[,4:ncol(data_4Aiii_csv)]
data_4Aiii_max_init <-
  apply(data_4Aiii_plotdat, 2, function(x)
    max(x, na.rm = TRUE))
data_4Aiii_max_max <- max(data_4Aiii_max_init)
data_4Aiii_max_corrected <-
  as.numeric(data_4Aiii_max_max - data_4Aiii_spont[cell_num,])
data_4Aiii_basearray <- c(rep(data_4Aiii_spont[cell_num,], nbins))
data_4Aiii_psth_letters <- seq(1:nbins)
data_4Aiii_spontplot <- as.numeric(data_4Aiii_spont[cell_num,])

data_4Aiii_plotz <- NULL
data_4Aiii_plt_subbase <- NULL
panel_4Ciii_datlist <- NULL
for (i in 1:nrow(data_4Aiii_plotdat)) {
  data_4Aiii_plt_subbase[[i]] <-
    data_4Aiii_plotdat[i, 1:nbins] #- data_4Aiii_basearray
  tmp <- as.vector(t(data_4Aiii_plt_subbase[[i]]))
  names(tmp) <- NULL
  panel_4Ciii_datlist[[i]] <- tibble(lets = data_4Aiii_psth_letters,
                                     psth_dat = tmp)
  data_4Aiii_plotz[[i]] <-
    ggplot(panel_4Ciii_datlist[[i]], aes(x = lets, y = psth_dat)) +
    # col_hb = "#ED0080"
    geom_col(fill = col_zf,
             color = col_zf,
             size = 0) +
    #scale_y_continuous(breaks = seq(0, data_4Aiii_max_corrected, by = 10)) +
    scale_x_continuous(breaks = c(0, 100, 200), 
                       labels = c("0", "1", "2")) +
    coord_cartesian(ylim = 
                      c(0, #-1*data_4Aiii_spontplot,
                        170 #data_4Aiii_max_corrected
                        ),
                    expand = FALSE) +
    geom_hline(yintercept = data_4Aiii_spontplot,
               linetype = "dotted",
               color = "grey20",
               ) +
    ylab("spikes/s") +
    xlab("time (s)") +
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
}

panel_4Ciii_iv_psth <- plot_grid(plotlist = rev(data_4Aiii_plotz),
                                 nrow = 6)
panel_4Ciii <- data_4Aiii_plotz[[32]] ## strong initial transient
panel_4Civ <- data_4Aiii_plotz[[13]] ## weak IT; consistent throughout
