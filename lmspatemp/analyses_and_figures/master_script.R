## Note from VBB
## This is the master script used to recreate all figures and analyses. Specific
## scripts pertaining to each figure and corresponding analyses can be found in
## each subfolder within this directory. The master script sources all specific
## scripts and then generates multi-panel plots for each figure. For each
## figure, some tuning of aesthetics may also have been done in Adobe
## Illustrator, none of which affected patterns in data. Insertions of 
## illustrations (e.g., Fig 1A, F, bird heads throughout all figures) were 
## handled in Illustrator.
## 
## It is strongly recommended to load all necessary packages (see next section) 
## prior to running anything else.


############################ Load necessary packages ###########################
## Specify the packages of interest
packages <-
  c("tidyverse", "ggplot2", "tidyr", "dplyr", "tibble", "readr", "readxl",
    "ggthemes", "magick", "grid", "cowplot", "simr", "MCMCglmm", "tidybayes", 
    "reshape2", "lme4", "splines", "circular", "metR", "viridis", "R.matlab",
    "gaussplotR", "plotly", "rstatix")

## Load or install & load all packages
package.check <- lapply(
  packages,
  FUN = function(x)
  {
    if (!require(x, character.only = TRUE))
    {
      install.packages(x, dependencies = TRUE,
                       repos = "http://cran.us.r-project.org")
      library(x, character.only = TRUE)
    }
  }
)

############################## Set up other generics ###########################
## For 2-column formats (such as research articles and reviews), the sizes are
## 85 mm (1 column), 114 mm (1.5 columns), and 174 mm (full width of the page)
## max figure widths, units in inches
fig_width_sm <- 3.34 #85 mm
fig_width_md <- 4.47 #114 mm
fig_width_lg <- 6.84 #174 mm

## expressions for plot tick labels
two_neg_six   <- expression(paste(2^-6)) 
two_neg_five  <- expression(paste(2^-5)) 
two_neg_four  <- expression(paste(2^-4)) 
two_neg_three <- expression(paste(2^-3)) 
two_neg_two   <- expression(paste(2^-2)) 
two_neg_one   <- expression(paste(2^-1)) 
two_zero      <- expression(paste(2^0)) 
two_six       <- expression(paste(2^6))
two_five      <- expression(paste(2^5))
two_four      <- expression(paste(2^4))
two_three     <- expression(paste(2^3))
two_two       <- expression(paste(2^2))
two_one       <- expression(paste(2^1))
ten_neg_two   <- expression(paste(10^-2)) 
ten_neg_one   <- expression(paste(10^-1)) 
ten_zero      <- expression(paste(10^0)) 
ten_three     <- expression(paste(10^3))
ten_two       <- expression(paste(10^2))
ten_one       <- expression(paste(10^1))
## r2 expression
r2 <- expression(R^2) #%>% as.character()


## Species-specific colors
col_hb <- "#ED0080"
# zb and zf are the same, just different abbreviations for finches that may
# occurr
col_zb <- "#F48D00"
col_zf <- "#F48D00" 
col_pg <- "#12bcff"

## other plot generics
plot_margs <- unit(c(0.1, 0.1, 0.1, 0.1), "cm")
contour_wds <- 0.04
panel_font_size = 12
geom_text_size = 1.5

################################### FIGURE 1 ###################################
## The following blocks suppressMessages() and suppressWarnings() to try to 
## minimize console printing, but there will be a lot of output regardless. 

## Import panels
suppressMessages(suppressWarnings(source(
  "./analyses_and_figures/figure 1/Panel_1_BCD.R"
)))
suppressMessages(suppressWarnings(source(
  "./analyses_and_figures/figure 1/Panel_1_E.R"
)))
suppressMessages(suppressWarnings(source(
  "./analyses_and_figures/figure 1/Panel_1_GHI.R"
)))
suppressMessages(suppressWarnings(source(
  "./analyses_and_figures/figure 1/Panel_1_JKL.R"
)))

suppressMessages(suppressWarnings(source(
  "./analyses_and_figures/figure 1/Panel_1_M.R"
)))

panel_1BCD <- plot_grid(NULL, NULL, NULL,
                        panel_1B, panel_1C, panel_1D,
                        NULL, NULL, NULL,
                        nrow = 3,
                        rel_heights = c(0.1, 1, 0.05),
                        labels = c("B", "C", "D",
                                   " ", " ", " ", 
                                   " ", " ", " "))

fig1_top_row <-
  plot_grid(NULL, panel_1BCD, panel_1E,
            nrow = 1,
            labels = c("A", " ", "E"),
            rel_widths = c(1, 3, 2))

panel_1Gh <- plot_grid(NULL, panel_1G, NULL,
                       ncol = 3,
                       rel_widths = c(0.1, 1, 0.01)
)
## Next vertically
panel_1Gm <- plot_grid(panel_1Gh, NULL,
                       ncol = 1,
                       rel_heights = c(1, 0.1)
)

## Add padding to each legend
panel_1H_legpad <- plot_grid(NULL, panel_1H_legend, NULL,
                             ncol = 1,
                             rel_heights = c(0.05, 0.80, 0.2))

panel_1I_legpad <- plot_grid(NULL, panel_1I_legend, NULL,
                             ncol = 1,
                             rel_heights = c(0.05, 0.80, 0.2))

fig1_mid_row <- plot_grid(panel_1Gm, 
                          panel_1H, NULL, panel_1H_legend, NULL,
                          panel_1I, NULL, panel_1I_legend, NULL,
                          nrow = 1,
                          labels = c("G", 
                                     "H", "", "", "", 
                                     "I", "", "", ""),
                          label_size = panel_font_size,
                          label_x = c(0, 
                                      -0.03, 0, 0, 0, 
                                      -0.035, 0, 0, 0),
                          label_fontfamily = "sans",
                          rel_widths = c(0.40, 
                                         0.26, 0.005, 0.03, 0.0225, 
                                         0.24, 0.005, 0.03, 0.0075),
                          scale = c(0.9, 
                                    1, 1, 1, 1,
                                    1, 1, 1, 1)
)

fig1_gaussian_plots <- plot_grid(panel_1JKL,
                                 NULL,
                                 ncol = 1,
                                 rel_heights = c(1, 0.05),
                                 labels = c("J"," ", " ", ""),
                                 label_x = -0.055,
                                 label_size = panel_font_size,
                                 label_fontfamily = "sans"
)

cow_1M <- plot_grid(
  panel_1M,
  ncol = 1,
  labels = c("M"),
  label_x = 0.21,
  label_size = panel_font_size,
  label_fontfamily = "sans"
)

fig1_bot_row <- plot_grid(NULL, fig1_gaussian_plots, NULL, cow_1M,
                         nrow = 1,
                         rel_widths = c(0.04, 0.76, 0.01, 0.25),
                         label_size = panel_font_size,
                         label_fontfamily = "sans"
)

Figure_1 <- 
  plot_grid(fig1_top_row,
            fig1_mid_row,
            fig1_bot_row,
            rel_heights = c(1, 1, 1),
            ncol = 1)

pdf(file = "./analyses_and_figures/Figure 1.pdf",
    width = fig_width_lg, height = fig_width_lg*0.8,
    title = "Figure 1", paper = "letter", bg = "white",
    pagecentre = TRUE, colormodel = "srgb")
Figure_1
dev.off()

################################## FIGURE S1 ###################################
suppressMessages(suppressWarnings(source(
  "./analyses_and_figures/figure S1/Panel_S1AB.R"
)))

Figure_S1 <-
  plot_grid(panel_S1A, panel_S1B,
            nrow = 1,
            labels = c("A", "B"),
            rel_widths = c(1.2, 0.75),
            label_size = panel_font_size,
            label_x = c(0.152, 0.228),
            label_fontfamily = "sans"
            )

pdf(file = "./analyses_and_figures/Figure S1.pdf",
    width = fig_width_md, height = fig_width_sm*0.6,
    title = "Figure S1", paper = "letter", bg = "white",
    pagecentre = TRUE, colormodel = "srgb")
Figure_S1
dev.off()

################################### FIGURE 2 ###################################
suppressMessages(suppressWarnings(source(
  "./analyses_and_figures/figure 2/Panel_2_A2D.R"
)))

cow_2A <-
  plot_grid(
    NULL,
    Panel_2A_complete,
    nrow = 1,
    rel_widths = c(0.05, 1),
    labels = c(" ", "A"),
    #label_y = 1.05,
    #rel_heights = c(0.05, 1),
    label_x = 0,
    label_size = panel_font_size,
    label_fontfamily = "sans"
  )

cow_2B2E <-
  plot_grid(
    Panel_2B_LM_vel, Panel_2C_LM_SF, Panel_2D_LM_TF,
    #rel_heights = c(1, -0.2, 1),
    labels = c("B", "C", "D"),
    #label_x = 0.18,
    label_size = panel_font_size,
    label_fontfamily = "sans",
    nrow = 1
  )

## stitch together
Figure_2 <-
  plot_grid(cow_2A, cow_2B2E,
            ncol = 1,
            rel_heights = c(2.1, 0.9)
            )
   
pdf(file = "./analyses_and_figures/Figure 2.pdf",
    width = fig_width_sm, height = fig_width_sm*1.05,
    title = "Figure 2", paper = "letter", bg = "white",
    pagecentre = TRUE, colormodel = "srgb")
Figure_2
dev.off()


################################### FIGURE 3 ###################################
suppressMessages(suppressWarnings(source(
  "./analyses_and_figures/figure 3/Panel_3_A2J.R"
)))
suppressMessages(suppressWarnings(source(
  "./analyses_and_figures/figure 3/Panel_3_KLM.R"
)))

cow_3A2E <-
  plot_grid(
    NULL, p_3A2E_i, p_3A2E_iii, p_3A2E_iv,NULL, panel_3D, NULL, panel_3E,
    nrow = 1,
    rel_widths = c(0.05, 0.75, 0.75, 0.75, 0.07, 1.45, -0.03, 0.95),
    labels = c("", "A", "B", "C", "", "D", "", "E"),
    label_x = c(0, -0.12, -0.01, -0.03, 0, 0.16, 0, 0.13),
    label_y = 1.05,
    label_size = panel_font_size,
    label_fontfamily = "sans"
  )

cow_3F2J <-
  plot_grid(
    NULL, p_3F2J_i, p_3F2J_iii, p_3F2J_iv,NULL, panel_3I, NULL, panel_3J,
    nrow = 1,
    rel_widths = c(0.05, 0.75, 0.75, 0.75, 0.07, 1.45, -0.03, 0.95),
    labels = c("", "F", "G", "H", "", "I", "", "J"),
    label_x = c(0, -0.12, -0.01, -0.03, 0, 0.16, 0, 0.13),
    label_y = 1.05,
    label_size = panel_font_size,
    label_fontfamily = "sans"
  )

cow_3LM <- 
  plot_grid(Panel_3L, NULL, Panel_3M, 
            nrow = 1,
            rel_widths = c(0.65, 0.35, 1), 
            labels = c("L", "" , "M"),
            label_size = panel_font_size,
            label_fontfamily = "sans")

cow_3LMs <-
  plot_grid(NULL, cow_3LM,
            nrow = 2, 
            rel_heights = c(0.2, 1))

cow_3KLM <-
  plot_grid(Panel_3K, cow_3LMs, 
            nrow = 1,
            rel_widths = c(1, 1.3),
            labels = c("K"," "),
            label_size = panel_font_size,
            label_fontfamily = "sans")

## stitch together
Figure_3 <-
  plot_grid(cow_3A2E, cow_3F2J, cow_3KLM,
            nrow = 3,
            rel_heights = c(1, 1, 1.75)
  )

pdf(file = "./analyses_and_figures/Figure 3.pdf",
    width = fig_width_lg, height = fig_width_lg*0.66,
    title = "Figure 3", paper = "letter", bg = "white",
    pagecentre = TRUE, colormodel = "srgb")
Figure_3
dev.off()


################################### FIGURE 4 ###################################
suppressMessages(suppressWarnings(source(
  "./analyses_and_figures/figure 4/Panel_4AB.R"
)))
suppressMessages(suppressWarnings(source(
  "./analyses_and_figures/figure 4/Panel_4C.R"
)))
suppressMessages(suppressWarnings(source(
  "./analyses_and_figures/figure 4/Panel_4D2F.R"
)))
suppressMessages(suppressWarnings(source(
  "./analyses_and_figures/figure 4/Panel_4G2L.R"
)))

panel_4AB <-
  plot_grid(panel_4A, NULL, panel_4B,
            labels = c("A", " ", "B"),
            rel_widths = c(1, 0.02, 1),
            nrow = 1)

panel_4C <-
  plot_grid(panel_4Ci,   panel_4Cii,
            panel_4Ciii, panel_4Civ,
            labels = c("C", c(rep("", 3))),
            ncol = 2)
            
panel_4D2F <-
  plot_grid(panel_4D, panel_4E, panel_4F,
            labels = c("D", "E", "F"),
            nrow = 1)

panel_4G <-
  plot_grid(NULL, p_4G_hb, NULL, p_4G_zf,
            nrow = 1,
            labels = c(" ","G", "", ""),
            rel_widths = c(0.05, 1, -0.05, 1),
            label_size = panel_font_size,
            label_fontfamily = "sans")
            
panel_4H_pans <-
  plot_grid(p_4H_vel_hb, p_4H_vel_zf,
            p_4H_SF_hb,  p_4H_SF_zf,
            p_4H_TF_hb,  p_4H_TF_zf,
            nrow = 1,
            labels = c("H", rep("", 5)),
            #rel_widths = c(1, 0.1, 1),
            label_size = panel_font_size,
            label_fontfamily = "sans")

panel_4H <-
  plot_grid(NULL, panel_4H_pans, NULL,
            ncol = 1,
            rel_heights = c(0.07, 1.12, 0.01))

panel_4I <-
  plot_grid(NULL, p_4I_hb, NULL, p_4I_zf,
            nrow = 1,
            labels = c(" ","I", "", ""),
            rel_widths = c(0.05, 1, -0.05, 1),
            label_size = panel_font_size,
            label_fontfamily = "sans")

panel_4J2L_pans <-
  plot_grid(p_4J_Rvel_hb, p_4J_Rvel_zf,
            p_4J_Rind_hb, p_4J_Rind_zf,
            Panel_4K, Panel_4L,
            nrow = 1,
            labels = c("J", rep("", 3),
                       "K", "L"),
            rel_widths = c(1, 1, 1, 1, 0.75, 0.75),
            label_size = panel_font_size,
            label_fontfamily = "sans")

panel_4J2L <-
  plot_grid(NULL, panel_4J2L_pans, NULL,
            ncol = 1,
            rel_heights = c(0.07, 1.06, 0.07))

panel_4CGI <-
  plot_grid(panel_4C, panel_4G, panel_4I,
            ncol = 1)
panel_4DHJ <-
  plot_grid(panel_4D2F, NULL, panel_4H, NULL, panel_4J2L,
            ncol = 1, 
            rel_heights = c(1, 0.02, 1, 0.02, 1)
            )

fig_4_bottom <-
  plot_grid(
    panel_4CGI,
    NULL,
    panel_4DHJ,
    ncol = 3,
    rel_widths = c(1.05, 0.05, 1.5),
    label_size = panel_font_size,
    label_fontfamily = "sans"
  )

Figure_4 <- 
  plot_grid(
    panel_4AB, 
    NULL,
    fig_4_bottom,
    ncol = 1,
    rel_heights = c(0.32, 0.04, 0.64)
  )

pdf(file = "./analyses_and_figures/Figure 4.pdf",
    width = fig_width_lg - 0.5, height = fig_width_lg + 0.12, 
    title = "Figure 4", paper = "letter", bg = "white",
    pagecentre = TRUE, colormodel = "srgb")
Figure_4
dev.off()

