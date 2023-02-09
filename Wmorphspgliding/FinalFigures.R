## ---------------------------------------------------------------------------------------
## -------------------------------- FIGURE 2: PCA FIGURE --------------------------------

### ---------- PANEL D PC1 vs PC2 ---------
plot_pca <- ggplot()+
  geom_point(data = subset(pc_results, dataset %in% "labman"), aes(x = PC1, y = PC2), col = col_lab, size = 0.8, alpha = 0.7) +
  labs(title="", y="", x="PC1") +
  geom_point(data = subset(pc_results, dataset %in% "invivo"), aes(x = PC1, y = PC2), 
             colour   = line_invivo, fill   = col_invivo, shape  = 24, size   = 5) +
  geom_point(data = subset(pc_results, dataset %in% "windtunnel" & orientation %in% "hand"), aes(x=PC1, y= PC2), 
             fill   = col_wt, colour = line_wt, shape  = 22, size = 5) +
  scale_colour_gradientn(name   = lab_elbow,
                         colours = cc, limits = c(30,150), breaks = c(45,65,85,105,125,145), labels =  c(45,65,85,105,125,145)) +
   th +
  scale_x_continuous(limits = pc1_limits, expand = c(0,0), breaks = pc1_breaks) +
  scale_y_continuous(limits = pc2_limits, expand = c(0,0), breaks = pc2_breaks) +
  theme(legend.position = "none",
        axis.text.y = element_blank())

### ---------- PANEL A ELBOW ANGLE VS MANUS ANGLE ---------
plot_angles <- ggplot(subset(pc_results, dataset %in% "labman"), 
                       aes(y = manus.angle, x = elbow.angle))+
  geom_point(pch = 16, col = col_lab, size = 0.8, alpha = 0.7)+
  geom_point(data = pc_results_invivo, aes(y=manus.angle, x=elbow.angle, fill = Wind.Speed), 
             colour = line_invivo, shape  = 24, size   = 5, alpha = 0.7) + 
  geom_point(data = subset(pc_results, dataset %in% "windtunnel" & orientation %in% "hand"), aes(y=manus.angle, x=elbow.angle), 
             fill   = col_wt, colour   = line_wt, shape  = 22, size   = 5) +
  labs(title="", x="",y=lab_manus) +
  scale_x_continuous(limits = elbow_limits, expand = c(0,0), breaks = elbow_breaks) +
  scale_y_continuous(limits = manus_limits, expand = c(0,0), breaks = manus_breaks) + 
  th + theme(legend.position = "none",
             axis.text.x = element_blank())

### ---------- PANEL C ELBOW ANGLE VS PC2 ---------
plot_elbow_den_invivo <- ggplot()  +
  geom_point(data = subset(pc_results, dataset %in% "labman" ), aes(x = PC2, y = elbow.angle),col = col_lab, size = 0.8)+
  geom_point(data = pc_results_invivo, aes(x = PC2, y = elbow.angle),
             fill   = col_invivo, colour = line_invivo, shape  = 24, size   = 5, alpha = 0.7)+
  geom_point(data = subset(pc_results, dataset %in% "windtunnel" & orientation %in% "hand"), aes(x=PC2, y=elbow.angle), 
             fill   = col_wt, colour   = line_wt, shape  = 22, size   = 5) +
  labs(title="", y=lab_elbow, x="PC2") +
  scale_x_continuous(limits = pc2_limits, expand = c(0,0), breaks = pc2_breaks) + 
  scale_y_continuous(limits = elbow_limits, expand = c(0,0), breaks = elbow_breaks) +
  coord_flip() + th + 
  theme(legend.position = "none") 

### ---------- PANEL B MANUS ANGLE VS PC1 ---------
plot_manus_den_invivo <- ggplot() +
  geom_point(data = subset(pc_results, dataset %in% "labman" ), aes(x = PC1, y = manus.angle), col = col_lab, size = 0.8)+
  geom_point(data = pc_results_invivo, aes(x = PC1, y = manus.angle),
             fill   = col_invivo, colour = line_invivo,shape  = 24, size   = 5, alpha = 0.7)+
  geom_point(data = subset(pc_results, dataset %in% "windtunnel" & orientation %in% "hand"), aes(x=PC1, y=manus.angle), 
             fill   = col_wt, colour   = line_wt, shape  = 22, size   = 5) +
  labs(title="", y="", x="") +
  scale_x_continuous(limits = pc1_limits, expand = c(0,0), breaks = pc1_breaks) + 
  scale_y_continuous(limits = manus_limits, expand = c(0,0), breaks = manus_breaks) +
  th + theme(legend.position = "none",
                 axis.text.x = element_blank(),
                 axis.text.y = element_blank(),
                 panel.background = element_rect(fill = "transparent"), # bg of the panel
                 plot.background = element_rect(fill = "transparent")) # bg of the plot) 

## ------------------ FINAL PLOTTING -------------------

g_manus_den_invivo <- ggplot_gtable(ggplot_build(plot_manus_den_invivo))
g_elbow_den_invivo <- ggplot_gtable(ggplot_build(plot_elbow_den_invivo))
g_pca       <- ggplot_gtable(ggplot_build(plot_pca))
g_angles      <- ggplot_gtable(ggplot_build(plot_angles))
lay_pca <- rbind(c(NA, 6, 6,8),
                 c(1 , 5, 5, 7),
                 c(2 , 5, 5, 7),
                 c(NA, 3, 4,NA))

grid.arrange(plot_PC2max, plot_PC2min,
             plot_PC1min, plot_PC1max,
             g_pca,g_manus_den_invivo, g_elbow_den_invivo, g_angles,
             layout_matrix = lay_pca, widths=c(1,1,1,2), heights=c(2,1,1,1))

library(cowplot)
plot_grid(plot_angles,plot_manus_den_invivo,plot_elbow_den_invivo,plot_pca, nrow = 2, ncol= 2,align = "hv",rel_widths = c(1, 1), rel_heights = c(0.66,0.76))
#save with 12x8.52


## ---------------------------------------------------------------------------------------
## -------------------------------- FIGURE 3: EFFICIENCY FIGURE --------------------------------

## ---------------------------------------------------------------------------------------
## ------------------------------------------AERO EFFICIENCY-----------------------------

xlim1=c(0,1.1); 
ylim1=c(-0.5,1.25);

## ----- CREATES THE LEGEND AS A SEPERATE ITEM
plotlegend <- ggplot(data = subset(dat[order(dat$WingID, dat$Angle.No),],Grid %in% grids[1]),
            aes(x=CD, y=CLift, colour = Angle.True)) + 
  geom_path(size = 0.5) +           #add y values
  labs(y=lab_lift_multi[1], x=lab_drag) +
  scale_colour_gradientn(name   = lab_elbow,
                         colours = cc, limits = elbow_x, breaks = elbow_breaks, labels =  elbow_breaks) +
  coord_fixed(xlim=xlim1,ylim=ylim1) + th +
  theme(legend.text=element_text(size=8), legend.title = element_text(size=10))
elbowlegend <- g_legend(plotlegend)

##LOOPS - Generates graphs for Supplemental figure 6
for (plotcount in 1:3){
  p <- ggplot(data = subset(dat[order(dat$WingID, dat$Angle.No),],Grid %in% grids[plotcount]),
              aes(x=CD, y=CLift, colour = Angle.True)) + 
    geom_point(shape = 18, size = 0.8) +
    geom_path(size = 0.45, alpha = 0.9) +           #add y values
     labs(title=title_ti[plotcount], y=lab_lift_multi[plotcount], x=lab_drag) +
     scale_colour_gradientn(colours = cc, limits = elbow_x, breaks = elbow_breaks, labels =  elbow_breaks) +
     coord_fixed(xlim=xlim1,ylim=ylim1)+
    th_noleg + 
    theme(axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          plot.margin= unit(c(1, 1, 1, -0.7), "lines")) +
    scale_x_continuous(limits = xlim1, expand = c(0,0), breaks = c(0.0,0.5,1.0)) +
    scale_y_continuous(limits = ylim1, expand = c(0,0), breaks = c(-0.50,-0.25,0.0,0.25,0.50,0.75, 1.00,1.25))
  plot_list[[plotcount]] = p
}
plot_list[[2]] <- plot_list[[2]] + theme(axis.text.y = element_blank())
plot_list[[3]] <- plot_list[[3]] + theme(axis.text.y = element_blank())
#SAVE AS GROBS

g_LD1     <- ggplotGrob(plot_list[[1]]) #PANEL A
g_LD2     <- ggplotGrob(plot_list[[2]])
g_LD3     <- ggplotGrob(plot_list[[3]])


## ---------------------------------------------------------------------------------------
## -------------------------------- LONGITUDINAL STABILITY  --------------------------------

xlim1=c(-0.5,1.27)
ylim1=c(-0.65,0.07)
##LOOPS - Generates graphs for Supplemental figure 6
for (plotcount in 1:3){
  p <- ggplot(subset(stall_subset[order(stall_subset$WingID, stall_subset$Angle.No),], Grid %in% grids[plotcount]), 
              aes(x=CLift, y=Cm, colour = Angle.True)) + 
    geom_point(shape = 18, size = 0.8)+
    geom_path(size = 0.45, alpha = 0.9) +     
    labs(y=lab_pitch_multi[plotcount], x=lab_lift) +
    scale_colour_gradientn(colours = cc, limits = elbow_x) +
    coord_fixed(xlim=xlim1,ylim=ylim1) + th_noleg +
    theme(axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          plot.margin= unit(c(1, -0.7, 1, 1), "lines")) +
  scale_x_continuous(limits = xlim1, expand = c(0,0), breaks = c(-0.5,-0.25,0.0,0.25,0.5,0.75,1.0,1.25)) +
    scale_y_continuous(limits = ylim1, expand = c(0,0), breaks = c(-0.8,-0.6,-0.4,-0.2,0.0))
  plot_list[[plotcount]] = p
}

plot_list[[2]] <- plot_list[[2]] + theme(axis.text.y = element_blank())
plot_list[[3]] <- plot_list[[3]] + theme(axis.text.y = element_blank())
#SAVE AS GROBS
g_PS1     <- ggplotGrob(plot_list[[1]]) #PANEL E
g_PS2     <- ggplotGrob(plot_list[[2]])
g_PS3     <- ggplotGrob(plot_list[[3]])



## ---------------------------------------------------------------------------------------
## ------------------------------------------DERIVED EFFICIENCY RESULTS-----------------------------

#------ PANEL B + C : Maximum Lift Coefficient and Minimum Drag Coefficient ------
plot_Clmax <- ggplot(max_results) + 
  geom_point(aes(x=Angle.True, y=CLmax, fill = as.factor(Grid), colour = as.factor(Grid)), 
             shape =21, size = 4, alpha = 0.9) +#add y values 
  geom_line(aes(x=Angle.True, y = median_clmax, 
                colour = as.factor(Grid)), alpha = 0.7, linetype = 2, lwd = 1)+
  geom_errorbar(aes(x=Angle.True, ymin = CLmax-CLmaxsd, ymax = CLmax+CLmaxsd, colour = as.factor(Grid)),width = 3,  alpha = 0.3) +
  geom_point(aes(x=Angle.True, y=CDmin, fill = as.factor(Grid), colour = as.factor(Grid)), 
             shape =21, size = 4, alpha = 0.9) +#add y values 
  geom_line(aes(x=Angle.True, y = median_cdmin, 
                colour = as.factor(Grid)), alpha = 0.7, linetype = 2, lwd = 1)+
  geom_errorbar(aes(x=Angle.True, ymin = CDmin-CDminsd, ymax = CDmin+CDminsd, colour = as.factor(Grid)),width = 3,  alpha = 0.3) +
  labs(y=lab_CL_max, x="") +
  scale_fill_manual(name   = paste("Turbulence \n Intensity"),
                      breaks = c("NoGrid","Rd38","Large"),
                      labels = title_ti,
                      values = TIcolors) +
  scale_colour_manual(name   = paste("Turbulence \n Intensity"),
                    breaks = c("NoGrid","Rd38","Large"),
                    labels = title_ti,
                    values = TIlinecolors) +
  coord_cartesian(xlim=elbow_x,ylim=c(0.085,1.4)) + th +
  scale_x_continuous(limits = elbow_x, breaks = elbow_breaks) +
  scale_y_continuous(limits = c(0.085,1.4), breaks = c(0.09,0.135,0.18,0.85,1.05,1.25)) +
  theme(legend.text=element_text(size=8), legend.title = element_text(size=10), legend.title.align = 0.5) 


TIlegend <- g_legend(plot_Clmax)
plot_Clmax <- plot_Clmax + theme(legend.position = "none", 
                                 axis.text.x = element_blank(),
                                 plot.margin= unit(c(1, 1, -0.7, 0.5), "lines"))


#------ PANEL D: Maximum Lift/Drag Ratio ------
plot_LD <- ggplot(max_results) + 
    geom_point(aes(x=Angle.True, y=MaxL_D, fill = as.factor(Grid), colour = as.factor(Grid)), 
               shape =21, size = 4, alpha = 0.9) +
  geom_line(aes(x=Angle.True, y = median_ld, 
            colour = as.factor(Grid)), alpha = 0.7, linetype = 2, lwd = 1)+
  geom_errorbar(aes(x=Angle.True, ymin = MaxL_D-MaxL_Dsd, ymax = MaxL_D+MaxL_Dsd, colour = as.factor(Grid)),
                width = 3,  alpha = 0.3) +
    labs(y=lab_LD_max, x="") +
  scale_fill_manual(name   = "Turbulence Intensity",
                    breaks = c("NoGrid","Rd38","Large"),
                    labels = title_ti,
                    values = TIcolors) +
  scale_colour_manual(name   = "Turbulence Intensity",
                      breaks = c("NoGrid","Rd38","Large"),
                      labels = title_ti,
                      values = TIlinecolors) +
    coord_cartesian(xlim=elbow_x,ylim=c(2,5)) + th_noleg +
    scale_x_continuous(limits = elbow_x, breaks = elbow_breaks) +
    theme(axis.title.x = element_text(),
          axis.text.x = element_blank(),
          plot.margin= unit(c(1, 1, -0.7, 0.5), "lines"))

#------ PANEL G: Pitching Slope ------
ylim1=c(-0.9,0)
plot_cmcl <- ggplot(fit_results) + 
  geom_point(aes(x=Angle.True, y=slopemean, fill = as.factor(Grid), colour = as.factor(Grid)), 
             shape =21, size = 4, alpha = 0.9) + 
  geom_line(aes(x=Angle.True, y = median_cmcl, 
                colour = as.factor(Grid)), alpha = 0.7, linetype = 2, lwd = 1)+
  geom_errorbar(aes(x=Angle.True, ymin = slope2.5, ymax = slope97.5, colour = as.factor(Grid)),
                width = 3,  alpha = 0.3) +
  labs(y=lab_cm, x=lab_elbow) +
  scale_fill_manual(name   = lab_ti,
                      breaks = c("NoGrid","Rd38","Large"),
                      labels = title_ti,
                      values = TIcolors) +
  scale_colour_manual(name   = lab_ti,
                    breaks = c("NoGrid","Rd38","Large"),
                    labels = title_ti,
                    values = TIlinecolors) +
  coord_cartesian(xlim=elbow_x,ylim=ylim1) + th_noleg +
  scale_x_continuous(limits = elbow_x, breaks = elbow_breaks) +
  scale_y_continuous(breaks = c(-0.8,-0.6,-0.4,-0.2,0.0)) +
  theme(plot.margin= unit(c(1, 1, -0.7, 0.5), "lines"))

#------ PANEL F: Zero-Lift Pitching Moment ------
plot_cm0 <- ggplot(fit_results) + 
  geom_point(aes(x=Angle.True, y=Cm_AC, fill = as.factor(Grid), colour = as.factor(Grid)), 
             shape =21, size = 4, alpha = 0.9) + 
  geom_line(aes(x=Angle.True, y = median_cm0, 
                colour = as.factor(Grid)), alpha = 0.7, linetype = 2, lwd = 1)+
  geom_errorbar(aes(x=Angle.True, ymin = Cm_AC-Cm_AC_sd, ymax = Cm_AC+Cm_AC_sd, colour = as.factor(Grid)),
                width = 3,  alpha = 0.3) +
  labs(y=lab_cm0, x=lab_elbow) +
  scale_fill_manual(name   = lab_ti,
                      breaks = c("NoGrid","Rd38","Large"),
                      labels = title_ti,
                      values = TIcolors) +
  scale_colour_manual(name   = lab_ti,
                    breaks = c("NoGrid","Rd38","Large"),
                    labels = title_ti,
                    values = TIlinecolors) +
  coord_cartesian(xlim=elbow_x,ylim=c(-0.075,1.24)) + th_noleg +
  scale_x_continuous(limits = elbow_x, breaks = elbow_breaks) +
  scale_y_continuous(limits = c(-0.075,1.24), breaks = c(-0.067,-0.02,0.027)) +
  theme(plot.margin= unit(c(0, 1, -0.7, 0.5), "lines"))

#SAVE AS GROBS
g_Clmax    <- ggplotGrob(plot_Clmax)
g_LD       <- ggplotGrob(plot_LD)
g_cmcl     <- ggplotGrob(plot_cmcl)
g_cm0      <- ggplotGrob(plot_cm0)

## ---------------------------------------------------------------------------------------
## -------------------------------- ASSEMBLE FIGURE 3  --------------------------------

#paste graphs together
g_derived   = gridExtra::rbind.gtable(g_Clmax, g_cm0,size = "max")
g_derslope  = gridExtra::rbind.gtable(g_LD,    g_cmcl, size = "max")
g_LDplots   = gridExtra::cbind.gtable(g_LD1,   g_LD2, g_LD3,size = "max")
g_PSplots   = gridExtra::cbind.gtable(g_PS1,   g_PS2, g_PS3,size = "max")
g_allLD     = gridExtra::rbind.gtable(g_LDplots,g_PSplots,size = "max")

lay_fig4 <- rbind(c(1,NA,2,3,NA),
                  c(1,5 ,2,3,4 ),
                  c(1,NA,2,3,NA))

grid.arrange(g_allLD, g_derived,g_derslope,
             TIlegend,elbowlegend,
             layout_matrix = lay_fig4, widths = c(4,0.5,1.5,1.5,0.5))



## ---------------------------------------------------------------------------------------
## -------------------------------- Figure 4: Wing Camber --------------------------------
## ---------------------------------------------------------------------------------------

# - Spanwise camber - PANEL A
# This is just for one of the wind tunnel wings to build up figure 4 this was run for each wing and assembled in Adobe Illustrator
camber <- camber[order(camber$ID,camber$path.no),]
camber_plot <- ggplot()+
  geom_path(data = camber[camber$ID == Wing_Number[11],],aes(x = x_final, y = y_final, col = Angle.True)) + 
  geom_point(data = dat[dat$WingID == Wing_Number[11],],aes(x = -as.numeric(COG.span), y = as.numeric(COG.z), col = Angle.True), shape = 10, size = 9)+ 
  scale_colour_gradientn(colours = cc, limits = c(30,150))  + coord_fixed(xlim = c(-0.62,0), ylim = c(-0.15,0.09)) + th

max_results$L_Dper = ((max_results$MaxL_D - min(max_results$MaxL_D))/(max(max_results$MaxL_D) - min(max_results$MaxL_D)))*100
fit_results$cmclper = ((fit_results$slopemean - max(fit_results$slopemean))/(min(fit_results$slopemean) - max(fit_results$slopemean)))*100
                                                                     ##--- Create Gradients -----
# - Aerodynamic efficiency - PANEL B
lifteff <- ggplot() +
  geom_point(data = max_results, aes(y = Grid, x = Angle.True, size = L_Dper, fill = as.factor(Grid), col = as.factor(Grid)), shape = 21, 
             alpha = 0.7) + 
  scale_fill_manual(name   = lab_ti,
                     breaks = c("NoGrid","Rd38","Large"),
                     labels = title_ti,
                     values = TIcolors) +
  scale_colour_manual(name   = lab_ti,
                    breaks = c("NoGrid","Rd38","Large"),
                    labels = title_ti,
                    values = TIlinecolors) +
  th_noleg + coord_flip() + scale_radius(range = c(2,36)) +
                scale_x_continuous(limits = c(30,150), breaks = round(unique(max_results$Angle.True))) 

# - Passive pitch stability - PANEL C
longstab <- ggplot() +
  geom_point(data = fit_results, aes(y=Grid, x = Angle.True, size = cmclper, fill = as.factor(Grid), col = as.factor(Grid)), shape = 21, 
             alpha = 0.7) + 
  scale_fill_manual(name   = lab_ti,
                    breaks = c("NoGrid","Rd38","Large"),
                    labels = title_ti,
                    values = TIcolors) +
  scale_colour_manual(name   = lab_ti,
                      breaks = c("NoGrid","Rd38","Large"),
                      labels = title_ti,
                      values = TIlinecolors) +
  th_noleg + coord_flip() + scale_radius(range = c(2.5,36)) +
  scale_x_continuous(limits = c(30,150), breaks = round(unique(max_results$Angle.True))) 

# - Behaviour - PANEL D
wspeed <- ggplot(pc_results_invivo) +
  geom_point(aes(x = Wind.Speed, y = elbow.angle, col = elbow.angle)) + 
  scale_colour_gradientn(colours = cc_wind, limits = c(80,170))  + th_noleg +
  scale_y_continuous(limits = c(89,155), breaks =  c(round(unique(max_results$Angle.True)))) +
  scale_x_continuous(limits = c(2,26), breaks = c(2,6,10,14,18,22,26)) 

wgust <- ggplot(pc_results_invivo) +
  geom_point(aes(x = Max.Wind.Gust, y = elbow.angle, col = elbow.angle)) + 
  scale_colour_gradientn(colours = cc_wind, limits = c(80,170))  + th_noleg +
  scale_y_continuous(limits = c(89,155), breaks = (round(unique(max_results$Angle.True)))) +
  scale_x_continuous(limits = c(2,26), breaks = c(2,6,10,14,18,22,26)) 

wgust  <- ggplotGrob(wgust)
wspeed <- ggplotGrob(wspeed)
g_wind   = gridExtra::cbind.gtable(wspeed, wgust, size = "max")

## ---------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------
## ----------------------------- Supplemental Figures ------------------------------------
## ---------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------

## ---------------------------------------------------------------------------------------
## -------------------------------- Figure S1 --------------------------------------------
## ---------------------------------------------------------------------------------------
beta_plot <- ggplot(curve_results[order(curve_results$Angle.True),], aes(x = Angle.True, y = beta, col = Angle.True)) +
             geom_smooth(formula = y~poly(x,2), method = "lm", linetype = 2, col = "black") +
             geom_point()+
             scale_x_continuous(limits = elbow_x, breaks = elbow_breaks) +
             labs(x = lab_elbow, y = "Maximum span-wise camber (% of span)") +
             scale_colour_gradientn(colours = cc, limits = c(30,150))  + th


## ---------------------------------------------------------------------------------------
## -------------------------------- Figure S2 --------------------------------------------
## ---------------------------------------------------------------------------------------

errmodelplot <- ggplot(modelerr) +
  geom_point(aes(x = ID, y = rmserr_lab, group = 1, shape = as.factor(modeltype)), col = "black")+
  geom_point(aes(x = ID, y = rmserr_wt, group = 1, shape = as.factor(modeltype)), col = "forestgreen")+
  geom_point(aes(x = ID, y = rmserr_wtvivo, group = 1, shape = as.factor(modeltype)), col = "blue")+
  th_noleg + labs(title="", x="Model Identifier", y="RMS error (°)") +
  scale_y_continuous(limits = c(0,35)) + scale_x_discrete(limits=modelerr$ID[order(-modelerr$score)])

elbrangeplot <- ggplot(elbow.range.long) +
  geom_point(aes(x = ID, y = elbow.angle, col = elbow.angle, shape = as.factor(modeltype)), alpha = 0.3, position = position_jitter(w = 0.2)) +
  scale_colour_gradientn(colours = cc_wind, limits = c(70,170)) +
  scale_y_continuous(limits = c(70,180), breaks = c(70,80,90,100,110,120,130,140,150,160,170,180)) + 
  scale_x_discrete(limits=modelerr$ID[order(-modelerr$score)]) +
  th_noleg + labs(title="", x="Model Identifier", y="Predicted Elbow Angle (°)")

coeffplot <- ggplot(modelerr) +
  geom_point(aes(x = ID, y = coeff.wind, shape = as.factor(modeltype)), col = "darkblue") +
  geom_point(aes(x = ID, y = coeff.gust, shape = as.factor(modeltype)), col = "lightblue") +
  geom_errorbar(aes(x=ID, ymin = confint2.5.wind, ymax = confint97.5.wind), 
                col = "darkblue", width = 0.2, alpha = 0.8)+
  geom_errorbar(aes(x=ID, ymin = confint2.5.gust, ymax = confint97.5.gust), 
                col = "lightblue", width = 0.2, alpha = 0.8) +
  scale_x_discrete(limits=modelerr$ID[order(-modelerr$score)]) +
  scale_y_continuous(limits = c(-1.7,0.3)) + 
  th_noleg + labs(title="", x="Model Identifier", y="Slope of Regression (°/m/s)")

g_errmodelplot   <- ggplotGrob(errmodelplot)
g_elbrangeplot   <- ggplotGrob(elbrangeplot)
g_coeffplot     <- ggplotGrob(coeffplot)
g_4  = gridExtra::rbind.gtable(g_errmodelplot,g_elbrangeplot,g_coeffplot,size = "max")


## ---------------------------------------------------------------------------------------
## -------------------------------- Figure S3 --------------------------------------------
## ---------------------------------------------------------------------------------------
predict.elb1 <-  predict(rf.4pc,    newdata = pc_results_wt[,7:10])
predict.elb2 <-  predict(rf.4pc,    newdata = pc_results_labman[,7:10]) 

err1 <- abs(predict.elb1 - pc_results_wt$elbow.angle)
err2 <- abs(predict.elb2 - pc_results_labman$elbow.angle)

err.elb <- ggplot()+
  geom_point(aes(x = pc_results_wt$elbow.angle, y = err1, col = pc_results_wt$elbow.angle)) + 
  labs(y = "Absolute error of predicted elbow angle (°)", x = "True elbow angle (°)")+
  coord_cartesian(xlim=elbow_limits,ylim=c(0,80)) + 
  scale_colour_gradientn(colours = cc, limits = c(30,150))  +
  scale_y_continuous(breaks = c(0,10,20,30), expand = c(0,0))+ 
  scale_x_continuous(breaks = elbow_breaks, expand = c(0,0)) + th_noleg

predict.elb1.man <-  predict(rf.4pc.man,    newdata = pc_results_wt[,7:10])
predict.elb2.man <-  predict(rf.4pc.man,    newdata = pc_results_labman[,7:10])

err1.man <- abs(predict.elb1.man - pc_results_wt$manus.angle)
err2.man <- abs(predict.elb2.man - pc_results_labman$manus.angle)


err.man <- ggplot()+
  geom_point(aes(x = pc_results_wt$manus.angle, y = err1.man, col = pc_results_wt$elbow.angle)) +
  labs(y = "Absolute error of predicted manus angle (°)", x = "True manus angle (°)")+
  coord_fixed(xlim=c(120,170),ylim=c(0,80)) + 
  scale_colour_gradientn(colours = cc, limits = c(30,150))  +
  scale_y_continuous(breaks = c(0,10,20,30), expand = c(0,0)) + 
  scale_x_continuous(breaks = c(120,130,140,150,160,170), expand = c(0,0)) +th_noleg

g_err.elb   <- ggplotGrob(err.elb)
g_err.man    <- ggplotGrob(err.man)

g_1  = gridExtra::cbind.gtable(g_err.elb,g_err.man,size = "max")

sens.elb.x <- ggplot(plotdf.y) +
  geom_point(aes(x = rot.x, y = err.elb.abs, col = true.elbow.angle)) + 
  labs(y = "Absolute error of predicted elbow angle (°)", x = "Rotation about the body axis (°)")+ th_noleg +
  coord_cartesian(ylim=c(0,80), xlim = c(-85,25)) + 
  scale_colour_gradientn(colours = cc, limits = c(30,150))  +
  scale_y_continuous(breaks = c(0,20,40,60,80), expand = c(0,0)) +
  scale_x_continuous(breaks = c(-80,-60,-40,-20,0,20), expand = c(0,0)) 

sens.elb.y <- ggplot(plotdf.x) +
  geom_point(aes(x = rot.y, y = err.elb.abs, col = true.elbow.angle)) + 
  labs(y = "Absolute error of predicted elbow angle (°)", x = "Rotation about the span axis (°)")+ th_noleg +
  coord_cartesian(ylim=c(0,80), xlim = c(-75,25)) + 
  scale_colour_gradientn(colours = cc, limits = c(30,150))  +
  scale_y_continuous(breaks = c(0,20,40,60,80), expand = c(0,0)) +
  scale_x_continuous(breaks = c(-75,-50,-25,0,25), expand = c(0,0))  

g_sens.elb.x   <- ggplotGrob(sens.elb.x)
g_sens.elb.y    <- ggplotGrob(sens.elb.y)

g_2  = gridExtra::cbind.gtable(g_sens.elb.x,g_sens.elb.y,size = "max")

sens.man.x <- ggplot(plotdf.y) +
  geom_point(aes(x = rot.x, y = err.man.abs, col = true.elbow.angle)) + 
  labs(y = "Absolute error of predicted manus angle (°)", x = "Rotation about the body axis (°)")+ th_noleg +
  scale_colour_gradientn(colours = cc, limits = c(30,150))  +
  coord_cartesian(ylim=c(0,80), xlim = c(-85,25)) + 
  scale_y_continuous(breaks = c(0,10,20), expand = c(0,0)) +
  scale_x_continuous(breaks = c(-80,-60,-40,-20,0,20), expand = c(0,0)) 

sens.man.y <- ggplot(plotdf.x) +
  geom_point(aes(x = rot.y, y = err.man.abs, col = true.elbow.angle)) + 
  labs(y = "Absolute error of predicted manus angle (°)", x = "Rotation about the span axis (°)")+ th_noleg +
  scale_colour_gradientn(colours = cc, limits = c(30,150))  +
  coord_cartesian(ylim=c(0,80), xlim = c(-75,25)) + 
  scale_y_continuous(breaks = c(0,10,20), expand = c(0,0)) +
  scale_x_continuous(breaks = c(-75,-50,-25,0,25), expand = c(0,0))  


g_sens.man.x   <- ggplotGrob(sens.man.x)
g_sens.man.y    <- ggplotGrob(sens.man.y)

g_3  = gridExtra::cbind.gtable(g_sens.man.x,g_sens.man.y,size = "max")

## ---------------------------------------------------------------------------------------
## -------------------------------- Figure S4 --------------------------------------------
## ---------------------------------------------------------------------------------------

chordplot = ggplot(dat[!duplicated(dat$Angle.True),]) +
  geom_point(aes(x = Angle.True, y= Root.Chord, col = Angle.True)) + 
  labs(y = "Root Chord (m)", x = "") + th_noleg +
  scale_colour_gradientn(colours = cc, limits = c(30,150)) +
  scale_y_continuous(limits = c(0.16,0.26), breaks = c(0.16,0.19,0.22,0.25), expand = c(0,0)) + 
  theme(axis.text.x = element_blank())

areaplot = ggplot(dat[!duplicated(dat$Angle.True),]) +
  geom_point(aes(x = Angle.True, y= Wing.Area, col = Angle.True)) + 
  labs(y = "Projected Area (m )", x = "") + th_noleg +
  scale_colour_gradientn(colours = cc, limits = c(30,150)) +
  scale_y_continuous(limits = c(0.06,0.1), breaks = c(0.06,0.08,0.1), expand = c(0,0)) +
  theme(axis.text.x = element_blank())


Spanplot = ggplot(dat[!duplicated(dat$Angle.True),]) +
  geom_point(aes(x = Angle.True, y= Wing.Eff.Span, col = Angle.True)) + 
  labs(y = "Effective Wing Span (m)", x = lab_elbow) + th_noleg +
  scale_colour_gradientn(colours = cc, limits = c(30,150)) +
  scale_y_continuous(limits = c(0.42,0.62), breaks = c(0.42,0.52,0.62), expand = c(0,0)) 


ARplot = ggplot(dat[!duplicated(dat$Angle.True),]) +
  geom_point(aes(x = Angle.True, y= AR, col = Angle.True)) + 
  labs(y = "Aspect Ratio", x = lab_elbow) + th_noleg +
  scale_colour_gradientn(colours = cc, limits = c(30,150))  +
  scale_y_continuous(limits = c(5.4,8.5), breaks = c(5.5,6.5,7.5,8.5), expand = c(0,0)) 

lay_fig4 <- rbind(c(1,2))
g_chordplot           <- ggplotGrob(chordplot)
g_areaplot    <- ggplotGrob(areaplot)
g_ARplot        <- ggplotGrob(ARplot)
g_Spanplot        <- ggplotGrob(Spanplot)

g_1  = gridExtra::rbind.gtable(g_chordplot,g_Spanplot,size = "max")
g_2  = gridExtra::rbind.gtable(g_areaplot,g_ARplot,size = "max")
grid.arrange(g_1,g_2, layout_matrix = lay_fig4)

## ---------------------------------------------------------------------------------------
## -------------------------------- Figure S5 --------------------------------------------
## ---------------------------------------------------------------------------------------
indplot = ggplot(data = max_results)+
  geom_point(aes(x = Angle.True, y = MaxL_D)) +
  geom_point(aes(x = Angle.True, y = MaxL_Dtrue), col = "red") + 
  scale_x_continuous(limits = c(25,160), breaks = elbow_breaks) +
  labs(x = lab_elbow, y = lab_LD) + th

errindplot = ggplot(data = max_results)+
  geom_point(aes(x = Angle.True, y = (MaxL_Dtrue-MaxL_D)*100/MaxL_Dtrue)) +
  scale_x_continuous(limits = c(25,160), breaks = elbow_breaks) +
  labs(x = "", y = "Error (%)") + th +
  theme(axis.text.x = element_blank())

sideoffsetplot = ggplot(data = dat)+
  geom_point(aes(x = Angle.True, y = (side_offset - ((3+(1/16))*25.4)/1000))) +
  scale_x_continuous(limits = c(25,160), breaks = elbow_breaks) +
  scale_y_continuous(breaks = c(0.06,0.09,0.12,0.15))+
  labs(x = lab_elbow, y = "Distance from Body (m)") + th

lay_fig4 <- rbind(c(1,2),
                  c(1,2))
g_indplot           <- ggplotGrob(indplot)
g_sideoffsetplot    <- ggplotGrob(sideoffsetplot)
g_errindplot        <- ggplotGrob(errindplot)

g_supp  = gridExtra::rbind.gtable(g_errindplot,g_sideoffsetplot,size = "max")
grid.arrange(g_indplot,g_supp, layout_matrix = lay_fig4)


