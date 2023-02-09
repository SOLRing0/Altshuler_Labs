## ---------------------------------------------------------------------------------------
## --------------------------------- DEFINE CONSTANTS FOR ALL PLOTS   ------------------------
## ---------------------------------------------------------------------------------------
library(plyr)
library(ggplot2)


grids              <- c("NoGrid","Rd38","Large")
Wing_Number        <- c("Gull4", "17_0353", "17_0340", "17_0294", "Gull1", "17_0373_2", "17_0069", "Gull2", "17_0264", "16_1080","17_0243","Gull3", "17_0365", "17_0285")


#---------------------------- Colours ----------------------------
#colour for PCA plots
col_lab     = "#CCCCCC"
col_invivo  = "#E9E9E9"
line_invivo = "#3C4049"
col_wt      = "#7E7F80"
line_wt     = "#1D1D21"

#---- Blue Scale
cc1           <- scales::seq_gradient_pal("#42B3D5", "#1A237E", "Lab")(seq(0,1,length.out=6))
cc2           <- scales::seq_gradient_pal("#DCEDC8", "#42B3D5", "Lab")(seq(0,1,length.out=7))
cc_TI         <- c(cc2,cc1)
TIcolors      <- c(cc_TI[13],cc_TI[4],cc_TI[10])
TIlinecolors <-c("#060C4A","#418B88","#0E4A83")

#---- Purple scale Graphiq
cc4  <- scales::seq_gradient_pal("#6A1B9A", "black", "Lab")(seq(0,1,length.out=16))
cc3  <- scales::seq_gradient_pal("#E85285", "#6A1B9A", "Lab")(seq(0,1,length.out=42))
cc2  <- scales::seq_gradient_pal("#FFECB3", "#E85285", "Lab")(seq(0,1,length.out=60))
cc1  <- scales::seq_gradient_pal("white", "#FFECB3", "Lab")(seq(0,1,length.out=29))
cc_full   <- c(cc1[26:29],cc2[2:60],cc3[2:40],cc4[2:15])
#Done based on the differences between the angles
cc <- c(cc_full[1],cc_full[16],cc_full[19],cc_full[42],cc_full[47],cc_full[49],cc_full[68],cc_full[80],cc_full[82],cc_full[105],cc_full[106],cc_full[116])
title_ti = c("0.14%","1.47%","4.94%")


cc_wind <- c("70" = cc_full[36],
             "80" = cc_full[46],
             "90" = cc_full[56],
             "100" = cc_full[66],
             "110" = cc_full[76],
             "120" = cc_full[86],
             "130" = cc_full[96],
             "140" = cc_full[106],
             "150" = cc_full[116],
             "160" = "black",
             "170" = "black")

#---------------------------- Themes ----------------------------

th_noleg = theme_classic() +  
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=12), 
        legend.position = "none") #removes legend so I can put my own

th = theme_classic() +  
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=12))

th_blank = theme(axis.line=element_blank(),axis.text.x=element_blank(),
                 axis.text.y=element_blank(),axis.ticks=element_blank(),
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank(),legend.position="none",
                 panel.background=element_blank(),panel.border=element_blank(),
                 panel.grid.major=element_blank(),
                 panel.grid.minor=element_blank())

#---------------------------- Limits and Breaks ----------------------------

plot_list = list()

elbow_x      = c(30,150)
elbow_limits = c(30,170)
manus_limits = c(100,185)
pc1_limits   = c(-0.30,0.20)
pc2_limits   = c(-0.18,0.18)

elbow_breaks = c(30,50,70,90,110,130,150,170)
manus_breaks = c(100,120,140,160,180)
pc1_breaks   = c(-0.30,-0.20,-0.10,0,0.10,0.20)
pc2_breaks   = c(-0.18,0,0.18)


#------- All LABELS FOR THE FUNCTIONS -------
lab_lift       = "Coefficient of Lift"
lab_lift_multi = c("Coefficient of Lift","","")
lab_lift_multi_sp = c(paste("Specific speed \n coefficient of lift"),"","")
lab_drag       = "Coefficient of Drag"
lab_drag_multi = c("","","Coefficient of Drag")
lab_drag_sp       = paste("Specific speed \n coefficient of drag")
lab_LD         = expression(paste("C" [L], "/C" [D] ))
lab_LD_multi   = c(lab_LD ,"","")
lab_LD_max     = expression(paste("(C" [L], "/C" [D],")"['max'] ))
lab_side       = c("Coefficient of Side Force","","")
lab_ti         = paste("Turbulence Intensity (%)")
lab_elbow      = "Elbow Angle (°)"
lab_manus      = "Manus Angle (°)"
lab_pitch      = "Coefficient of Pitching Moment"
lab_pitch_multi = c("Coefficient ofPitching Moment","","")
lab_cm         = expression(paste(delta,"C" [M],"/", delta, "C" [L] ))
lab_cm0        = expression(paste("C" [M0]))
lab_alpha      = expression(paste(alpha,"(°)"))
lab_yaw        = "Coefficient of Yaw"
lab_yaw_multi  = c(lab_yaw ,"","")
lab_roll       = "Coefficient of Roll"
lab_camber     = expression(paste('Camber Factor (',beta,")"))
lab_CL_max     = expression(paste("C" [L['max']]))
lab_CD_min     = expression(paste("C" [D['min']]))
lab_dydr       = expression(paste("(C" [n],"/C" [l],")"['min'] ))
lab_cnmin      = expression(paste("C" [n['min']]))

## ---------------------------------------------------------------------------------------
## --------------------------------- DEFNIE ALL NECESSARY FUNCTIONS   ------------------------
## ---------------------------------------------------------------------------------------

#create function to retrieve the legend as an object
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
