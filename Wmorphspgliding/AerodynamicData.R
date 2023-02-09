## ---------------------------------------------------------------------------------------
## --------------------------------- INITIALIZE TO PREP FOR PLOTTING SCRIPTS   ------------------------
## ---------------------------------------------------------------------------------------

## --------- Updated: 08-Sep-17
library(ggplot2)
library(gridExtra)
library(grid) #needed to change title size in gridarrange
library(data.table)
library(geomorph)

setwd(working_directory)
dat      <- read.csv('ProcessedData.csv',    stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c(""))
gullinfo <- read.csv('2017_11_24_gullinfo.csv',   stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
wtangles <- read.csv('windtunnelangles.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )

## ---------------------------------------------------------------------------------------
## --------------------------------- UPDATE data   ------------------------
## ---------------------------------------------------------------------------------------
#------ Remove bad runs
dat <- dat[-c(which(dat$WingID =="17_0243")),]                        #broke mid test
dat <- dat[-c(which(dat$WingID =="17_0353")),]                        #appears to have a strange reading at high angles of attack
dat <- dat[-c(which(dat$WingID =="16_1080" & dat$Grid == "Large")),]  #broke mid test
dat <- dat[-c(which(dat$WingID =="Gull2"   & dat$Grid == "Large")),]  #non-constant pressure measurements during this run

##----- Effective Span saved in ProcessedData.csv is incorrect - update here and fix AR to solve
dat$Wing.Eff.Span  = gullinfo$Eff.Span[match(dat$WingID, gullinfo$WingID)]*0.01 #cm to m
dat$Wing.Weight    = gullinfo$Wing.Weight[match(dat$WingID, gullinfo$WingID)]
dat$AR             = 2*(dat$Wing.Eff.Span)^2/(dat$Wing.Area)

for (i in 55:127){dat[,i] <- as.numeric(dat[,i]);}

## ----- Fix Turbulence intensity values - hard coded these values are from the turbulence matlab file
dat$TI[which(dat$Grid =="NoGrid")] = 0.04
dat$TI[which(dat$Grid =="Rd38")]   = 1.42
dat$TI[which(dat$Grid =="Large")]  = 4.61

## ---------------------------------------------------------------------------------------
## --------------------------------- create stall subset   ------------------------
## ---------------------------------------------------------------------------------------
#this subset inclues only the linear portion of all curves (pre-stall) - to be used for slope fitting
# Updated:  05-Sep-17 
stall_subset = dat[c(which(dat$WingID == "Gull1" & dat$Grid == "NoGrid" & dat$curr.angle <= -4),
                     which(dat$WingID == "Gull1" & dat$Grid == "Rd38"   & dat$curr.angle <= -4),
                     which(dat$WingID == "Gull1" & dat$Grid == "Large"  & dat$curr.angle <= -4),
                     
                     which(dat$WingID == "Gull2" &  dat$Grid == "NoGrid" & (dat$curr.angle <= 0 & dat$curr.angle >= -20)),
                     which(dat$WingID == "Gull2" &  dat$Grid == "Rd38"   & (dat$curr.angle <= 0 & dat$curr.angle >= -20)),
                     which(dat$WingID == "Gull2" &  dat$Grid == "Large"  & (dat$curr.angle <= 0 & dat$curr.angle >= -20)),
                     
                     which(dat$WingID == "Gull3" & dat$Grid == "NoGrid" & (dat$curr.angle <= 37 & dat$curr.angle >= 4)),
                     which(dat$WingID == "Gull3" & dat$Grid == "Rd38"   & (dat$curr.angle <= 37 & dat$curr.angle >= 4)),
                     which(dat$WingID == "Gull3" & dat$Grid == "Large"  & (dat$curr.angle <= 37 & dat$curr.angle >= 4)),
                     
                     which(dat$WingID == "Gull4" & dat$Grid == "NoGrid" & (dat$curr.angle <= 28 & dat$curr.angle >= -5)),
                     which(dat$WingID == "Gull4" & dat$Grid == "Rd38"   & (dat$curr.angle <= 28 & dat$curr.angle >= -5)),
                     which(dat$WingID == "Gull4" & dat$Grid == "Large"  & (dat$curr.angle <= 28 & dat$curr.angle >= -5)),            
                     
                     #which(dat$WingID == "17_0243" & dat$Grid == "NoGrid" & (dat$curr.angle <= 34 & dat$curr.angle >= -15)), 
                     
                     which(dat$WingID == "17_0264" & dat$Grid == "NoGrid" & (dat$curr.angle <= 4 & dat$curr.angle >= -20)),
                     which(dat$WingID == "17_0264" & dat$Grid == "Rd38"   & (dat$curr.angle <= 4 & dat$curr.angle >= -20)),
                     which(dat$WingID == "17_0264" & dat$Grid == "Large"  & (dat$curr.angle <= 4 & dat$curr.angle >= -20)),
                     
                     which(dat$WingID == "17_0340" & dat$Grid == "NoGrid" & (dat$curr.angle <= 10 & dat$curr.angle >= -20)),
                     which(dat$WingID == "17_0340" & dat$Grid == "Rd38"   & (dat$curr.angle <= 10 & dat$curr.angle >= -20)),
                     which(dat$WingID == "17_0340" & dat$Grid == "Large"  & (dat$curr.angle <= 10 & dat$curr.angle >= -25)),
                     
                     which(dat$WingID == "17_0353" & dat$Grid == "NoGrid" & (dat$curr.angle <= 20 & dat$curr.angle >= -5)),
                     which(dat$WingID == "17_0353" & dat$Grid == "Rd38"   & (dat$curr.angle <= 20 & dat$curr.angle >= -5)),
                     
                     which(dat$WingID == "16_1080" & dat$Grid == "NoGrid" & (dat$curr.angle <= 20 & dat$curr.angle >= 0)),
                     which(dat$WingID == "16_1080" & dat$Grid == "Rd38"   & (dat$curr.angle <= 20 & dat$curr.angle >= 0)),
                     #which(dat$WingID == "16_1080" & dat$Grid == "Large"  & (dat$curr.angle <= 28 & dat$curr.angle >= 0)),
                     
                     which(dat$WingID == "17_0365" & dat$Grid == "NoGrid" & (dat$curr.angle <= 37 & dat$curr.angle >= 6)),
                     which(dat$WingID == "17_0365" & dat$Grid == "Rd38"   & (dat$curr.angle <= 37 & dat$curr.angle >= 6)),
                     which(dat$WingID == "17_0365" & dat$Grid == "Large"  & (dat$curr.angle <= 37 & dat$curr.angle >= 6)),
                     
                     which(dat$WingID == "17_0373_2" & dat$Grid == "NoGrid" & (dat$curr.angle <= 2 & dat$curr.angle >= -30)),
                     which(dat$WingID == "17_0373_2" & dat$Grid == "Rd38"   & (dat$curr.angle <= 2 & dat$curr.angle >= -30)),
                     which(dat$WingID == "17_0373_2" & dat$Grid == "Large"  & (dat$curr.angle <= 2 & dat$curr.angle >= -30)),
                     
                     which(dat$WingID == "17_0294" & dat$Grid == "NoGrid" & (dat$curr.angle <= 10 & dat$curr.angle >= -10)),
                     
                     which(dat$WingID == "17_0069" & dat$Grid == "NoGrid" & (dat$curr.angle <= 4 & dat$curr.angle >= -25)),
                     which(dat$WingID == "17_0069" & dat$Grid == "Rd38"   & (dat$curr.angle <= 4 & dat$curr.angle >= -25)),
                     which(dat$WingID == "17_0069" & dat$Grid == "Large"  & (dat$curr.angle <= 4 & dat$curr.angle >= -25)),
                     
                     which(dat$WingID == "17_0285" & dat$Grid == "NoGrid" & (dat$curr.angle <= 37 & dat$curr.angle >= 18)),
                     which(dat$WingID == "17_0285" & dat$Grid == "Rd38"   & (dat$curr.angle <= 37 & dat$curr.angle >= 18))),]

## ---------------------------------------------------------------------------------------
## --------------------------------- CREATE FIT MODEL DATA FRAME   ------------------------
## ---------------------------------------------------------------------------------------

#create fit models for the linear model
fit_results        <- data.frame(matrix(nrow = 31, ncol = 9))
names(fit_results) <- c("Bin","Grid","TI","slopeup","slopedown","slopemean","interceptup","interceptdown","interceptmean")

count = 1
for (i in 1:3){
  for (j in 1:14){
    new_dat = subset(stall_subset, Grid %in% grids[i] & up_or_down %in% c("up") & WingID %in%  c(Wing_Number[j]))
    new_dat2 = subset(stall_subset, Grid %in% grids[i] & up_or_down %in% c("down") & WingID %in%  c(Wing_Number[j]))
    new_dat3  = subset(dat, Grid %in% grids[i] & WingID %in%  c(Wing_Number[j]))
    new_dat3 = new_dat3[order(as.numeric(abs(new_dat3$Lift))),]
    if (nrow(new_dat)==0) {next}
    
    fit3     <- lm(Cm~CLift, data = new_dat3) #pooled cycles
    fit_sum_pooled <- summary(fit3) #pooled cycles
    
    fit_results$Grid[count]       = grids[i]
    fit_results$Bin[count]        = unique(new_dat3$Bin)
    fit_results$Angle.True[count] = unique(new_dat3$Angle.True)
    fit_results$TI[count]         = unique(new_dat3$TI)
    fit_results$WingID[count]     = unique(new_dat3$WingID)
    
    fit_results$slopemean[count]  <- fit_sum_pooled$coefficients[2]
    fit_results$slopesd[count]    <- fit_sum_pooled$coefficients[4]
    fit_results$slope2.5[count]      <- confint(fit3)[2]
    fit_results$slope97.5[count]     <- confint(fit3)[4]
    
    fit_results$interceptmean[count]  <- fit_sum_pooled$coefficients[1]
    fit_results$interceptsd[count]    <- fit_sum_pooled$coefficients[3]
    fit_results$intercept2.5[count]      <- confint(fit3)[1]
    fit_results$intercept97.5[count]     <- confint(fit3)[3]
    
    fit_results$xinterceptmean[count]  <- -fit_sum_pooled$coefficients[1]/fit_sum_pooled$coefficients[2]
    fit_results$xinterceptsd[count]    <- -fit_sum_pooled$coefficients[3]/fit_sum_pooled$coefficients[4]
    
    ### ------- Calculate the moment about the aerodynamic center at zero-lift -----------
    # loops until it finds the index which represents the smallest lift value on the other side of zero
    for (m in 1:length(new_dat3$WingID)){
      if (new_dat3$Lift[1] > 0 & new_dat3$Lift[1] > new_dat3$Lift[m+1]){
        ind = m+1;
        break}
      if (new_dat3$Lift[1] < 0 & new_dat3$Lift[1] < new_dat3$Lift[m+1]){
        ind = m+1;
        break}}
    
    #Interpolates to find values when Lift = 0
    fit_results$Drag_0[count]  = (-new_dat3$Lift[1]/(new_dat3$Lift[ind]-new_dat3$Lift[1]))*(new_dat3$Drag[ind]-new_dat3$Drag[1]) + new_dat3$Drag[1]
    vertoff                    = (-new_dat3$Lift[1]/(new_dat3$Lift[ind]-new_dat3$Lift[1]))*(new_dat3$vert_offset[ind]-new_dat3$vert_offset[1]) + new_dat3$vert_offset[1]
    fit_results$Pitch_0[count] = (-new_dat3$Lift[1]/(new_dat3$Lift[ind]-new_dat3$Lift[1]))*(new_dat3$Pitch0[ind]-new_dat3$Pitch0[1]) + new_dat3$Pitch0[1]
   
     #Aerodynamic center pitching moment
    fit_results$Cm_AC[count]      = (-new_dat3$Lift[1]/(new_dat3$Lift[ind]-new_dat3$Lift[1]))*(new_dat3$Cm[ind]-new_dat3$Cm[1]) + new_dat3$Cm[1]
    
    #COMPUTE AERODYNAMIC CENTER AT THE MAXIMAL AERODYNAMIC EFFICIENCY
    CL1                       = which(abs(new_dat3$L_D)==max(abs(new_dat3$L_D)))
   
    fit_results$maxCdi[count] = new_dat3$CLift[CL1]^2/(pi*new_dat3$AR[CL1])
    fit_results$Cd[count]     = new_dat3$CD[CL1]
    fit_results$Cdtrue[count] = new_dat3$CD[CL1] - fit_results$maxCdi[count]
    fit_results$L_Dtrue[count]= new_dat3$CLift[CL1]/fit_results$Cdtrue[count]
    
    #------------------ ERROR PROP0GATION ------------------------
    fit_results$Drag_0_sd[count]   = sqrt((abs((-new_dat3$Lift[1]/(new_dat3$Lift[ind]-new_dat3$Lift[1]))*(new_dat3$Drag[ind]-new_dat3$Drag[1]))*
                                             sqrt((new_dat3$Lift_sd[1]/new_dat3$Lift[1])^2+
                                                    (sqrt(new_dat3$Lift_sd[ind]^2 + new_dat3$Lift_sd[1]^2)/(new_dat3$Lift[ind]-new_dat3$Lift[1]))^2+
                                                    (sqrt(new_dat3$Drag_sd[ind]^2 + new_dat3$Drag_sd[1]^2)/(new_dat3$Drag[ind]-new_dat3$Drag[1]))^2))^2+
                                            (new_dat3$Drag_sd[1])^2)
    #need to add allowance for if vertical offset doesn't change
    if ((new_dat3$vert_offset[ind]-new_dat3$vert_offset[1])==0){
      vertoff_sd = new_dat3$vert_offset_sd[1]
    }else{
    vertoff_sd                     = sqrt((abs((-new_dat3$Lift[1]/(new_dat3$Lift[ind]-new_dat3$Lift[1]))*(new_dat3$vert_offset[ind]-new_dat3$vert_offset[1]))*
                                             sqrt((new_dat3$Lift_sd[1]/new_dat3$Lift[1])^2+
                                                    (sqrt(new_dat3$Lift_sd[ind]^2 + new_dat3$Lift_sd[1]^2)/(new_dat3$Lift[ind]-new_dat3$Lift[1]))^2+
                                                    (sqrt(new_dat3$vert_offset_sd[ind]^2 + new_dat3$vert_offset_sd[1]^2)/(new_dat3$vert_offset[ind]-new_dat3$vert_offset[1]))^2))^2+
                                            (new_dat3$vert_offset_sd[1])^2)
    }
    fit_results$Pitch_0_sd[count]  = sqrt((abs((-new_dat3$Lift[1]/(new_dat3$Lift[ind]-new_dat3$Lift[1]))*(new_dat3$Pitch0[ind]-new_dat3$Pitch0[1]))*
                                             sqrt((new_dat3$Lift_sd[1]/new_dat3$Lift[1])^2+
                                                    (sqrt(new_dat3$Lift_sd[ind]^2 + new_dat3$Lift_sd[1]^2)/(new_dat3$Lift[ind]-new_dat3$Lift[1]))^2+
                                                    (sqrt(new_dat3$Pitch0_sd[ind]^2 + new_dat3$Pitch0_sd[1]^2)/(new_dat3$Pitch0[ind]-new_dat3$Pitch0[1]))^2))^2+
                                            (new_dat3$Pitch0_sd[1])^2)
    
    fit_results$Cm_AC_sd[count]   = sqrt((abs((-new_dat3$Lift[1]/(new_dat3$Lift[ind]-new_dat3$Lift[1]))*(new_dat3$Cm[ind]-new_dat3$Cm[1]))*
                                            sqrt((new_dat3$Lift_sd[1]/new_dat3$Lift[1])^2+
                                                   (sqrt(new_dat3$Lift_sd[ind]^2+new_dat3$Lift_sd[1]^2)/(new_dat3$Lift[ind]-new_dat3$Lift[1]))^2+
                                                   (sqrt(new_dat3$Cm_sd[ind]^2+new_dat3$Cm_sd[1]^2)/(new_dat3$Cm[ind]-new_dat3$Cm[1]))^2))^2+
                                           (new_dat3$Cm_sd[1])^2)
    
    count = count + 1
  }
}

dat <- merge(dat, fit_results[,c("WingID","Grid","xac","yac","zac","xac_sd","yac_sd","zac_sd")], by = c("WingID", "Grid"), all.x = TRUE, all.y = FALSE)

#-------- calculate the error due to additional induced drag over root
dat$Cdi            = dat$CLift^2/(pi*dat$AR)
dat$Cdtrue         = dat$CD - dat$Cdi 
dat$L_Dtrue        = dat$CLift/dat$Cdtrue

## ---------------------------------------------------------------------------------------
## --------------------------------- CREATE MAXIMUM RESULTS DATA FRAME   ------------------------
## ---------------------------------------------------------------------------------------

max_results <- data.frame(matrix(nrow = 31, ncol = 14))
names(max_results) <- c("Bin","Grid","TI","MaxL_D", "MaxL_Dup","MaxL_Ddown","MaxL_Dsd","CLmax", "CLmaxup","CLmaxdown","CDmin","CDminup","CDmindown","GlideAngle")

count = 1
for (i in 1:3){
  for (j in 1:14){
    new_dat  = subset(dat, Grid %in% grids[i] & up_or_down %in% c("up") & WingID %in%  c(Wing_Number[j]))
    new_dat2 = subset(dat, Grid %in% grids[i] & up_or_down %in% c("down") & WingID %in%  c(Wing_Number[j]))
    if (nrow(new_dat)==0) { next}#stop if this configuration doesnt exist
    
    max_results$Grid[count]       = grids[i]
    max_results$Bin[count]        = unique(new_dat$Bin)
    max_results$Angle.True[count] = unique(new_dat$Angle.True)
    max_results$TI[count]         = unique(new_dat$TI)
    max_results$WingID[count]     = unique(new_dat$WingID)
    
    #----- Lift to Drag Ratio ------ - COMBINES RESULTS OF UP AND DOWN
    maxup                         = max(new_dat$L_D)
    maxdown                       = max(new_dat2$L_D)
    max_results$MaxL_D[count]     = max(c(maxup,maxdown))
    max_results$MaxL_Dup[count]   = maxup
    max_results$MaxL_Ddown[count] = maxdown
    upsd                          = new_dat$L_D_sd[which.max(new_dat$L_D)]
    downsd                        = new_dat2$L_D_sd[which.max(new_dat2$L_D)]
    max_results$MaxL_Dsd[count]   = c(upsd,downsd)[which.max(c(maxup,maxdown))]

    max_results$dYdR_LDup[count]    =  new_dat$dYdR[which.max(new_dat$L_D)]
    max_results$dYdR_LDdown[count]  =  new_dat2$dYdR[which.max(new_dat2$L_D)]
    max_results$dYdR_LD[count]      = c(max_results$dYdR_LDup[count],max_results$dYdR_LDdown[count])[which.max(c(maxup,maxdown))]
    
    #----- Lift to Drag Ratio ------ - COMBINES RESULTS OF UP AND DOWN
    maxup                             = max(new_dat$L_Dtrue)
    maxdown                           = max(new_dat2$L_Dtrue)
    max_results$MaxL_Dtrue[count]     = max(c(maxup,maxdown))
    max_results$MaxL_Dtrueup[count]   = maxup
    max_results$MaxL_Dtruedown[count] = maxdown
    upsd                              = new_dat$L_D_sd[which.max(new_dat$L_Dtrue)]
    downsd                            = new_dat2$L_D_sd[which.max(new_dat2$L_Dtrue)]
    max_results$MaxL_Dtruesd[count]   = c(upsd,downsd)[which.max(c(maxup,maxdown))]
    
    # ----- Maximum Lift Coefficient -----
    max_results$CLmaxup[count]   = max(new_dat$CLift)
    max_results$CLmaxdown[count] = max(new_dat2$CLift)
    max_results$CLmax[count]     = max(c(max_results$CLmaxup[count],max_results$CLmaxdown[count]))
    upsd                         = new_dat$CLift_sd[which.max(new_dat$CLift)]
    downsd                       = new_dat2$CLift_sd[which.max(new_dat2$CLift)]
    max_results$CLmaxsd[count]   = c(upsd,downsd)[which.max(c(max_results$CLmaxup[count],max_results$CLmaxdown[count]))]
    
    # ----- Minimum Drag Coefficient -----
    max_results$CDminup[count]   = min(new_dat$CD)
    max_results$CDmindown[count] = min(new_dat2$CD)
    max_results$CDmin[count]     = min(c(max_results$CDminup[count],max_results$CDmindown[count]))
    upsd                         = new_dat$CD_sd[which.min(new_dat$CD)]
    downsd                       = new_dat2$CD_sd[which.min(new_dat2$CD)]
    max_results$CDminsd[count]   = c(upsd,downsd)[which.min(c(max_results$CDminup[count],max_results$CDmindown[count]))]
    
    
    max_results$Dragminup[count]   = min(new_dat$Drag)
    max_results$Dragmindown[count] = min(new_dat2$Drag)
    max_results$Dragmin[count]     = min(c(max_results$Dragminup[count],max_results$Dragmindown[count]))
    upsd                           = new_dat$CD_sd[which.min(new_dat$CD)]
    downsd                         = new_dat2$CD_sd[which.min(new_dat2$CD)]
    max_results$Dragminsd[count]   = c(upsd,downsd)[which.min(c(max_results$CDminup[count],max_results$CDmindown[count]))]
    
    count = count + 1
  }
}

max_results$GlideAngle = (180/pi)*atan(1/max_results$MaxL_D)
