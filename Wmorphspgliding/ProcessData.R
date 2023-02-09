## ----------- Main Processing Data Script ------
## ------------ Date: 08-Sep-17 ------------------
library(R.matlab)
working_directory = "ENTERYOURDIRECTORYHERE"
# ----------- Set the working directory 
setwd(working_directory)

# --------------- Read in all Data Files ---------------
data          <- read.csv('2017_11_24_gullinfo.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
varying.data  <- read.csv('ConcatenatedRunData.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
static.data   <- read.csv('ConcatenatedStaticData.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
wtangles      <- read.csv('windtunnelangles.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
COGnum        <- read.csv('COGcalcs.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )

setwd(working_directory)
source("tareremoval.R")
source("rotate2WTaxis.R")
source("volts2forces.R")
#------------- set variables ----------------
N                    <- 40960 #number of samples taken during all force sensor and pitot tube readings
len_data             <- length(varying.data$WingID)
VtoPa                <- 2*133.322; #VERIFIED FROM CODE from pressure transducer

#----------------------------CLEAN UP DATA COMING IN ----------------------------
#-------------------- 1: Zero runs that weren't properly zeroed ----------------
setwd(working_directory)
staticdata_17_0264  <- readMat("2017_05_03_125deg_17_0264_Test3_Large_10msstatic1.mat")
for (l in 1:6){varying.data[(varying.data$Test.No == 3 & varying.data$WingID == "17_0264"),(l+14)] = varying.data[(varying.data$Test.No == 3 & varying.data$WingID == "17_0264"),(l+14)] - staticdata_17_0264$data.mean[1,l+1]}
remove(staticdata_17_0264) #cleans up variable list

#---------------------2: Match up Wing Specific information from the GullInfo file -------
varying.data$Angle.True    <- wtangles$elbow.angle[match(varying.data$WingID, wtangles$WingID)] #true measured elbow angle
varying.data$Wing.Weight   <- (as.numeric(data$Wing.Weight[match(varying.data$WingID, data$WingID)]) - 4.1 - 2.8)*0.001; #removes the weight of the nuts and test rod and changes from g to kg
varying.data$COG.span      <- as.numeric(COGnum$COG.span[match(varying.data$WingID, COGnum$WingID)])*0.001; #mm to m
varying.data$COG.chord     <- as.numeric(COGnum$COG.chord[match(varying.data$WingID, COGnum$WingID)])*0.001; #mm to m
varying.data$COG.z         <- as.numeric(COGnum$COG.z[match(varying.data$WingID, COGnum$WingID)])*0.001; #mm to m
varying.data$COG.span.std  <- as.numeric(COGnum$COG.span.std[match(varying.data$WingID, COGnum$WingID)])*0.001; #mm to m
varying.data$COG.chord.std <- as.numeric(COGnum$COG.chord.std[match(varying.data$WingID, COGnum$WingID)])*0.001; #mm to m
varying.data$COG.z.std     <- as.numeric(COGnum$COG.z.std[match(varying.data$WingID, COGnum$WingID)])*0.001; #mm to m
varying.data$Wing.Area     <- as.numeric(data$Wing.Area.True[match(varying.data$WingID, data$WingID)])*0.0001 #cm^2 to m^2
varying.data$Wing.Area.sd  <- as.numeric(data$Wing.Area.sd[match(varying.data$WingID, data$WingID)])*0.0001 #cm^2 to m^2
varying.data$Wing.Eff.Span <- as.numeric(data$Eff.Span[match(varying.data$WingID, data$WingID)])*0.01 #cm to m
varying.data$Root.Chord    <- as.numeric(data$Root.Chord.True[match(varying.data$WingID, data$WingID)])*0.01  #cm to m
varying.data$Root.Chord.sd <- as.numeric(data$Root.Chord.sd[match(varying.data$WingID, data$WingID)])*0.01  #cm to m
varying.data               <- merge(varying.data, data[,c("WingID","Grid","I","H","A","F", "Bin", "Support.Bin", "Angle.Set")],
                                    by = c("WingID", "Grid"), all.x = TRUE, all.y = FALSE)
varying.data$TI                             = 0
varying.data$TI[which(varying.data$Grid =="NoGrid")] = 0.04
varying.data$TI[which(varying.data$Grid =="Rd38")]   = 1.42
varying.data$TI[which(varying.data$Grid =="Large")]  = 4.61

varying.data$up_or_down = "up"
varying.data$check_direction = 0
for (i in 1:(length(varying.data$WingID)-1)){
  varying.data$check_direction[i] = varying.data$curr.angle[i+1] - varying.data$curr.angle[i]; #find the difference between
  if (varying.data$check_direction[i] < 0){
    varying.data$up_or_down[i] = "down"}}

#-------------------- 2: Update variable data types ----------------
varying.data$curr.angle           = as.integer(varying.data$curr.angle); 
for (i in 6:8){varying.data[,i]   <- as.numeric(varying.data[,i]);}
for (i in 10:23){varying.data[,i] <- as.numeric(varying.data[,i]);}
varying.data$I = as.numeric(varying.data$I)
varying.data$H = as.numeric(varying.data$H)
varying.data$A = as.numeric(varying.data$A)
varying.data$F = as.numeric(varying.data$F)

#------------ Precision uncertainty calculations --------------
varying.data$sig_Vfx = varying.data$Vfx_std/(sqrt(N/(varying.data$Vfx_indsamp)))
varying.data$sig_Vfy = varying.data$Vfy_std/(sqrt(N/(varying.data$Vfy_indsamp)))
varying.data$sig_Vfz = varying.data$Vfz_std/(sqrt(N/(varying.data$Vfz_indsamp)))
varying.data$sig_Vmx = varying.data$Vmx_std/(sqrt(N/(varying.data$Vmx_indsamp)))
varying.data$sig_Vmy = varying.data$Vmy_std/(sqrt(N/(varying.data$Vmy_indsamp)))
varying.data$sig_Vmz = varying.data$Vmz_std/(sqrt(N/(varying.data$Vmz_indsamp)))

# ----------- Calculate the value of volts create at each every angle of attack just due to the support system ----------
support_results      <- tareremoval(varying.data,N,len_data)
revised.support.data <- support_results[[1]] 
sig                  <- support_results[[2]]
### ADJUSTED VOLTAGES - SUPPORT DATA REMOVED
varying.data$Vfx_adj = varying.data$Vfx_mean-revised.support.data[,1]
varying.data$Vfy_adj = varying.data$Vfy_mean-revised.support.data[,2]
varying.data$Vfz_adj = varying.data$Vfz_mean-revised.support.data[,3]

varying.data$Vmx_adj = varying.data$Vmx_mean-revised.support.data[,4]
varying.data$Vmy_adj = varying.data$Vmy_mean-revised.support.data[,5]
varying.data$Vmz_adj = varying.data$Vmz_mean-revised.support.data[,6]

### ABSOLUTE UNCERTAINTY OF THE ADJUSTED VOLTAGES
varying.data$Vfx_adj_std = sqrt((varying.data$sig_Vfx^2)+(sig[,1]^2))
varying.data$Vfy_adj_std = sqrt((varying.data$sig_Vfy^2)+(sig[,2]^2))
varying.data$Vfz_adj_std = sqrt((varying.data$sig_Vfz^2)+(sig[,3]^2))

varying.data$Vmx_adj_std = sqrt((varying.data$sig_Vmx^2)+(sig[,4]^2))
varying.data$Vmy_adj_std = sqrt((varying.data$sig_Vmy^2)+(sig[,5]^2))
varying.data$Vmz_adj_std = sqrt((varying.data$sig_Vmz^2)+(sig[,6]^2))

# ------ Change Volts to Forces (N) and Moments (Nm) -----------------------------------
varying.data <- volts2forces(varying.data, len_data)

# ------ Transofrm the coordinate system to aerodynamic variables
varying.data <- rotate2WTaxis(varying.data, len_data)

# ---------- Airspeed and Dynamic Pressure Calculation -----------
# bring in the p offset for each run by wing id and grid run
varying.data            <- merge(varying.data, data[,c("WingID","Grid","P.offset")], by = c("WingID", "Grid"), all.x = TRUE, all.y = FALSE)
varying.data$U          <- sqrt(((varying.data$pt_mean-varying.data$P.offset)*VtoPa*2)/varying.data$rho.air) 
varying.data$dynrho     <- (varying.data$pt_mean-varying.data$P.offset)*VtoPa
# uncertainty
varying.data$pt_uncert  <- varying.data$pt_std/(sqrt(N/(varying.data$pt_indsamp)))
static.data$P.offset_sd <- static.data$pt_std/sqrt(N/(static.data$pt_indsamp))
varying.data            <- merge(varying.data, static.data[,c("WingID","Grid","P.offset_sd")],by = c("WingID", "Grid"), all.x = TRUE, all.y = FALSE)
varying.data$dynrho_sd  <- (sqrt(varying.data$pt_uncert^2+varying.data$P.offset_sd^2 + 2*((0.003*5)^2)))*VtoPa # includes both bias and random error

# ------ NON-DIMENSIONALIZATION
varying.data$CLift    = varying.data$Lift/(varying.data$dynrho*varying.data$Wing.Area)
varying.data$CD       = varying.data$Drag/(varying.data$dynrho*varying.data$Wing.Area)
varying.data$Cm       = varying.data$Pitch/(varying.data$dynrho*varying.data$Wing.Area*varying.data$Root.Chord)
varying.data$L_D      = varying.data$Lift/varying.data$Drag
varying.data$AR       = 2*(varying.data$Wing.Eff.Span)^2/(varying.data$Wing.Area)
varying.data$Re       = varying.data$rho.air*varying.data$U*varying.data$Root.Chord/varying.data$visc.air

#------ absolute uncertainty
varying.data$CLift_sd = abs(varying.data$CLift)*sqrt((varying.data$Lift_sd/varying.data$Lift)^2  + (varying.data$dynrho_sd/varying.data$dynrho)^2 + (varying.data$Wing.Area.sd/varying.data$Wing.Area)^2)
varying.data$CD_sd    = abs(varying.data$CD)*sqrt((varying.data$Drag_sd/varying.data$Drag)^2     + (varying.data$dynrho_sd/varying.data$dynrho)^2 + (varying.data$Wing.Area.sd/varying.data$Wing.Area)^2)
varying.data$Cm_sd    = abs(varying.data$Cm)*sqrt((varying.data$Pitch_sd/varying.data$Pitch)^2   + (varying.data$dynrho_sd/varying.data$dynrho)^2 + (varying.data$Wing.Area.sd/varying.data$Wing.Area)^2+ (varying.data$Root.Chord.sd/varying.data$Root.Chord)^2)
varying.data$L_D_sd   = abs(varying.data$L_D)*sqrt((varying.data$Lift_sd/varying.data$Lift)^2    + (varying.data$Drag_sd/varying.data$Drag)^2)

## varying.data = ProcessedData.csv