LC_cal_fit <- function(){
  # ------------ Read in calibration files from Matlab -----------
  setwd(working_directory)
  LC_cal_Fx_My        <- read.csv('LC_Cal_Fx_My.csv', header = FALSE, stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
  LC_cal_Fy_Mx        <- read.csv('LC_Cal_Fy_Mx.csv', header = FALSE, stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
  LC_cal_Fz           <- read.csv('LC_Cal_Fz.csv', header = FALSE, stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
  LC_cal_Mz           <- read.csv('LC_Cal_Mz.csv', header = FALSE, stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
  names(LC_cal_Fx_My) <-c("Vfx","Vfy","Vfz","Vmx","Vmy","Vmz","Force","Moment")
  names(LC_cal_Fy_Mx) <-c("Vfx","Vfy","Vfz","Vmx","Vmy","Vmz","Force","Moment")
  names(LC_cal_Fz)    <-c("Vfx","Vfy","Vfz","Vmx","Vmy","Vmz","Force","Moment")
  names(LC_cal_Mz)    <-c("Vfx","Vfy","Vfz","Vmx","Vmy","Vmz","Force","Moment")
  
  LC_forces_cal      <- data.frame(matrix(nrow = 3, ncol = 3))
  LC_moments_cal     <- data.frame(matrix(nrow = 3, ncol = 3))
  LC_forces_std      <- data.frame(matrix(nrow = 3, ncol = 3))
  LC_moments_std     <- data.frame(matrix(nrow = 3, ncol = 3))
  # ------ Want to find the best fit line for the calibration -------
  
  # ------ Loading due to Force on X axis -----------
  vfx_byFx_fit <- lm(LC_cal_Fx_My$Vfx ~ LC_cal_Fx_My$Force)
  vfx_byFx     <- summary(vfx_byFx_fit)
  vfy_byFx_fit <- lm(LC_cal_Fx_My$Vfy ~ LC_cal_Fx_My$Force)
  vfy_byFx     <- summary(vfy_byFx_fit)
  vfz_byFx_fit <- lm(LC_cal_Fx_My$Vfz ~ LC_cal_Fx_My$Force)
  vfz_byFx     <- summary(vfz_byFx_fit)
  # ------ Loading due to Force on Y axis -----------
  vfy_byFy_fit <- lm(LC_cal_Fy_Mx$Vfy ~ LC_cal_Fy_Mx$Force)
  vfy_byFy     <- summary(vfy_byFy_fit)
  vfx_byFy_fit <- lm(LC_cal_Fy_Mx$Vfx ~ LC_cal_Fy_Mx$Force)
  vfx_byFy     <- summary(vfx_byFy_fit)
  vfz_byFy_fit <- lm(LC_cal_Fy_Mx$Vfz ~ LC_cal_Fy_Mx$Force)
  vfz_byFy     <- summary(vfz_byFy_fit) 
  # ------ Loading due to Force on Z axis -----------
  # note: not calculating effect on x or y due to possible difference between calibration and test condition
  vfz_byFz_fit <- lm(LC_cal_Fz$Vfz ~ LC_cal_Fz$Force)
  vfz_byFz     <- summary(vfz_byFz_fit)
  # ------ Loading due to Moment about X axis -----------
  vmx_byMx_fit <- lm(LC_cal_Fy_Mx$Vmx ~ LC_cal_Fy_Mx$Moment)
  vmx_byMx     <- summary(vmx_byMx_fit)
  vmy_byMx_fit <- lm(LC_cal_Fy_Mx$Vmy ~ LC_cal_Fy_Mx$Moment)
  vmy_byMx     <- summary(vmy_byMx_fit)
  vmz_byMx_fit <- lm(LC_cal_Fy_Mx$Vmz ~ LC_cal_Fy_Mx$Moment) 
  vmz_byMx     <- summary(vmz_byMx_fit)#weight was offset from main axis therefore should not include this effects
  # ------ Loading due to Moment about Y axis -----------
  vmy_byMy_fit <- lm(LC_cal_Fx_My$Vmy ~ LC_cal_Fx_My$Moment)
  vmy_byMy     <- summary(vmy_byMy_fit)
  vmx_byMy_fit <- lm(LC_cal_Fx_My$Vmx ~ LC_cal_Fx_My$Moment) 
  vmx_byMy     <- summary(vmx_byMy_fit)
  vmz_byMy_fit <- lm(LC_cal_Fx_My$Vmz ~ LC_cal_Fx_My$Moment)
  vmz_byMy     <- summary(vmz_byMy_fit) #weight was offset from main axis therefore should not include this effects
  
  # ------ Loading due to  Moment about Z axis -----------
  # note: not calculating effect on x or y due to this measurement was done at an offset that by neccessity affects moments about x and y
  vmz_byMz_fit <- lm(LC_cal_Mz$Vmz ~ LC_cal_Mz$Moment)
  vmz_byMz     <- summary(vmz_byMz_fit)
  ## ------------- Create the Force based Matrix and Std Error Values ---------
  # Note: from the note above, we did not account for interavtion betwen x-z and y-z due to the orientation when loading. This was assumed to be negligable.
  LC_forces_cal[1,1]     =  vfx_byFx$coefficients[2]
  LC_forces_cal[2,1]     =  vfy_byFx$coefficients[2]
  LC_forces_cal[3,1]     =  0 #vfz_byFx$coefficients[2]
  
  LC_forces_cal[1,2]     =  vfx_byFy$coefficients[2]
  LC_forces_cal[2,2]     =  vfy_byFy$coefficients[2]
  LC_forces_cal[3,2]     =  0 #vfz_byFy$coefficients[2]
  
  LC_forces_cal[1,3]     =  0
  LC_forces_cal[2,3]     =  0
  LC_forces_cal[3,3]     =  vfz_byFz$coefficients[2]
  
  LC_forces_std[1,1]     =  vfx_byFx$coefficients[4]
  LC_forces_std[2,1]     =  vfy_byFx$coefficients[4]
  LC_forces_std[3,1]     =  vfz_byFx$coefficients[4]
  
  LC_forces_std[1,2]     =  vfx_byFy$coefficients[4]
  LC_forces_std[2,2]     =  vfy_byFy$coefficients[4]
  LC_forces_std[3,2]     =  0 #vfz_byFy$coefficients[4]
  
  LC_forces_std[1,3]     =  0
  LC_forces_std[2,3]     =  0
  LC_forces_std[3,3]     =  vfz_byFz$coefficients[4]
  
  ## ------------- Create the Moment based Matrix and Std Error Values ---------
  LC_moments_cal[1,1]     =  vmx_byMx$coefficients[2]
  LC_moments_cal[2,1]     =  vmy_byMx$coefficients[2]
  LC_moments_cal[3,1]     =  0 #vmz_byMx$coefficients[2]
  
  LC_moments_cal[1,2]     =  vmx_byMy$coefficients[2]
  LC_moments_cal[2,2]     =  vmy_byMy$coefficients[2]
  LC_moments_cal[3,2]     =  0 #vmz_byMy$coefficients[2]
  
  LC_moments_cal[1,3]     =  0
  LC_moments_cal[2,3]     =  0
  LC_moments_cal[3,3]     =  vmz_byMz$coefficients[2]
  
  LC_moments_std[1,1]     =  vmx_byMx$coefficients[4]
  LC_moments_std[2,1]     =  vmy_byMx$coefficients[4]
  LC_moments_std[3,1]     =  0 #vmz_byMx$coefficients[4]
  
  LC_moments_std[1,2]     =  vmx_byMy$coefficients[4]
  LC_moments_std[2,2]     =  vmy_byMy$coefficients[4]
  LC_moments_std[3,2]     =  0 #vmz_byMy$coefficients[4]
  
  LC_moments_std[1,3]     =  0
  LC_moments_std[2,3]     =  0
  LC_moments_std[3,3]     =  vmz_byMz$coefficients[4]
  
  return(list(LC_forces_cal,LC_moments_cal,LC_forces_std,LC_moments_std))
  
}
