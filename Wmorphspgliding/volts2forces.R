volts2forces <-function(varying.data, len_data){
  
  # ------------ Translate Volts to Forces -----------
  # ------------ Date: 08-Sep-17 ---------------------
  
  # ------------ Read in calibration files -----------
  setwd(working_directory)
  source("LC_cal_fit.R")
  calibrationresults <- LC_cal_fit()
  ForceSlope         <- calibrationresults[[1]]
  MomentSlope        <- calibrationresults[[2]]
  ForceSlope_stderr  <- calibrationresults[[3]]
  MomentSlope_stderr <- calibrationresults[[4]]
  forceinverse       <- list()
  momentinverse      <- list()
  data.frame(matrix(nrow = len_data, ncol = 1))
  
  alphaf = data.frame(matrix(nrow = len_data, ncol = 1))
  betaf  = data.frame(matrix(nrow = len_data, ncol = 1))
  gamf   = data.frame(matrix(nrow = len_data, ncol = 1))
  delf   = data.frame(matrix(nrow = len_data, ncol = 1))
  
  sig_alphaf = data.frame(matrix(nrow = len_data, ncol = 1))
  sig_betaf  = data.frame(matrix(nrow = len_data, ncol = 1))
  sig_gamf   = data.frame(matrix(nrow = len_data, ncol = 1))
  sig_delf   = data.frame(matrix(nrow = len_data, ncol = 1))
  
  alpham = data.frame(matrix(nrow = len_data, ncol = 1))
  betam  = data.frame(matrix(nrow = len_data, ncol = 1))
  gamm   = data.frame(matrix(nrow = len_data, ncol = 1))
  delm   = data.frame(matrix(nrow = len_data, ncol = 1))
  
  sig_alpham = data.frame(matrix(nrow = len_data, ncol = 1))
  sig_betam  = data.frame(matrix(nrow = len_data, ncol = 1))
  sig_gamm   = data.frame(matrix(nrow = len_data, ncol = 1))
  sig_delm   = data.frame(matrix(nrow = len_data, ncol = 1))  
  
  # ------ Mean Voltages -> Forces and moments slopes and intercepts calculated as mV/N or mV/Nm hence then multiplication by 1000

  for (i in 1:len_data){
    # -------- Forces -------
    forces <- solve(ForceSlope,c((varying.data$Vfx_adj[i]*1000),(varying.data$Vfy_adj[i]*1000),(varying.data$Vfz_adj[i]*1000)))
    varying.data$Fx[i] <- forces[1]
    varying.data$Fy[i] <- forces[2] 
    varying.data$Fz[i] <- forces[3] 
    ## -------- Forces Uncertainty Propagation -------
    forceinverse[[i]]  <- solve(ForceSlope)
    
    sig_a = ForceSlope_stderr[1,1]
    sig_b = ForceSlope_stderr[1,2]
    sig_c = ForceSlope_stderr[2,1]
    sig_d = ForceSlope_stderr[2,2]
    
    alphaf[i,1] = forceinverse[[1]][1]
    betaf[i,1]  = forceinverse[[1]][4]
    gamf[i,1]   = forceinverse[[1]][2]
    delf[i,1]   = forceinverse[[1]][5]
    
    sig_alphaf[i,1] = sqrt((alphaf[i,1]^4*sig_a^2)              + (alphaf[i,1]^2*gamf[i,1]^2*sig_b^2) + (alphaf[i,1]^2*betaf[i,1]^2*sig_c^2) + (betaf[i,1]^2*gamf[i,1]^2*sig_d^2))
    sig_betaf[i,1]  = sqrt((alphaf[i,1]^2*betaf[i,1]^2*sig_a^2) + (alphaf[i,1]^2*delf[i,1]^2*sig_b^2) + (betaf[i,1]^4*sig_c^2)               + (betaf[i,1]^2*delf[i,1]^2*sig_d^2))
    sig_gamf[i,1]   = sqrt((alphaf[i,1]^2*gamf[i,1]^2*sig_a^2)  + (gamf[i,1]^4*sig_b^2)               + (alphaf[i,1]^2*delf[i,1]^2*sig_c^2)  + (delf[i,1]^2*gamf[i,1]^2*sig_d^2))
    sig_delf[i,1]   = sqrt((betaf[i,1]^2*gamf[i,1]^2*sig_a^2)   + (gamf[i,1]^2*delf[i,1]^2*sig_b^2)   + (betaf[i,1]^2*delf[i,1]^2*sig_c^2)   + (delf[i,1]^4*sig_d^2))
    
    # -------- MOMEMTS -------
    moments <- solve(MomentSlope,c((varying.data$Vmx_adj[i]*1000),(varying.data$Vmy_adj[i]*1000),(varying.data$Vmz_adj[i]*1000)))
    varying.data$Mx[i] <- moments[1]
    varying.data$My[i] <- moments[2]
    varying.data$Mz[i] <- moments[3] 
    ## -------- Forces Uncertainty Propagation -------
    momentinverse[[i]]  <- solve(MomentSlope)
    
    sig_a = MomentSlope_stderr[1,1]
    sig_b = MomentSlope_stderr[1,2]
    sig_c = MomentSlope_stderr[2,1]
    sig_d = MomentSlope_stderr[2,2]
    
    alpham[i,1] = momentinverse[[1]][1]
    betam[i,1]  = momentinverse[[1]][4]
    gamm[i,1]   = momentinverse[[1]][2]
    delm[i,1]   = momentinverse[[1]][5]
    
    sig_alpham[i,1] = sqrt((alpham[i,1]^4*sig_a^2)            + (alpham[i,1]^2*gamm[i,1]^2*sig_b^2) + (alpham[i,1]^2*betam[i,1]^2*sig_c^2) + (betam[i,1]^2*gamm[i,1]^2*sig_d^2))
    sig_betam[i,1]  = sqrt((alpham[i,1]^2*betam[i,1]^2*sig_a^2) + (alpham[i,1]^2*delm[i,1]^2*sig_b^2) + (betam[i,1]^4*sig_c^2)             + (betam[i,1]^2*delm[i,1]^2*sig_d^2))
    sig_gamm[i,1]   = sqrt((alpham[i,1]^2*gamm[i,1]^2*sig_a^2)  + (gamm[i,1]^4*sig_b^2)             + (alpham[i,1]^2*delm[i,1]^2*sig_c^2)  + (delm[i,1]^2*gamm[i,1]^2*sig_d^2))
    sig_delm[i,1]   = sqrt((betam[i,1]^2*gamm[i,1]^2*sig_a^2)   + (gamm[i,1]^2*delm[i,1]^2*sig_b^2)   + (betam[i,1]^2*delm[i,1]^2*sig_c^2)   + (delm[i,1]^4*sig_d^2))
    
  }

  # ------ compute the absolute uncertainty 

  #includes uncertainty of calibration line fit, random and bias error
  varying.data$Fx_std <- sqrt((varying.data$Vfx_adj*1000*sig_alphaf)^2 + (alphaf*varying.data$Vfx_adj_std*1000)^2 + (varying.data$Vfy_adj*1000*sig_betaf)^2 + (varying.data$Vfy_adj_std*1000*betaf)^2)
  varying.data$Fy_std <- sqrt((varying.data$Vfx_adj*1000*sig_gamf)^2   + (gamf*varying.data$Vfx_adj_std*1000)^2   + (varying.data$Vfy_adj*1000*sig_delf)^2  + (varying.data$Vfy_adj_std*1000*delf)^2)
  varying.data$Fz_std <- abs(varying.data$Fz)*sqrt(((varying.data$Vfz_adj_std)/(varying.data$Vfz_adj))^2 + (ForceSlope_stderr[3,3]/ForceSlope[3,3])^2)
  
  varying.data$Mx_std <- sqrt((varying.data$Vmx_adj*1000*sig_alpham)^2 + (alpham*varying.data$Vmx_adj_std*1000)^2 + (varying.data$Vmy_adj*1000*sig_betam)^2 + (varying.data$Vmy_adj_std*1000*betam)^2)
  varying.data$My_std <- sqrt((varying.data$Vmx_adj*1000*sig_gamm)^2   + (gamm*varying.data$Vmx_adj_std*1000)^2   + (varying.data$Vmy_adj*1000*sig_delm)^2 + (varying.data$Vmy_adj_std*1000*delm)^2)
  varying.data$Mz_std <- abs(varying.data$Mz)*sqrt(((varying.data$Vmz_adj_std)/(varying.data$Vmz_adj))^2 + (MomentSlope_stderr[3,3]/MomentSlope[3,3])^2)
  
  return(varying.data)
}
