#--------------Curvature Code ----------------
##Info: This code reads in digitized camber points, rotates them to the wind tunnel orientaion of the wing 
## and finally computes the spanwise camber of the wings.

library(reshape2)

#-------- Read in info and set working directory
setwd(working_directory)
gullinfo <- read.csv('GullInfo_wRuns.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
wtangles <- read.csv('windtunnelangles.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c(""))

setwd(working_directory)
camber            <- read.csv('camber.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
camber$Angle.True <- gullinfo$Angle.Set[match(camber$ID, gullinfo$WingID)]

#Step 1: Make WRIST the beginning location
for (i in 1:length(unique(camber$ID))){
  camber$wristx[camber$ID == unique(camber$ID)[i]] = camber$X[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 3]
  camber$wristy[camber$ID == unique(camber$ID)[i]] = camber$Y[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 3]
  camber$wristz[camber$ID == unique(camber$ID)[i]] = camber$Z[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 3]
}

camber$x_adj = camber$X - camber$wristx
camber$y_adj = camber$Y - camber$wristy
camber$z_adj = camber$Z - camber$wristz

camber$path.no = 0
camber$path.no[camber$Pt.No == 1] = 1 #humerus
camber$path.no[camber$Pt.No == 8] = 2
camber$path.no[camber$Pt.No == 5] = 3
camber$path.no[camber$Pt.No == 3] = 4 #wrist
camber$path.no[camber$Pt.No == 6] = 5
camber$path.no[camber$Pt.No == 4] = 6 #carpo
camber$path.no[camber$Pt.No == 7] = 7
camber$path.no[camber$Pt.No == 2] = 8 #P10


#- Step 2: Calculate angle between x axis and 3->4 (rotate about z)
for (i in 1:length(unique(camber$ID))){  
  #no subtraction as pt1 = 0 - project on z-y plane to calculate angle
  i2 = camber$x_adj[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 4] 
  j2 = camber$y_adj[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 4] 
  k2 = 0
  
  dotprod = abs((i2*1)+(j2*0)+(k2*0))
  len1    = sqrt(i2^2+j2^2+k2^2)
  len2    = sqrt(1^2+0^2+0^2)
  interior = dotprod/(len1*len2)
  camber$thetax[camber$ID == unique(camber$ID)[i]]   = - acos(interior)
  
  if (j2 < 0){
    camber$thetax[camber$ID == unique(camber$ID)[i]] = acos(interior)
  }
  
}

camber$x_adj1 = cos(camber$thetax)*camber$x_adj-sin(camber$thetax)*camber$y_adj
camber$y_adj1 = sin(camber$thetax)*camber$x_adj+cos(camber$thetax)*camber$y_adj
camber$z_adj1 = camber$z_adj


#Step 3: Calculate angle between z axis and 3->4 (rotate about y)
for (i in 1:length(unique(camber$ID))){  
  #--- Vectors of Plane 1 - create vector from two 3D pts 
  #no subtraction as pt1 = 0
  i1 = camber$x_adj1[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 4]
  j1 = 0
  k1 = camber$z_adj1[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 4]
  
  #--- Dot Product of Orthogonal Vectors
  dotprod = abs((i1*1)+(j1*0)+(k1*0))
  len1    = sqrt(i1^2+j1^2+k1^2)
  len2    = sqrt(1^2+0^2+0^2)
  interior = dotprod/(len1*len2)
  camber$thetaz[camber$ID == unique(camber$ID)[i]]   = acos(interior)
  
  if (k1 < 0){
    camber$thetaz[camber$ID == unique(camber$ID)[i]] = - acos(interior)
  }
  
}
camber$x_adj2 = cos(camber$thetaz)*camber$x_adj1+sin(camber$thetaz)*camber$z_adj1
camber$y_adj2 = camber$y_adj1
camber$z_adj2 = -sin(camber$thetaz)*camber$x_adj1+cos(camber$thetaz)*camber$z_adj1


#- Step 4: Calculate angle between y axis and 3->5 (rotate about x)
for (i in 1:length(unique(camber$ID))){  
  #no subtraction as pt1 = 0 - project on x-y plane to calculate angle
  i2 = 0
  j2 = camber$y_adj2[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 6] 
  k2 = camber$z_adj2[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 6] 
  
  dotprod = abs((i2*0)+(j2*1)+(k2*0))
  len1    = sqrt(i2^2+j2^2+k2^2)
  len2    = sqrt(0^2+1^2+0^2)
  interior = dotprod/(len1*len2)
  
  camber$thetay[camber$ID == unique(camber$ID)[i]]   =  acos(interior) 
  if (j2>0){
    camber$thetay[camber$ID == unique(camber$ID)[i]]   = -acos(interior)
  }
  #Depending on initial angle of elbow align either pt 6 or 7 or straighten the vector between 6 and 7
  if(camber$Angle.True[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 6] == 125){
    i2 = 0
    j2 = camber$y_adj2[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 6]
    k2 = camber$z_adj2[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 6]
    
    dotprod = abs((i2*0)+(j2*1)+(k2*0))
    len1    = sqrt(i2^2+j2^2+k2^2)
    len2    = sqrt(0^2+1^2+0^2)
    interior = dotprod/(len1*len2)
    
    camber$thetay[camber$ID == unique(camber$ID)[i]]   =  acos(interior) 
    if (j2>0){
      camber$thetay[camber$ID == unique(camber$ID)[i]]   = -acos(interior)}
  }
  
  if(camber$Angle.True[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 6] == 155){
    i2 = 0
    j2 = camber$y_adj2[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 7] 
    k2 = camber$z_adj2[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 7] 
    
    dotprod = abs((i2*0)+(j2*0)+(k2*1))
    len1    = sqrt(i2^2+j2^2+k2^2)
    len2    = sqrt(0^2+1^2+0^2)
    interior = dotprod/(len1*len2)
    
    camber$thetay[camber$ID == unique(camber$ID)[i]]   =  acos(interior) + pi 
    if (j2>0){
      camber$thetay[camber$ID == unique(camber$ID)[i]]   = -acos(interior)}
  }
  
  if(camber$ID[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 6] == "Gull2"){
    i2 = 0
    j2 = camber$y_adj2[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 7] 
    k2 = camber$z_adj2[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 7] 
    
    dotprod = abs((i2*0)+(j2*1)+(k2*0))
    len1    = sqrt(i2^2+j2^2+k2^2)
    len2    = sqrt(0^2+1^2+0^2)
    interior = dotprod/(len1*len2)
    
    camber$thetay[camber$ID == unique(camber$ID)[i]]   =  acos(interior) 
  }
  
  if(camber$ID[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 6] == "Gull1"){
    i2 = 0
    j2 = camber$y_adj2[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 6] - camber$y_adj2[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 7]
    k2 = camber$z_adj2[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 6] - camber$z_adj2[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 7]
    
    dotprod = abs((i2*0)+(j2*1)+(k2*0))
    len1    = sqrt(i2^2+j2^2+k2^2)
    len2    = sqrt(0^2+1^2+0^2)
    interior = dotprod/(len1*len2)
    
    camber$thetay[camber$ID == unique(camber$ID)[i]]   =  acos(interior)
  }
  
  if(camber$ID[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 6] == "16_1080"){
    i2 = 0
    j2 = camber$y_adj2[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 6] 
    k2 = camber$z_adj2[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 6] 
    
    dotprod = abs((i2*0)+(j2*1)+(k2*0))
    len1    = sqrt(i2^2+j2^2+k2^2)
    len2    = sqrt(0^2+1^2+0^2)
    interior = dotprod/(len1*len2)
    
    camber$thetay[camber$ID == unique(camber$ID)[i]]   =  acos(interior) 
  }
  
    
}
#-Rotate about y by the calculated angle
camber$x_adj3 = camber$x_adj2
camber$y_adj3 = cos(camber$thetay)*camber$y_adj2-sin(camber$thetay)*camber$z_adj2
camber$z_adj3 = sin(camber$thetay)*camber$y_adj2+cos(camber$thetay)*camber$z_adj2

#Step 5: Make humerus the beginning location
for (i in 1:length(unique(camber$ID))){
  camber$humx[camber$ID == unique(camber$ID)[i]] = camber$x_adj3[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 1]
  camber$humy[camber$ID == unique(camber$ID)[i]] = camber$y_adj3[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 1]
  camber$humz[camber$ID == unique(camber$ID)[i]] = camber$z_adj3[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 1]
}
camber$x_adj4 = camber$x_adj3 - camber$humx
camber$y_adj4 = camber$y_adj3 - camber$humy
camber$z_adj4 = camber$z_adj3 - camber$humz

#- Step 6: Flip about z 
camber$x_final = cos(-pi)*camber$x_adj4-sin(-pi)*camber$y_adj4
camber$y_final = sin(-pi)*camber$x_adj4+cos(-pi)*camber$y_adj4
camber$z_final = camber$z_adj4

## ------ FOR BETA CALCULATION:    ONLY FLATTEN OUT THE WINGS -----

#- Step 7: Calculate angle between x axis and 1->2 (rotate about z)
for (i in 1:length(unique(camber$ID))){  
  #no subtraction as pt1 = 0 - project on z-y plane to calculate angle
  i2 = camber$x_final[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 2] 
  j2 = camber$y_final[camber$ID == unique(camber$ID)[i] & camber$Pt.No == 2] 
  k2 = 0
  
  dotprod = abs((i2*1)+(j2*0)+(k2*0))
  len1    = sqrt(i2^2+j2^2+k2^2)
  len2    = sqrt(1^2+0^2+0^2)
  interior = dotprod/(len1*len2)
  camber$thetaflat[camber$ID == unique(camber$ID)[i]]   = acos(interior)
  
  if (j2 < 0){
    camber$thetaflat[camber$ID == unique(camber$ID)[i]] = - acos(interior)
  }
  
}

#- Rotate about y by the calculated angle
camber$x_flat = cos(camber$thetaflat)*camber$x_final-sin(camber$thetaflat)*camber$y_final
camber$y_flat = sin(camber$thetaflat)*camber$x_final+cos(camber$thetaflat)*camber$y_final
camber$z_flat = camber$z_final


## ------- Fit Function to Leading Edge -------
curve_results        <- data.frame(matrix(nrow = 14, ncol = 11))
names(curve_results) <- c("WingID","Angle.True","a","b","c","a_ci_up","a_ci_down","b_ci_up","b_ci_down","c_ci_up","c_ci_down")

fit     <-list()
fit_sum <-list()
for (j in 1:length(unique(camber$ID))){
  new_dat = subset(camber, ID %in%  c(Wing_Number[j]))
  
  curve_results$WingID[j]  <- Wing_Number[j]
  curve_results$Wing.Eff.Span[j] <- min(new_dat$x_final)-max(new_dat$x_final)
  
  fit[[j]]  <- lm(y_final~poly(x_final,2), data = new_dat)
  fit_sum[[j]]   <- summary(fit[[j]])
  
  curve_results$a[j]  <- fit_sum[[j]]$coefficients[3]
  curve_results$b[j]  <- fit_sum[[j]]$coefficients[2]
  curve_results$c[j]  <- fit_sum[[j]]$coefficients[1]
  
  curve_results$a_ci_up[j]     <- confint(fit[[j]], level = 0.95)[3]
  curve_results$a_ci_down[j]   <- confint(fit[[j]], level = 0.95)[6]
  
  curve_results$b_ci_up[j]     <- confint(fit[[j]], level = 0.95)[2]
  curve_results$b_ci_down[j]   <- confint(fit[[j]], level = 0.95)[5]
  
  curve_results$c_ci_up[j]     <- confint(fit[[j]], level = 0.95)[1]
  curve_results$c_ci_down[j]   <- confint(fit[[j]], level = 0.95)[4]
  
  curve_results$beta[j]     <- (max(new_dat$y_flat)-min(new_dat$y_flat))/abs((new_dat$x_flat[2]))*100 #highest point divied by effective wingspan
  curve_results$beta_loc[j] <- (abs(new_dat$x_flat)[which(new_dat$y_flat == max(new_dat$y_flat))]/abs(new_dat$x_flat[2]))*100
}

betamodel <- lm(curve_results$beta_loc~curve_results$Angle.True)
curve_results$Angle.True      <- wtangles$elbow.angle[match(curve_results$WingID, wtangles$WingID)]
curve_results$vertex_loc      <- curve_results$b/(2*curve_results$a)
curve_results$vertex          <- (curve_results$a*curve_results$vertex_loc^2)+(curve_results$b*curve_results$vertex_loc) + (curve_results$c)

camber$Angle.True             <- wtangles$elbow.angle[match(camber$ID, wtangles$WingID)]
max_results$beta              <- curve_results$beta[match(max_results$WingID, curve_results$WingID)]
max_results$beta_loc          <- curve_results$beta_loc[match(max_results$WingID, curve_results$WingID)]

#------ Remove wings without wind tunnel results ------
curve_results                 <- curve_results[-c(which(curve_results$WingID =="17_0243")),]  
curve_results                 <- curve_results[-c(which(curve_results$WingID =="17_0353")),]  
camber                        <- camber[-c(which(camber$ID =="17_0243")),]  
camber                        <- camber[-c(which(camber$ID =="17_0353")),]  