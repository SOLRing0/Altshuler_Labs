tareremoval <- function(varying.data, N, len_data){

  # -----------------Read in support data -----------------------------
  setwd(working_directory)
  support.data  <- read.csv('SupportInfo.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )

  #---------- initialize matrices ------------
  revised.support.data        <- data.frame(matrix(nrow = len_data, ncol = 8)) #create dataframe that will allow me to verify the interpolation
  sig                         <- data.frame(matrix(nrow = len_data, ncol = 6)) 
  names(revised.support.data) <- c("Vfx_sup","Vfy_sup","Vfz_sup","Vmx_sup","Vmy_sup","Vmz_sup","type","WingID")
  names(sig)                  <- c("Vfx_sup_sig","Vfy_sup_sig","Vfz_sup_sig","Vmx_sup_sig","Vmy_sup_sig","Vmz_sup_sig")
  # --------------------- CLEAN UP INCOMING DATA ---------------------------------------------
  # ------------ 1: Zero Runs 7 and 9 which were not properly zeroed prior to run ---------------
  setwd(working_directory)
  staticdata_17mm  <- readMat("2017_05_06_support_I-17mm_Large_10msstatic1.mat")
  for (l in 1:6){support.data[support.data$Test.No == 7,(l+14)] = support.data[support.data$Test.No == 7,(l+14)] - staticdata_17mm$data.mean[1,l+1]}

  setwd(working_directory)
  staticdata_25mm  <- readMat("2017_05_06_support_I-25mm_Large_10msstatic1.mat")
  for (l in 1:6){support.data[support.data$Test.No == 9,(l+14)] = support.data[support.data$Test.No == 9,(l+14)] - staticdata_25mm$data.mean[1,l+1]}

  # ------------ 2: Update Variable data types ---------------
  for (i in 8:12){support.data[,i]  <- as.numeric(support.data[,i]);}
  for (i in 14:27){support.data[,i] <- as.numeric(support.data[,i]);}
  support.data$curr.angle           = as.integer(support.data$curr.angle); 
  
  support.data$up_or_down = "up"
  support.data$check_direction = 0
  for (i in 1:(length(support.data$WingID)-1)){
    support.data$check_direction[i] = support.data$curr.angle[i+1] - support.data$curr.angle[i]; #find the difference between
    if (support.data$check_direction[i] < 0){
      support.data$up_or_down[i] = "down"}}

  #------------ 3: Calculate the precision error from the readings ---------------
  support.data$sig_Vfx = support.data$Vfx_std/(sqrt(N/(support.data$Vfx_indsamp)))
  support.data$sig_Vfy = support.data$Vfy_std/(sqrt(N/(support.data$Vfy_indsamp)))
  support.data$sig_Vfz = support.data$Vfz_std/(sqrt(N/(support.data$Vfz_indsamp)))
  support.data$sig_Vmx = support.data$Vmx_std/(sqrt(N/(support.data$Vmx_indsamp)))
  support.data$sig_Vmy = support.data$Vmy_std/(sqrt(N/(support.data$Vmy_indsamp)))
  support.data$sig_Vmz = support.data$Vmz_std/(sqrt(N/(support.data$Vmz_indsamp)))
  
  # --------------- ITERATE THROUGH EACH ENTRY TO REMOVE THE LIFT, DRAG AND WEIGHT OF THE SUPPORT  ------------------------
  for (i in 1:len_data) {
    #--------- Step 1: Check to see if the current AOA was tested for the support data -------------
    #creates a working files for the support data that matches up to the current AOA
    workingsupp.data <- subset(support.data,(Support.Bin == varying.data$Support.Bin[i] & Grid == varying.data$Grid[i]))
    #create a list indices of where the current tested angle of attack was tested in the support data 
    # if this current angle was not tested for the support data row_exist = NA
    row_exist = which(!is.na(match(workingsupp.data$curr.angle, varying.data$curr.angle[i]))); 
    
    #--------- Step 1a: Finds the exact row which corresponds  -------------
    #checks whether or not the current position was caused by the wing moving up or down
    #will be positive if moving up and negative if moving down
    #if at the last AOA=0 then check direction will equal 0; at the maximum value check_direction with be a negative number
    check_direction = varying.data$curr.angle[i+1] - varying.data$curr.angle[i]; #find the difference between
    if (i == len_data)      {check_direction = 0}#adds a direction check for the last possible angle
    
    # adds a direction to the angle and saves the indice of the value in the workingsupp.data frame for positive ranges
    if (check_direction > 0){
      if(varying.data$curr.angle[i] >= 0){ #the equal is so that the middle zero takes the 2nd option
        row_final = row_exist[1]} else{
          row_final = row_exist[2]}}else{
            #now give downwards moving scenario
            if(varying.data$curr.angle[i] >= 0){ #the equal is so that the middle zero takes the 2nd option
              row_final = row_exist[2]} else{
                row_final = row_exist[1]}}
    
    #overwrites the previous line for only the final reading of each wing
    #allows the selction of the 3rd zero for the final reading
    if(check_direction == 0){row_final = row_exist[3]}
    #check if it is the maximum value if so it gives the only available value
    if (abs(varying.data$curr.angle[i]) == max(abs(workingsupp.data$curr.angle))){row_final = row_exist[1]}
    
    # --------------- STEP 2: INTERPOLATION IF THE ANGLE DOES NOT EXIST -----------------
    if (is.na(row_final)){ 
      
      #-----Step 2b: Finds the row number in the subset support data for the closest angle to the current tested one
      #finds the row number of the smallest closest support AOA to the current wing AOA that exists in the vector
      check_vector  = abs(workingsupp.data$curr.angle - varying.data$curr.angle[i]) # get an array of the differences between all tested values and the current angle
      closest_angle = which(check_vector == min(check_vector)) # find all minimums ## can't use which.min as this will only return one of the possible option
      closest_angle = closest_angle[which(workingsupp.data$up_or_down[closest_angle] == varying.data$up_or_down[i])] #limit to the angles that are in the same motion
      closest_angle = closest_angle[1]
      #-----Step 2c: tells if need to increment up or down in the interpolation 
      if(workingsupp.data$curr.angle[closest_angle] < varying.data$curr.angle[i] & varying.data$up_or_down[i] == "up"){
        inc = 1 }else{
          if(workingsupp.data$curr.angle[closest_angle] < varying.data$curr.angle[i] & varying.data$up_or_down[i] == "down"){
            inc = -1}else{
              if(workingsupp.data$curr.angle[closest_angle] > varying.data$curr.angle[i] & varying.data$up_or_down[i] == "up"){
                inc = -1}else{
                  inc = 1}}}
      
      ###---- Step 2d: CALCULATE ONE SIDE OF THE INTERPOLATION EQUATION (ANGLE RATIO)
      ang_top  = (workingsupp.data$curr.angle[closest_angle+inc]-varying.data$curr.angle[i])
      ang_bot  = (workingsupp.data$curr.angle[closest_angle+inc]-workingsupp.data$curr.angle[closest_angle])
      ratio    = ang_top/ang_bot
      ###---- Uncertainty of stepper ----- 
      ang_sig  = 1/360 #of a degree 100 arc seconds accuracy
      rd_ratio = sqrt((ang_sig/ang_top)^2+(ang_sig/ang_bot)^2) #relative uncertainty of this ratio calculation
      
      for (j in 1:6){
        vf_vi = (workingsupp.data[closest_angle+inc,j+14]-workingsupp.data[closest_angle,j+14])
        ### -- Uncertainty
        d_sigvf = workingsupp.data[closest_angle+inc,j+36] #absolute uncertainty of vf
        d_sigvi = workingsupp.data[closest_angle,j+36]     #absolute uncertainty of vi
        d_vf_vi = sqrt((d_sigvf)^2+(d_sigvi)^2)            #absolute uncertainty of vf - vi         
        d_vf_vi_ratio = vf_vi*ratio*(sqrt((d_vf_vi/vf_vi)^2+(rd_ratio)^2))  #absolute uncertainty of (vf-vi)*ratio
        
        ###---- FINAL INTERPOLATED RESULTS ------
        revised.support.data[i,j] = workingsupp.data[closest_angle+inc,j+14] - (vf_vi*ratio)
        revised.support.data[i,7] = "interpolated"
        revised.support.data[i,8] = varying.data$WingID[i]
        sig[i,j]                  = sqrt((d_sigvf)^2+(d_vf_vi_ratio)^2) #absolute uncertainty of the interpolated voltage
      }
      
    } else{
      #-------------------- BASIC SUBSTRACTION IF ANGLE EXISTS -----------------------------
      revised.support.data[i,1] = workingsupp.data$Vfx_mean[row_final]
      revised.support.data[i,2] = workingsupp.data$Vfy_mean[row_final]
      revised.support.data[i,3] = workingsupp.data$Vfz_mean[row_final]
      revised.support.data[i,4] = workingsupp.data$Vmx_mean[row_final]
      revised.support.data[i,5] = workingsupp.data$Vmy_mean[row_final]
      revised.support.data[i,6] = workingsupp.data$Vmz_mean[row_final]   
      revised.support.data[i,7] = "exact"
      revised.support.data[i,8] = varying.data$WingID[i]
      
      ### ABSOLUTE UNCERTAINTY OF THE ADJUSTED VOLTAGES
      sig[i,1] = workingsupp.data$sig_Vfx[row_final]
      sig[i,2] = workingsupp.data$sig_Vfy[row_final]
      sig[i,3] = workingsupp.data$sig_Vfz[row_final]
      sig[i,4] = workingsupp.data$sig_Vmx[row_final]
      sig[i,5] = workingsupp.data$sig_Vmy[row_final]
      sig[i,6] = workingsupp.data$sig_Vmz[row_final]
    }
  }

return(list(revised.support.data,sig))
}
