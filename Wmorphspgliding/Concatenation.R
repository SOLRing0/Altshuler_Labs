# Info: This script concates all .mat files and measured values saved in excel documents
# only use if necessary to remake the data.frames used for the analysis portion

#----call libraries
library(R.matlab)
library(dirmcmc)

# ----- Data Processing for First Wind Tunnel Experiment
setwd(working_directory)

#Read in gull specific data from csv
data   <- read.csv('Gull_Info.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )

# ------ Define Constants
RunNum   = length(data$WingID)
count = 0

varying.data <- data.frame(matrix(nrow = 2500, ncol = 22))
names(varying.data) <-c("WingID","Test.No","Grid","Angle.No","temp.atm","rho.air","visc.air","curr.angle","pt_mean","Vfx_mean","Vfy_mean","Vfz_mean","Vmx_mean","Vmy_mean","Vmz_mean","pt_std","Vfx_std","Vfy_std","Vfz_std","Vmx_std","Vmy_std","Vmz_std")
foldername   = paste(data$Date.of.Run, "_", data$Angle.Set, "deg_", data$WingID, '_Test', data$Test.No, '_', data$Grid, '_10ms', sep="")

# --- Iterate through every run to save all individual data and mean values
for (i in 1:(RunNum)){
  #setsthe working directory to the foldername for each run
  setwd(working_directory)
  setwd(paste(working_directory,foldername[i],sep=""))
  
  #lists all the files inside that folder with AOA readings
  filelist <- list.files(pattern = paste(foldername[i],'_AOA_','*',sep=""))
   
  #calculates the number of angle of attacks
  AOAamount = length(filelist)
  
  # --- Iterate through every tested angle of attack for each wing
  for (j in 1:AOAamount){
  
    count = count + 1
    filename <- list.files(pattern = paste(foldername[i],'_AOA_',j,'_','_','*',sep="")) #use two underscores because the * ignores one for some reason?
    
    #import data into temp frame 
    rundata  <- readMat(filename)
    
    #save wing ID and test number
    varying.data$WingID[count]   = data$WingID[i]
    varying.data$Test.No[count]  = data$Test.No[i]
    varying.data$Grid[count]     = data$Grid[i]
    varying.data$Angle.No[count] = j #number for each angle tested
    
    #save fluid properties that varies with every angle of attack
    varying.data$temp.atm[count]   = rundata$T.atm
    varying.data$rho.air[count]    = rundata$rho.air
    varying.data$visc.air[count]   = rundata$visc.air
    varying.data$curr.angle[count] = rundata$curr.angle[j]
    
    #--- Save Pitot Tube data mean and std
    varying.data$pt_mean[count]    = rundata$data.mean[j,1]
    varying.data$pt_std[count]     = rundata$data.std[j,1]
    varying.data$pt_indsamp[count] = iact(rundata$data.raw[,1])
    
    #--- Save Force data mean and std
    varying.data$Vfx_mean[count]   = rundata$data.mean[j,2]
    varying.data$Vfx_std[count]    = rundata$data.std[j,2]
    varying.data$Vfx_indsamp[count]= iact(rundata$data.raw[,2])
    
    varying.data$Vfy_mean[count]   = rundata$data.mean[j,3]
    varying.data$Vfy_std[count]    = rundata$data.std[j,3]
    varying.data$Vfy_indsamp[count]= iact(rundata$data.raw[,3])
    
    varying.data$Vfz_mean[count]   = rundata$data.mean[j,4] 
    varying.data$Vfz_std[count]    = rundata$data.std[j,4]
    varying.data$Vfz_indsamp[count]= iact(rundata$data.raw[,4])
    
    #--- Save Moment data mean and std
    varying.data$Vmx_mean[count]   = rundata$data.mean[j,5]
    varying.data$Vmx_std[count]    = rundata$data.std[j,5]
    varying.data$Vmx_indsamp[count]= iact(rundata$data.raw[,5])
    
    varying.data$Vmy_mean[count]   = rundata$data.mean[j,6] 
    varying.data$Vmy_std[count]    = rundata$data.std[j,6]
    varying.data$Vmy_indsamp[count]= iact(rundata$data.raw[,6])
    
    varying.data$Vmz_mean[count]   = rundata$data.mean[j,7] 
    varying.data$Vmz_std[count]    = rundata$data.std[j,7]
    varying.data$Vmz_indsamp[count]= iact(rundata$data.raw[,7])
    
  }
  
  # This will save all of the one time data - at this point j should equal AOA amount and rundata will be the same
  data$P.atm[i]       = rundata$P.atm
  data$U.set[i]       = rundata$U.test
  data$P.offset[i]    = rundata$P.offset.v.fs
  data$countoffset[i] = count #will let me find the beginning of a run
}

len_data = min(grep("NA", varying.data$WingID)) - 1; # finds the location of the first NA and moves back a spot
varying.data = varying.data[1:len_data,]; #cuts the data to size
#update types
varying.data$curr.angle = as.integer(varying.data$curr.angle); 
for (i in 6:8){
  varying.data[,i]    <- as.numeric(varying.data[,i]);
}
for (i in 10:23){
  varying.data[,i]    <- as.numeric(varying.data[,i]);
}


### -------------------- STATIC DATA OFFSET -------------------------

#CREATE FILE TO SAVE ALL THE STATIC OFFSET INFO#Read in gull specific data from csv
data   <- read.csv('Gull_Info.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
static.data <- data.frame(matrix(nrow = 36, ncol = 9))
names(static.data) <-c("WingID","Test.No","Grid","temp.atm","rho.air","visc.air","pt_mean","pt_std","pt_indsamp")
RunNum   = length(data$WingID)
foldername   = paste(data$Date.of.Run, "_", data$Angle.Set, "deg_", data$WingID, '_Test', data$Test.No, '_', data$Grid, '_10ms', sep="")

# --- Iterate through every run to save all individual data and mean values
for (i in 1:(RunNum)){
  #setsthe working directory to the foldername for each run
  setwd(working_directory)
  setwd(paste(working_directory,foldername[i],sep=""))
  
  #grabs the static file
  filename <- list.files(pattern = paste('*',"static1.mat",sep=""))
  staticdata  <- readMat(filename)
  
  #save wing ID and test number
  static.data$WingID[i]   = data$WingID[i]
  static.data$Test.No[i]  = data$Test.No[i]
  static.data$Grid[i]     = data$Grid[i]
  
  #save fluid properties that varies with every angle of attack
  static.data$temp.atm[i]   = staticdata$T.atm
  static.data$rho.air[i]    = staticdata$rho.air
  static.data$visc.air[i]   = staticdata$visc.air
  
  #--- Save Pitot Tube data mean and std
  static.data$pt_mean[i]    = staticdata$data.mean[,1]
  static.data$pt_std[i]     = sd(staticdata$data.raw[,1])
  static.data$pt_indsamp[i] = iact(staticdata$data.raw[,1])
}

write.csv(data, file="ENTER YOUR WORKING DIRECTORY/GullInfo_wRuns.csv")     
write.csv(varying.data, file="ENTER YOUR WORKING DIRECTORY/ConcatenatedRunData.csv")     
write.csv(static.data, file="ENTER YOUR WORKING DIRECTORY/ConcatenatedStaticData.csv")     
