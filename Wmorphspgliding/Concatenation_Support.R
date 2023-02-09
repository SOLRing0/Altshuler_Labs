# Info: This script concates all .mat files and measured values from the support runs saved in excel documents
# only use if necessary to remake the data.frames used for the analysis portion

## -----------Support Data Concatenation
### - CAUTION IS YOU RE-RUN THIS, THE OUTPUT NEEDS TO BE UPDATED TO INCLUDE A SUPPORT.BIN AND GRID COLUMN MANUALLY
#----call libraries
library(R.matlab)
working_directory = "ENTERYOURDIRECTORYHERE"
# ----- Data Processing for First Wind Tunnel Experiment
setwd(working_directory)

count = 0

varying.data <- data.frame(matrix(nrow = 2500, ncol = 25))
names(varying.data) <-c("WingID","Test.No","Grid","Angle.No","U.test","P.atm","P.offset.v.fs","temp.atm","rho.air","visc.air","curr.angle","pt_mean","Vfx_mean","Vfy_mean","Vfz_mean","Vmx_mean","Vmy_mean","Vmz_mean","pt_std","Vfx_std","Vfy_std","Vfz_std","Vmx_std","Vmy_std","Vmz_std")

foldername <- list.files(pattern = '*')
RunNum = length(foldername)-1 # Note: removes last to not include the wing offset

# --- Iterate through every run to save all individual data and mean values
for (i in 1:RunNum){
  
  #setsthe working directory to the foldername for each run
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
    varying.data$WingID[count]      = rundata$WingNum
    varying.data$Test.No[count]     = i
    varying.data$Grid[count]        = 0 #put this in by hand
    varying.data$Angle.No[count]    = j #number for each angle tested
    # This will save all of the one time data - at this point j should equal AOA amount and rundata will be the same
    varying.data$P.atm[count]       = rundata$P.atm
    varying.data$U.test[count]       = rundata$U.test
    varying.data$P.offset.v.fs[count]    = rundata$P.offset.v.fs

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
}


write.csv(varying.data, file="ENTER YOUR WORKING DIRECTORY/SupportInfo_test.csv")     

