# custom functions used for data import & processing in this study
# all written by Vikram B. Baliga (vbaliga@zoology.ubc.ca)
# 
# please also see the workloopR package for peer-reviewed code to handle muscle
# physiology experiments, developed as this study was prepped: 
# devtools::install_github("ropensci/workloopR")

library(pracma)
library(signal)
library(reshape2)


########################## read ddf files - work loops #########################
read.ddf.WL<- 
  function(filename, 
           renameCols=list(c(2,3),c("Position","Force")),
           skipCols=4:11)
{
  #first extract metadata from header
  header<-readLines(filename,n=25)
  cyc.tmp<-(strsplit(header[18],"\t",fixed=TRUE))
  Cycle_Frequency<-as.numeric(strsplit(cyc.tmp[[1]][4],",")[[1]][1])
  Amplitude<-(as.numeric(strsplit(cyc.tmp[[1]][4],",")[[1]][2]))/2
  TotalCycles<-as.numeric(strsplit(cyc.tmp[[1]][4],",")[[1]][3])-1
  p.tmp<-(strsplit(header[17],"\t",fixed=TRUE))
  Pulses<-as.numeric(strsplit(p.tmp[[1]][4],",")[[1]][4])
  Sample_Frequency<-as.numeric(sub("Sample Frequency \\(Hz\\): ","",header[2]))
  Reference_Area<-as.numeric(sub("sq. mm","",sub("Reference Area: ","",
                                                 header[3])))
  Reference_Force<-as.numeric(sub("mN","",sub("Reference Force: ","",
                                              header[4])))
  Reference_Length<-as.numeric(sub("mm","",sub("Reference Length: ","",
                                               header[5])))
  calibration<-strsplit(header[7],"\t",fixed=TRUE)[[1]][-1]
  units<-strsplit(header[8],"\t",fixed=TRUE)[[1]][-1]
  scale<-as.numeric(strsplit(header[9],"\t",fixed=TRUE)[[1]][-1]) #units per V
  offset<-as.numeric(strsplit(header[10],"\t",fixed=TRUE)[[1]][-1]) #volts
  tads<-strsplit(header[11],"\t",fixed=TRUE)[[1]][-1] #volts
  protocol<-sapply(header[17:20],function(x) paste0(x,"\n"))
  columns<-strsplit(header[25],"\t",fixed=TRUE)[[1]]
  
  #now extract data
  body<-scan(filename,skip=25,what="numeric",sep="\t",quiet=TRUE)
  bodymat<-matrix(as.numeric(body),ncol=length(columns),byrow=TRUE)
  colnames(bodymat)<-columns
  
  #rename columns, if desired
  if(!is.null(renameCols))
    columns[renameCols[[1]]]<-renameCols[[2]]
  
  #rescale the data from volts to units
  bodyUnits<-sapply(1:length(calibration), function(i)
    (bodymat[,calibration[i]]+offset[i])*scale[i])
  colnames(bodyUnits)<-columns[-c(1,length(columns))]
  bodydf<-data.frame(Time=(1:nrow(bodyUnits))/Sample_Frequency,
                     bodyUnits,
                     Stim=bodymat[,"Stim"])
  bodydf<-bodydf[,-skipCols]
  units <-c("s",units[-skipCols+1],"TTL") # +1 b/c no Sample column
  
  attr(bodydf,"Sample_Frequency")<-Sample_Frequency
  attr(bodydf,"Pulses")<-Pulses
  attr(bodydf,"Total_cycles_p2p")<-TotalCycles
  attr(bodydf,"Amplitude")<-Amplitude
  attr(bodydf,"Reference_Area")<-Reference_Area
  attr(bodydf,"Reference_Force")<-Reference_Force
  attr(bodydf,"Reference_Length")<-Reference_Length
  attr(bodydf,"Cycle_Frequency")<-Cycle_Frequency
  attr(bodydf,"units")<-units
  
  bodydf
}


########################## read ddf files - isometric ##########################
read.ddf.isometric<- 
  function(filename, 
           renameColumns=list(c(2,3),c("Position","Force")),
           deleteColumns=4:11)
{
  #first extract metadata from header
  header<-readLines(filename,n=25)
  Sample_Frequency<-as.numeric(sub("Sample Frequency \\(Hz\\): ","",
                                   header[2]))
  Reference_Area<-as.numeric(sub(" sq. mm","", sub( "Reference Area: ","",
                                                    header[3])))
  Reference_Force<-as.numeric(sub(" mN","", sub( "Reference Force: ","",
                                                 header[4] )))
  Reference_Length<-as.numeric(sub(" mm","", sub( "Reference Length: ","",
                                                  header[5])))
  calibration<-strsplit(header[7],"\t",fixed=TRUE)[[1]][-1]
  units<-strsplit(header[8],"\t",fixed=TRUE)[[1]][-1]
  scale<-as.numeric(strsplit(header[9],"\t",fixed=TRUE)[[1]][-1]) #units per V
  offset<-as.numeric(strsplit(header[10],"\t",fixed=TRUE)[[1]][-1]) #volts
  tads<-strsplit(header[11],"\t",fixed=TRUE)[[1]][-1] #volts
  protocol<-sapply(header[17:18], function(x) paste0(x,"\n"))
  
  #last line of header is the column names of the data
  columns<-strsplit(header[23],"\t",fixed=TRUE)[[1]]
  
  #now extract data
  body<-scan(filename,skip=23,what="numeric",sep="\t",quiet=TRUE)
  bodymat<-matrix(as.numeric(body),ncol=length(columns),byrow=TRUE)
  colnames(bodymat) <-columns
  
  #rename columns, if desired
  if(!is.null(renameColumns))
    columns[renameColumns[[1]]]<-renameColumns[[2]]
  
  #rescale the data from volts to units
  bodyUnits<-sapply(1:length(calibration), function(i)
    (bodymat[,calibration[i]]+offset[i])*scale[i])
  colnames(bodyUnits)<-columns[-c(1,length(columns))]
  bodydf<-data.frame(Time=(1:nrow(bodyUnits))/Sample_Frequency,
                     bodyUnits,
                     Stim=bodymat[,"Stim"])
  bodydf<-bodydf[,-deleteColumns]
  units <-c("s",units[-deleteColumns+1],"TTL") #+1 b/c no Sample column
  
  attr(bodydf,"Sample_Frequency")<-Sample_Frequency
  attr(bodydf,"Reference_Area")<-Reference_Area
  attr(bodydf,"Reference_Force")<-Reference_Force
  attr(bodydf,"Reference_Length")<-Reference_Length
  attr(bodydf,"units")<-units
  
  bodydf
  }


############################ select cycles using L0 ############################
# from a DFF file, select cycles of interest and return a new data frame
# uses L0 as reference point for defining what a "cycle" is
# x is a data frame of dff data
# bworth_freq is the critical frequency for low-pass butterworth filtering
# nups is the minimum number of increasing steps before "peaks" are reached
# keep.peaks determines which cycles are retained
# cycles must be selected in a continuous set, e.g. 3-7 not 3, 5, 7
selectCycles_lo<-function(x,bworth_freq=0.05,keep.cycles=4:6)
{
  library(signal)
  library(pracma)
  # get cycle frequency and sample frequency
  cyc_freq<-attr(x,"Cycle_Frequency")
  samp_freq<-attr(x,"Sample_Frequency")
  
  # filtering
  bworth<-butter(2,bworth_freq)
  # some artifacts; replacing leading and trailing 100 points seems safe
  smPos<-filtfilt(bworth,x$Position)
  smPos[c(1:100,(nrow(x)-100):nrow(x))]<-smPos[101]
  
  # use quarter cycle frequency to determine nups in findpeaks()
  nupz<-floor(0.25*(1/cyc_freq)*samp_freq)
  # find the peaks
  peaks<-findpeaks(smPos,nups=(nupz-1))
  
  # put it all together
  pk<-peaks
  kc<-keep.cycles
  qf<-nupz
  
  C<-data.frame(Time=x$Time[(pk[kc[1],3]+qf):(pk[kc[length(kc)],4]+qf)],
                Position=x$Position[(pk[kc[1],3]+qf):(pk[kc[length(kc)],4]+qf)],
                Force=x$Force[(pk[kc[1],3]+qf):(pk[kc[length(kc)],4]+qf)],
                Stim=x$Stim[(pk[kc[1],3]+qf):(pk[kc[length(kc)],4]+qf)])
  attr(C,"Sample_Frequency")<-attr(x,"Sample_Frequency")
  attr(C,"Cycle_Frequency")<-attr(x,"Cycle_Frequency")
  attr(C,"total_cycles")<-dim(peaks)[1]
  attr(C,"units")<-attr(x,"units")
  attr(C,"Pulses")<-attr(x,"Pulses")
  attr(C,"Amplitude")<-attr(x,"Amplitude")

  # label the cycles with letters
  numcycs<-length(keep.cycles)
  cycleslist<-NULL
  cyclesvec<-NULL
  for (i in 1:numcycs){
    cycleslist[[i]]<-rep(letters[i],length((pk[kc[i],3]+1):(pk[kc[i],4])))
  }
  cyclesvec<-melt(cycleslist)[,1]
  cyclesvec<-c(cyclesvec,numcycs)
  C$cycles<-cyclesvec
  return(C)
}


########################### select cycles peak2peak ############################
# from a DFF file, select cycles of interest and return a new data frame
# uses peak lengths as reference point for defining what a "cycle" is
# x is a data frame of dff data
# bworth_freq is the critical frequency for low-pass butterworth filtering
# nups is the minimum number of increasing steps before "peaks" are reached
# keep.peaks determines which cycles are retained
# cycles must be selected in a continuous set, e.g. 3-7 not 3, 5, 7
selectCycles_p2p<-function(x,bworth_freq=0.05,keep.cycles=3:5)
{
  library(signal)
  library(pracma)
  # get cycle frequency and sample frequency
  cyc_freq<-attr(x,"Cycle_Frequency")
  samp_freq<-attr(x,"Sample_Frequency")
  
  # filtering
  bworth<-butter(2,bworth_freq)
  # some artifacts; replacing leading and trailing 100 points seems safe
  smPos<-filtfilt(bworth,x$Position)
  smPos[c(1:100,(nrow(x)-100):nrow(x))]<-smPos[101]
  
  # use quarter cycle frequency to determine nups in findpeaks()
  nupz<-floor(0.25*(1/cyc_freq)*samp_freq)
  # find the peaks
  peaks<-findpeaks(smPos,nups=(nupz-1))
  
  # put it all together
  pk<-peaks
  kc<-keep.cycles
  d<-4 # mode; 2 = max to max; 4 = min to min

  C<-data.frame(Time=x$Time[(pk[kc[1],d]):(pk[kc[length(kc)]+1,d])],
                Position=x$Position[(pk[kc[1],d]):(pk[kc[length(kc)]+1,d])],
                Force=x$Force[(pk[kc[1],d]):(pk[kc[length(kc)]+1,d])],
                Stim=x$Stim[(pk[kc[1],d]):(pk[kc[length(kc)]+1,d])])
  attr(C,"Sample_Frequency")<-attr(x,"Sample_Frequency")
  attr(C,"Cycle_Frequency")<-attr(x,"Cycle_Frequency")
  attr(C,"Pulses")<-attr(x,"Pulses")
  attr(C,"Amplitude")<-attr(x,"Amplitude")
  attr(C,"total_cycles")<-dim(peaks)[1]-1
  attr(C,"units")<-attr(x,"units")
  
  # label the cycles with letters
  numcycs<-length(keep.cycles)
  cycleslist<-NULL
  cyclesvec<-NULL
  for (i in 1:numcycs){
    cycleslist[[i]]<-rep(letters[i],length((pk[kc[i],d]+1):(pk[(kc[i]+1),d])))
  }
  cyclesvec<-melt(cycleslist)[,1]
  cyclesvec<-c(cyclesvec,numcycs)
  C$cycles<-cyclesvec
  return(C)
}


############################ trapezoidal integration ###########################
trapezoidal.integration<-function(x,f)
{
  ### 3 checks to ensure that the arguments are numeric and of equal lengths
  # check if the variable of integration is numeric
  if (!is.numeric(x))
  {
    stop('The variable (first argument) is not numeric.')
  }
  # check if the integrand is numeric
  if (!is.numeric(f))
  {
    stop('The integrand (second argument) is not numeric.')
  }
  # check if the variable of integration and the integrand have equal lengths
  if (length(x) != length(f))
  {
    stop('The lengths of the variable and the integrand do not match.')
  }
  ### finish checks
  # obtain length of variable of integration and integrand
  n=length(x)
  # integrate using the trapezoidal rule
  integral=0.5*sum((x[2:n]-x[1:(n-1)])*(f[2:n]+f[1:(n-1)]))
  # return the definite integral
  return(integral)
}


########################### work loop data extraction ##########################
# please use one of the selectCycles_X() functions before running this
analyzeWorkLoop<-function(x,GR=2,M=-1,vel_bf=0.05){
  # GR = gear ratio of the motor arm; set it to 1 if unknown
  # M = velocity multiplier, set to -1 to make work loop positive
  # vel_bf = critical frequency for butterworth filter applied to velocity
  # selectCycles has cycle numbers.
  # first chop up the data by cycle:
  cycname<-levels(as.factor(x$cycles))
  times<-NULL
  for (i in 1:length(cycname)){
    times[[i]]<-(x[x$cycles==cycname[i],]$Time)
  }
  positions<-NULL
  for (i in 1:length(cycname)){
    positions[[i]]<-(x[x$cycles==cycname[i],]$Position)*(1/GR)
  }
  forces<-NULL
  for (i in 1:length(cycname)){
    forces[[i]]<-(x[x$cycles==cycname[i],]$Force)*(GR)
  }
  stims<-NULL
  for (i in 1:length(cycname)){
    stims[[i]]<-(x[x$cycles==cycname[i],]$Stim)
  }
  # lists to separate upper and lower parts of curves
  DS<-list()
  for (i in 1:length(cycname)){
    DS[[i]]<-1:round(length(positions[[i]])/2)
  }
  US<-list()
  for (i in 1:length(cycname)){
    US[[i]]<- (max(DS[[i]])+1):length(forces[[i]])
  }
  uppercurve<-list()
  lowercurve<-list()
  work<-list()
  velocity<-list()
  instantPower<-list()
  netPower<-list()
  
  # work is taken as the difference between the integral of the upper curve
  # and the integral of the lower curve
  # the last step is to divide by 1000 to get from mJ to J
  for (i in 1:length(cycname)){
    uppercurve[[i]]<-trapezoidal.integration(positions[[i]][rev(DS[[i]])],
                                             forces[[i]][rev(DS[[i]])])
    lowercurve[[i]]<-trapezoidal.integration(positions[[i]][US[[i]]],
                                             forces[[i]][US[[i]]])
    work[[i]]<-(uppercurve[[i]]-lowercurve[[i]])/1000
  }
 
  # velocity is the instantanous change in length (i.e. position) multiplied
  # by the sampling rate
  for (i in 1:length(cycname)){
    for (j in 1:(length(positions[[i]])-1)){
      velocity[[i]]<-M*(positions[[i]][j+1]-positions[[i]][j])*attr(x,"Sample_Frequency")
    }
  }
  # now overwrite with the actual velocity vectors
  # divide by 1000 to get from mm/sec to meters/sec
  for (i in 1:length(cycname)){
    for (j in 1:(length(positions[[i]])-1)){
      velocity[[i]][j]<-(M*(positions[[i]][j+1]-positions[[i]][j])*attr(x,"Sample_Frequency")/1000)
      }
  }
  
  # apply a butterworth filter to velocity to smooth it out a bit
  buttah<-butter(2,vel_bf)
  filt_velocity<-NULL
  for (i in 1:length(cycname)){
    filt_velocity[[i]]<-signal::filtfilt(buttah,velocity[[i]])
  }
 
  # instantaneous power is calculated as each measurement of instantaneous force
  # multiplied by its corresponding instantaneous velocity
  # then divided by 1000 to get from mW to Watts
  for (i in 1:length(cycname)){
    instantPower[[i]]<-(forces[[i]][1:(length(forces[[i]])-1)]*filt_velocity[[i]])/1000
  }
  
  # net power is simply the mean of all instantaneous power
  for (i in 1:length(cycname)){
    netPower[[i]]<-mean(instantPower[[i]])
  }
  
  # combine everything into one useful object
  X<-list(Time=times,Position=positions,Force=forces,Stim=stims,
          Velocity=velocity,Velocity_filtered=filt_velocity,Work=work,
          InstantPower=instantPower,NetPower=netPower,DS=DS,US=US,
          rawPosition=x$Position,rawForce=x$Force)
  attr(X,"Sample_Frequency")<-attr(x,"Sample_Frequency")
  attr(X,"Cycle_Frequency")<-attr(x,"Cycle_Frequency")
  attr(X,"Pulses")<-attr(x,"Pulses")
  attr(X,"Amplitude")<-attr(x,"Amplitude")
  attr(X,"Gear_ratio_used")<-GR
  attr(X,"total_cycles")<-attr(x,"total_cycles")
  attr(X,"Original_units")<-attr(x,"units")
  attr(X,"Current_units")<-c("S","mm","mN","TTL","m/s","m/s","J","W","W","n/a",
                             "n/a","mm","mN")
  return(X)
}


###################### work loop reading and data extraction ###################
## all-in-one function
readAnalyzeWL<-function(filepath,bworth_freq=0.05,keep.cycles=3:5,GR=2){
  bw<-bworth_freq
  cyclestokeep<-keep.cycles
  gearratio<-GR
  fulldata<-read.ddf.WL(filename=filepath)
  retainedcycles<-selectCycles_p2p(fulldata,bworth_freq=bw,
                                   keep.cycles=cyclestokeep)
  X<-analyzeWorkLoop(retainedcycles,GR=gearratio)
  return(X)
}


###################### file info for sequence of work loops ####################
# more details here
WLfileinfo<-function(filepath,pattern="*.ddf"){
  exp_list<-file.info(list.files(path=filepath,pattern=pattern,
                                 full.names=TRUE,recursive=TRUE))
  exp_list$exp_names<-rownames(exp_list)
  # re-order by run order, using time stamps
  exp_list<-exp_list[with(exp_list, order(as.POSIXct(mtime))), ]
  return(exp_list)
}


######################### summarize sequence of work loops #####################
# more details here
WLtrialsSummary<-function(WLlist,WLtrials){
  exp_mat<-data.frame()
  for (i in 1:dim(WLlist)[1]){
    exp_mat[i,1]<-ldply(strsplit(WLlist$exp_names,"/"))[,ncol(ldply(strsplit(WLlist$exp_names,"/")))][i] #filenames
    exp_mat[i,2]<-i #run order
    exp_mat[i,3]<-attr(WLtrials[[i]],"Cycle_Frequency")
    exp_mat[i,4]<-attr(WLtrials[[i]],"Amplitude")
    exp_mat[i,5]<-attr(WLtrials[[i]],"Pulses")
    exp_mat[i,6]<-"fillme" #fill in phases later (hard to get from files)
    exp_mat[i,7]<-WLlist$mtime[i] #time that the file was made
    exp_mat[i,8]<-mean(ldply(WLtrials[[i]]$Work)$V1) #meanWork
    exp_mat[i,9]<-mean(ldply(WLtrials[[i]]$NetPower)$V1) #meanPower
  }
  colnames(exp_mat)<-c("filename","runOrder","freq",
                       "amplitude","pulses","phase",
                       "time","meanWork","meanPower")
  return(exp_mat)
}


############################### time Correction ################################
# across a batch of trials run sequentially, correct for potential
# degradation of a muscle. Assumes that the stimulation parameters for the first
# and last trials are identical - do not use if this is not the case!! Decline 
# in power output is assumed to be a linear function of time. Accordingly, the 
# difference between the last and first trial's (absolute) power output is used
# to 'correct' trials in between, with consideration of run order and time 
# elapsed.
timeCorrect<-function(x){
  initialpower<-x$massspecPower[1]
  finalpower<-x$massspecPower[dim(x)[1]]
  powerdifference<-initialpower-finalpower
  timedifference<-x$time[dim(x)[1]]-x$time[1]
  overallcorrection<-powerdifference/as.numeric(timedifference)
  tmp<-NULL
  for (i in 1:dim(x)[1]){
    tmp[i]<-x$time[i]-x$time[1]
  }
  x$elapsed<-tmp
  x$timeCorWork<-(x$massspecPower + (tmp*overallcorrection))/x$freq
  x$timeCorPower<-x$massspecPower + (tmp*overallcorrection)
  attr(x,"Power_difference")<-powerdifference
  attr(x,"Time_difference")<-timedifference
  attr(x,"Time_correction_rate")<-overallcorrection
  return(x)
}

  

############################### twitch timing ################################
# calculates time to peak force (twitch rise time) along with twitch 50% and 
# 90% relaxation time for a given isometric twitch data file.
# the data should preferably be read via read.ddf.isometric()
twitchTiming<-function(x){
  peakForce_row<-which.max(x$Force)
  stimulation_row<-which.max(x$Stim)
  
  forcerise_10<- 0.1* (x$Force[peakForce_row] - x$Force[stimulation_row])
  forcerise_90<- 0.9* (x$Force[peakForce_row] - x$Force[stimulation_row])
  
  force10_row<-min(which(x$Force>(forcerise_10+x$Force[stimulation_row])))
  force90_row<-min(which(x$Force>(forcerise_90+x$Force[stimulation_row])))
  
  attr(x,"Time_to_peak_force_secs")<-
    x$Time[peakForce_row]-x$Time[stimulation_row] #in seconds
  
  attr(x,"Force_rise_10-90%_secs")<-
    x$Time[force90_row]-x$Time[force10_row] #in seconds
  
  fiftyp_rows<-
    which(x$Force>((x$Force[peakForce_row]+x$Force[stimulation_row])/2))
  attr(x,"50%_relaxation_time_secs")<-
    x$Time[fiftyp_rows[length(fiftyp_rows)]]-x$Time[stimulation_row] #in seconds
  
  ninetyp_rows<-
    which(x$Force>(x$Force[peakForce_row]-((x$Force[peakForce_row]-x$Force[stimulation_row])*0.9)))
  attr(x,"90%_relaxation_time_secs")<-
    x$Time[ninetyp_rows[length(ninetyp_rows)]] - x$Time[stimulation_row] #in seconds
  return(x)
}

