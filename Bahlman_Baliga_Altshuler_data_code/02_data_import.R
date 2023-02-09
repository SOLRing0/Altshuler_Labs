# please be sure to run 01_allFunctions prior to running
# anything in this script

# ensure your working directory is the raw data folder of this repo
# or use an R-project to automate this
setwd("./raw_data")
library(plyr);library(dplyr);library(viridis)


########## bird muscle masses #######
# masses of the birds.
# all in kilograms!!
bird01_pecmass<-0.00125; bird01_bodymass<-0.01523;
bird02_pecmass<-0.00107; bird02_bodymass<-0.01215;
bird03_pecmass<-0.00146; bird03_bodymass<-0.01360;
bird04_pecmass<-0.00125; bird04_bodymass<-0.01550;
bird05_pecmass<-0.00120; bird05_bodymass<-0.01350;
bird06_pecmass<-0.00119; bird06_bodymass<-0.01350;
bird07_pecmass<-0.00120; bird07_bodymass<-0.01290;
bird08_pecmass<-0.00112; bird08_bodymass<-0.01390;
bird09_pecmass<-0.00123; bird09_bodymass<-0.01540;
bird10_pecmass<-0.00127; bird10_bodymass<-0.01510;
bird11_pecmass<-0.00112; bird11_bodymass<-0.01310;
bird12_pecmass<-0.00103; bird12_bodymass<-0.01280;

birdinfo<-read.csv("./BirdInfo.csv",stringsAsFactors=TRUE)

# to calculate PCSA via equation 4 in main text
# divide muscle mass by 1000 to convert to kgs
# divide mean muscle length by 1000 to convert to m
# use 1060 kg/m3
# and then multiply by 10000 to convert m2 to cm2
# resulting value in cm2
birdinfo$PCSA<-10000*(birdinfo$MuscleMass/1000)/(1060*(birdinfo$MeanMuscleLength/1000))


################### source bird-specific data import scripts ###################
# data import scripts for each bird are stored in one folder
# scan the folder for these files and then source them:
ltpath<-"./data_import_scripts/"
ltlist<-list.files(path=ltpath,pattern="import_bird",
                   full.names=TRUE,recursive=TRUE)
for(i in 1:length(ltlist)){
  source(ltlist[i])
}


# go back to the root folder to ensure we make no changes to raw data 
setwd('..')

# collect time to peak force from twitches
birdinfo$TTPF<-c(bird01_TTPF,bird02_TTPF,bird03_TTPF,bird04_TTPF,
                 bird05_TTPF,bird06_TTPF,bird07_TTPF,bird08_TTPF,
                 bird09_TTPF,bird10_TTPF,bird11_TTPF,bird12_TTPF)

# peak tetanic force
birdinfo$peakTetanic<-c(bird01_peakTetanic,bird02_peakTetanic,bird03_peakTetanic,
                        bird04_peakTetanic,bird05_peakTetanic,bird06_peakTetanic,
                        bird07_peakTetanic,bird08_peakTetanic,bird09_peakTetanic,
                        bird10_peakTetanic,bird11_peakTetanic,bird12_peakTetanic)

# calculate peak isometric stress
# multiply by 10000 to get from cm2 to m2
# divide by 1000 to get from N to kN
# therefore multiply by 10
birdinfo$peakIsoStress<-10*birdinfo$peakTetanic/birdinfo$PCSA


################# collect all summary data into one data.frame #################
# 
allbirds_summary<-rbind(bird01_summary,bird02_summary,bird03_summary,
                        bird04_summary,bird05_summary,bird06_summary,
                        bird07_summary,bird08_summary,bird09_summary,
                        bird10_summary,bird11_summary,bird12_summary)

# combine with birdinfo
allsummarydat<-left_join(birdinfo,allbirds_summary,by="bird")

# calculcate mean velocity as a function of amplitude and frequency
# amplitude = 1/4 complete lengthening & shortening cycle, which means that 
# the actual length change in the muscle is 4*amplitude
# this will be in mm/sec
allsummarydat$shortening_velocity<-4*allsummarydat$amplitude*allsummarydat$freq
# calculate strain as a proportion of mean muscle length
allsummarydat$strain<-(2*allsummarydat$amplitude)/allsummarydat$MeanMuscleLength
# duration is 0.2 times pulses plus 3.33 ms intervals
allsummarydat$stim_duration<-((0.2* allsummarydat$pulses) + 
                                ((allsummarydat$pulses - 1) * 3.33))/1000
# duty cycle is proportion of duration over cycle duration
allsummarydat$duty_cycle<-allsummarydat$stim_duration / (1/allsummarydat$freq)



