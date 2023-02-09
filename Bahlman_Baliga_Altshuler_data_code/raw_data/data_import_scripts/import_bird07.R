####bird 07#######
#Bird 7: 2014-10-01

# twitch data
bird07_tetanus<-read.ddf.isometric("./2014-10-01/tetanus 001.ddf")
bird07_peakTetanic<-max(bird07_tetanus$Force)*2/1000 # mult by 2 because of gear ratio then divided by 1000 to get from mN to N

bird07_twitch<-read.ddf.isometric("./2014-10-01/twitch004.ddf")
twitchTiming(bird07_twitch)->bird07_twitch
bird07_TTPF<-attr(bird07_twitch,"Time_to_peak_force_secs") * 1000 # convert secs to MS

#involved in experiments 1-3 only

###### import baselines ######
#WL006 baseline before 26Hz (experiment 1)
bird07_baseWL_exp01<-readAnalyzeWL("./2014-10-01/WL006.ddf")
#WL007 baseline before 28Hz (experiment 2)
bird07_baseWL_exp02<-readAnalyzeWL("./2014-10-01/WL007.ddf")
#WL008 baseline before 30Hz (experiment 3)
bird07_baseWL_exp03<-readAnalyzeWL("./2014-10-01/WL008.ddf")
# assess degradation via net power output 
# ("meanPower" is mean of net power over 3 cycles)
bird07_baseWL_exp01_meanPower<-mean(ldply(bird07_baseWL_exp01$NetPower)$V1)
bird07_baseWL_exp02_meanPower<-mean(ldply(bird07_baseWL_exp02$NetPower)$V1)
bird07_baseWL_exp03_meanPower<-mean(ldply(bird07_baseWL_exp03$NetPower)$V1)
bird07_baseline_degradation<-data.frame(bird="bird07",
                                        exp=1:3,
                                        meanPower=c(
                                          mean(ldply(bird07_baseWL_exp01$NetPower)$V1),
                                          mean(ldply(bird07_baseWL_exp02$NetPower)$V1),
                                          mean(ldply(bird07_baseWL_exp03$NetPower)$V1))
)
bird07_baseline_degradation$percent<-c((100*(abs(bird07_baseWL_exp01_meanPower)/max(abs(bird07_baseline_degradation$meanPower)))),
                                       (100*(abs(bird07_baseWL_exp02_meanPower)/max(abs(bird07_baseline_degradation$meanPower)))),
                                       (100*(abs(bird07_baseWL_exp03_meanPower)/max(abs(bird07_baseline_degradation$meanPower))))
)
bird07_baseline_degradation$time<-c(file.info("./2014-10-01/WL006.ddf")$mtime,
                                    file.info("./2014-10-01/WL007.ddf")$mtime,
                                    file.info("./2014-10-01/WL008.ddf")$mtime)
bird07_baseline_degradation$corFac<-(max(bird07_baseline_degradation$meanPower))/bird07_baseline_degradation$meanPower
bird07_baseline_degradation<-bird07_baseline_degradation[order(bird07_baseline_degradation$time),]
plot(bird07_baseline_degradation$time,bird07_baseline_degradation$percent,pch=19,ylim=c(0,100))
lines(bird07_baseline_degradation$time,bird07_baseline_degradation$percent)
abline(h=80,lty=2,col='red')


#####import trials#####
#Bird 7: 2014-10-01
#import experiment data
#26Hz (experiment 1)
bird07_exp01_path<-"./2014-10-01/vary duration 26Hz/"
bird07_exp01_list<-WLfileinfo(bird07_exp01_path)
bird07_exp01_files<-rownames(bird07_exp01_list)
exp01_bird07_alltrials<-NULL
for(i in 1:length(bird07_exp01_files)){
  exp01_bird07_alltrials[[i]]<-readAnalyzeWL(bird07_exp01_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird07_exp01_df<-WLtrialsSummary(bird07_exp01_list,exp01_bird07_alltrials)
bird07_exp01_df$phase<- -25
bird07_exp01_df$massspecWork=bird07_exp01_df$meanWork/bird07_pecmass
bird07_exp01_df$massspecPower=bird07_exp01_df$meanPower/bird07_pecmass
bird07_exp01_df_tc<-timeCorrect(bird07_exp01_df)
exp01_bird07_summary<-data.frame(experiment="exp01",bird="bird07",
                                 exp_bird="exp01_bird07",pec_mass=bird07_pecmass,
                                 bird07_exp01_df_tc)


#28Hz (experiment 2)
bird07_exp02_path<-"./2014-10-01/vary duration 28Hz/"
bird07_exp02_list<-WLfileinfo(bird07_exp02_path)
bird07_exp02_files<-rownames(bird07_exp02_list)
exp02_bird07_alltrials<-NULL
for(i in 1:length(bird07_exp02_files)){
  exp02_bird07_alltrials[[i]]<-readAnalyzeWL(bird07_exp02_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird07_exp02_df<-WLtrialsSummary(bird07_exp02_list,exp02_bird07_alltrials)
bird07_exp02_df$phase<- -25
bird07_exp02_df$massspecWork=bird07_exp02_df$meanWork/bird07_pecmass
bird07_exp02_df$massspecPower=bird07_exp02_df$meanPower/bird07_pecmass
bird07_exp02_df_tc<-timeCorrect(bird07_exp02_df)
exp02_bird07_summary<-data.frame(experiment="exp02",bird="bird07",
                                 exp_bird="exp02_bird07",pec_mass=bird07_pecmass,
                                 bird07_exp02_df_tc)

#30Hz (experiment 3)
bird07_exp03_path<-"./2014-10-01/vary duration 30Hz/"
bird07_exp03_list<-WLfileinfo(bird07_exp03_path)
bird07_exp03_files<-rownames(bird07_exp03_list)
exp03_bird07_alltrials<-NULL
for(i in 1:length(bird07_exp03_files)){
  exp03_bird07_alltrials[[i]]<-readAnalyzeWL(bird07_exp03_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird07_exp03_df<-WLtrialsSummary(bird07_exp03_list,exp03_bird07_alltrials)
bird07_exp03_df$phase<- -25
bird07_exp03_df$massspecWork=bird07_exp03_df$meanWork/bird07_pecmass
bird07_exp03_df$massspecPower=bird07_exp03_df$meanPower/bird07_pecmass
bird07_exp03_df_tc<-timeCorrect(bird07_exp03_df)
exp03_bird07_summary<-data.frame(experiment="exp03",bird="bird07",
                                 exp_bird="exp03_bird07",pec_mass=bird07_pecmass,
                                 bird07_exp03_df_tc)


#now correct for degradation across experiments
exp01_bird07_summary$expCorWork<-exp01_bird07_summary$timeCorWork*bird07_baseline_degradation$corFac[1]
exp01_bird07_summary$expCorPower<-exp01_bird07_summary$timeCorPower*bird07_baseline_degradation$corFac[1]
exp02_bird07_summary$expCorWork<-exp02_bird07_summary$timeCorWork*bird07_baseline_degradation$corFac[2]
exp02_bird07_summary$expCorPower<-exp02_bird07_summary$timeCorPower*bird07_baseline_degradation$corFac[2]
exp03_bird07_summary$expCorWork<-exp03_bird07_summary$timeCorWork*bird07_baseline_degradation$corFac[3]
exp03_bird07_summary$expCorPower<-exp03_bird07_summary$timeCorPower*bird07_baseline_degradation$corFac[3]


#combine for a bird-specific summary object
bird07_summary<-rbind(exp01_bird07_summary,exp02_bird07_summary,exp03_bird07_summary)
#plot(bird07_summary$pulses,bird07_summary$expCorPower,col=as.factor(bird07_summary$experiment),pch=19)
#  abline(h=0)
  
  