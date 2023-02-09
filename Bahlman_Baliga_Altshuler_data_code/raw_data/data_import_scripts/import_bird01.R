########## bird 01 #######
#Bird 1: 2014-07-29

# twitch data
bird01_tetanus<-read.ddf.isometric("./2014-07-29/tetanus001.ddf")
bird01_peakTetanic<-max(bird01_tetanus$Force)*2/1000 # mult by 2 because of gear ratio then divided by 1000 to get from mN to N

bird01_twitch<-read.ddf.isometric("./2014-07-29/twitch003.ddf")
twitchTiming(bird01_twitch)->bird01_twitch
bird01_TTPF<-attr(bird01_twitch,"Time_to_peak_force_secs") * 1000 # convert secs to MS



#involved in experiments 1-3 only
#import relevant baselines
#WL006 baseline before all three experiments
bird01_baseWL_allexp<-readAnalyzeWL("./2014-07-29/WL006.ddf")
#WL007 baseline after all three experiments
bird01_baseWL_postexp<-readAnalyzeWL("./2014-07-29/WL007.ddf")
# assess degradation via net power output 
# ("meanPower" is mean of net power over 3 cycles)
bird01_baseWL_allexp_meanPower<-mean(ldply(bird01_baseWL_allexp$NetPower)$V1)
bird01_baseWL_postexp_meanPower<-mean(ldply(bird01_baseWL_postexp$NetPower)$V1)
bird01_baseline_degradation<-data.frame(bird="bird01",
                                        exp=c(0,100), # this differs from others 
                                        meanPower=c(
                                          mean(ldply(bird01_baseWL_allexp$NetPower)$V1),
                                          mean(ldply(bird01_baseWL_postexp$NetPower)$V1))
)
bird01_baseline_degradation$percent<-c((100*(abs(bird01_baseWL_allexp_meanPower)/max(abs(bird01_baseline_degradation$meanPower)))),
                                       (100*(abs(bird01_baseWL_postexp_meanPower)/max(abs(bird01_baseline_degradation$meanPower))))
)
bird01_baseline_degradation$time<-c(file.info("./2014-07-29/WL006.ddf")$mtime,
                                    file.info("./2014-07-29/WL007.ddf")$mtime)
bird01_baseline_degradation$corFac<-(max(bird01_baseline_degradation$meanPower))/bird01_baseline_degradation$meanPower
bird01_baseline_degradation<-bird01_baseline_degradation[order(bird01_baseline_degradation$time),]
plot(bird01_baseline_degradation$time,bird01_baseline_degradation$percent,pch=19,ylim=c(0,100))
lines(bird01_baseline_degradation$time,bird01_baseline_degradation$percent)
abline(h=80,lty=2,col='red')


#Bird 1: 2014-07-29
#import experiment data
#26Hz (experiment 1)
bird01_exp01_path<-"./2014-07-29/vary duration 26Hz/"
bird01_exp01_list<-WLfileinfo(bird01_exp01_path)
bird01_exp01_files<-rownames(bird01_exp01_list)
exp01_bird01_alltrials<-NULL
for(i in 1:length(bird01_exp01_files)){
  exp01_bird01_alltrials[[i]]<-readAnalyzeWL(bird01_exp01_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird01_exp01_df<-WLtrialsSummary(bird01_exp01_list,exp01_bird01_alltrials)
bird01_exp01_df$phase<- -25
bird01_exp01_df$massspecWork=bird01_exp01_df$meanWork/bird01_pecmass
bird01_exp01_df$massspecPower=bird01_exp01_df$meanPower/bird01_pecmass
bird01_exp01_df_tc<-timeCorrect(bird01_exp01_df)
exp01_bird01_summary<-data.frame(experiment="exp01",bird="bird01",
                                 exp_bird="exp01_bird01",pec_mass=bird01_pecmass,
                                 bird01_exp01_df_tc)

#28Hz (experiment 2)
bird01_exp02_path<-"./2014-07-29/vary duration 28Hz/"
bird01_exp02_list<-WLfileinfo(bird01_exp02_path)
bird01_exp02_files<-rownames(bird01_exp02_list)
exp02_bird01_alltrials<-NULL
for(i in 1:length(bird01_exp02_files)){
  exp02_bird01_alltrials[[i]]<-readAnalyzeWL(bird01_exp02_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird01_exp02_df<-WLtrialsSummary(bird01_exp02_list,exp02_bird01_alltrials)
bird01_exp02_df$phase<- -25
bird01_exp02_df$massspecWork=bird01_exp02_df$meanWork/bird01_pecmass
bird01_exp02_df$massspecPower=bird01_exp02_df$meanPower/bird01_pecmass
bird01_exp02_df_tc<-timeCorrect(bird01_exp02_df)
exp02_bird01_summary<-data.frame(experiment="exp02",bird="bird01",
                                 exp_bird="exp02_bird01",pec_mass=bird01_pecmass,
                                 bird01_exp02_df_tc)

#30Hz (experiment 3)
bird01_exp03_path<-"./2014-07-29/vary duration 30Hz/"
bird01_exp03_list<-WLfileinfo(bird01_exp03_path)
bird01_exp03_files<-rownames(bird01_exp03_list)
exp03_bird01_alltrials<-NULL
for(i in 1:length(bird01_exp03_files)){
  exp03_bird01_alltrials[[i]]<-readAnalyzeWL(bird01_exp03_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird01_exp03_df<-WLtrialsSummary(bird01_exp03_list,exp03_bird01_alltrials)
bird01_exp03_df$phase<- -25
bird01_exp03_df$massspecWork=bird01_exp03_df$meanWork/bird01_pecmass
bird01_exp03_df$massspecPower=bird01_exp03_df$meanPower/bird01_pecmass
bird01_exp03_df_tc<-timeCorrect(bird01_exp03_df)
exp03_bird01_summary<-data.frame(experiment="exp03",bird="bird01",
                                 exp_bird="exp03_bird01",pec_mass=bird01_pecmass,
                                 bird01_exp03_df_tc)


#now correct for degradation across experiments
exp01_bird01_summary$expCorWork<-exp01_bird01_summary$timeCorWork*bird01_baseline_degradation$corFac[1]
exp01_bird01_summary$expCorPower<-exp01_bird01_summary$timeCorPower*bird01_baseline_degradation$corFac[1]
exp02_bird01_summary$expCorWork<-exp02_bird01_summary$timeCorWork*mean(bird01_baseline_degradation$corFac)
exp02_bird01_summary$expCorPower<-exp02_bird01_summary$timeCorPower*mean(bird01_baseline_degradation$corFac)
exp03_bird01_summary$expCorWork<-exp03_bird01_summary$timeCorWork*bird01_baseline_degradation$corFac[2]
exp03_bird01_summary$expCorPower<-exp03_bird01_summary$timeCorPower*bird01_baseline_degradation$corFac[2]


#combine for a bird-specific summary object
bird01_summary<-rbind(exp01_bird01_summary,exp02_bird01_summary,exp03_bird01_summary)
#plot(bird01_summary$pulses,bird01_summary$expCorPower,col=as.factor(bird01_summary$experiment),pch=19)
#  abline(h=0)