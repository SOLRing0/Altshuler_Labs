####bird 09#######
#Bird 9: 2014-12-17

# twitch data
#tetanic not recorded
bird09_peakTetanic<-0

bird09_twitch<-read.ddf.isometric("./2014-12-17/twitch004.ddf")
twitchTiming(bird09_twitch)->bird09_twitch
bird09_TTPF<-attr(bird09_twitch,"Time_to_peak_force_secs") * 1000 # convert secs to MS

#involved in experiments 1-3 only
#WL006 baseline before 26Hz (experiment 1)
bird09_baseWL_exp10<-readAnalyzeWL("./2014-12-17/WL011.ddf")
#WL007 baseline before 28Hz (experiment 2)
bird09_baseWL_exp14<-readAnalyzeWL("./2014-12-17/WL012.ddf")
# assess degradation via net power output 
# ("meanPower" is mean of net power over 3 cycles)
bird09_baseWL_exp10_meanPower<-mean(ldply(bird09_baseWL_exp10$NetPower)$V1)
bird09_baseWL_exp14_meanPower<-mean(ldply(bird09_baseWL_exp14$NetPower)$V1)
bird09_baseline_degradation<-data.frame(bird="bird09",
                                        exp=c(10,14),
                                        meanPower=c(
                                          mean(ldply(bird09_baseWL_exp10$NetPower)$V1),
                                          mean(ldply(bird09_baseWL_exp14$NetPower)$V1))
)
bird09_baseline_degradation$percent<-c((100*(abs(bird09_baseWL_exp10_meanPower)/max(abs(bird09_baseline_degradation$meanPower)))),
                                       (100*(abs(bird09_baseWL_exp14_meanPower)/max(abs(bird09_baseline_degradation$meanPower))))
)
bird09_baseline_degradation$time<-c(file.info("./2014-12-17/WL011.ddf")$mtime,
                                    file.info("./2014-12-17/WL012.ddf")$mtime)
bird09_baseline_degradation<-bird09_baseline_degradation[order(bird09_baseline_degradation$time),]
bird09_baseline_degradation$corFac<-(max(bird09_baseline_degradation$meanPower))/bird09_baseline_degradation$meanPower
plot(bird09_baseline_degradation$time,bird09_baseline_degradation$percent,pch=19,ylim=c(0,100))
lines(bird09_baseline_degradation$time,bird09_baseline_degradation$percent)
abline(h=80,lty=2,col='red')


#####import trials#####
#Bird 9: 2014-12-17
#import experiment data
#experiment 10
bird09_exp10_path<-"./2014-12-17/vary freq-CP 20% 3 pulse/"
bird09_exp10_list<-WLfileinfo(bird09_exp10_path)
bird09_exp10_files<-rownames(bird09_exp10_list)
exp10_bird09_alltrials<-NULL
for(i in 1:length(bird09_exp10_files)){
  exp10_bird09_alltrials[[i]]<-readAnalyzeWL(bird09_exp10_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird09_exp10_df<-WLtrialsSummary(bird09_exp10_list,exp10_bird09_alltrials)
bird09_exp10_df$phase<- -20
bird09_exp10_df$massspecWork=bird09_exp10_df$meanWork/bird09_pecmass
bird09_exp10_df$massspecPower=bird09_exp10_df$meanPower/bird09_pecmass
bird09_exp10_df_tc<-timeCorrect(bird09_exp10_df)
exp10_bird09_summary<-data.frame(experiment="exp10",bird="bird09",
                                 exp_bird="exp10_bird09",pec_mass=bird09_pecmass,
                                 bird09_exp10_df_tc)


#experiment 14
bird09_exp14_path<-"./2014-12-17/vary freq-CP 20% 4pulse/"
bird09_exp14_list<-WLfileinfo(bird09_exp14_path)
bird09_exp14_files<-rownames(bird09_exp14_list)
exp14_bird09_alltrials<-NULL
for(i in 1:length(bird09_exp14_files)){
  exp14_bird09_alltrials[[i]]<-readAnalyzeWL(bird09_exp14_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird09_exp14_df<-WLtrialsSummary(bird09_exp14_list,exp14_bird09_alltrials)
bird09_exp14_df$phase<- -20
bird09_exp14_df$massspecWork=bird09_exp14_df$meanWork/bird09_pecmass
bird09_exp14_df$massspecPower=bird09_exp14_df$meanPower/bird09_pecmass
bird09_exp14_df_tc<-timeCorrect(bird09_exp14_df)
exp14_bird09_summary<-data.frame(experiment="exp14",bird="bird09",
                                 exp_bird="exp14_bird09",pec_mass=bird09_pecmass,
                                 bird09_exp14_df_tc)


#now correct for degradation across experiments
exp10_bird09_summary$expCorWork<-exp10_bird09_summary$timeCorWork*bird09_baseline_degradation$corFac[1]
exp10_bird09_summary$expCorPower<-exp10_bird09_summary$timeCorPower*bird09_baseline_degradation$corFac[1]
exp14_bird09_summary$expCorWork<-exp14_bird09_summary$timeCorWork*bird09_baseline_degradation$corFac[2]
exp14_bird09_summary$expCorPower<-exp14_bird09_summary$timeCorPower*bird09_baseline_degradation$corFac[2]


#combine for a bird-specific summary object
bird09_summary<-rbind(exp10_bird09_summary,exp14_bird09_summary)
#plot(bird09_summary$pulses,bird09_summary$expCorPower,col=as.factor(bird09_summary$experiment),pch=19)
#  abline(h=0)