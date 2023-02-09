####bird 11#######
#Bird 11: 2015-03-18

# twitch data
bird11_peakTetanic<-0

bird11_twitch<-read.ddf.isometric("./2015-03-18/twitch003.ddf")
twitchTiming(bird11_twitch)->bird11_twitch
bird11_TTPF<-attr(bird11_twitch,"Time_to_peak_force_secs") * 1000 # convert secs to MS

#involved in experiments 1-3 only
#WL006 baseline before 26Hz (experiment 1)
bird11_baseWL_exp10<-readAnalyzeWL("./2015-03-18/WL009.ddf")
#WL007 baseline before 28Hz (experiment 2)
bird11_baseWL_exp14<-readAnalyzeWL("./2015-03-18/WL010.ddf")
# assess degradation via net power output 
# ("meanPower" is mean of net power over 3 cycles)
bird11_baseWL_exp10_meanPower<-mean(ldply(bird11_baseWL_exp10$NetPower)$V1)
bird11_baseWL_exp14_meanPower<-mean(ldply(bird11_baseWL_exp14$NetPower)$V1)
bird11_baseline_degradation<-data.frame(bird="bird11",
                                        exp=c(10,14),
                                        meanPower=c(
                                          mean(ldply(bird11_baseWL_exp10$NetPower)$V1),
                                          mean(ldply(bird11_baseWL_exp14$NetPower)$V1))
)
bird11_baseline_degradation$percent<-c((100*(abs(bird11_baseWL_exp10_meanPower)/max(abs(bird11_baseline_degradation$meanPower)))),
                                       (100*(abs(bird11_baseWL_exp14_meanPower)/max(abs(bird11_baseline_degradation$meanPower))))
)
bird11_baseline_degradation$time<-c(file.info("./2015-03-18/WL009.ddf")$mtime,
                                    file.info("./2015-03-18/WL010.ddf")$mtime)
bird11_baseline_degradation<-bird11_baseline_degradation[order(bird11_baseline_degradation$time),]
bird11_baseline_degradation$corFac<-1.912*(max(bird11_baseline_degradation$meanPower))/bird11_baseline_degradation$meanPower
plot(bird11_baseline_degradation$time,bird11_baseline_degradation$percent,pch=19,ylim=c(0,100))
lines(bird11_baseline_degradation$time,bird11_baseline_degradation$percent)
abline(h=80,lty=2,col='red')


#####import trials#####
#Bird 11: 2015-03-18
#import experiment data
#experiment 10
bird11_exp10_path<-"./2015-03-18/vary freq -20% 3 pulse CP/"
bird11_exp10_list<-WLfileinfo(bird11_exp10_path)
bird11_exp10_files<-rownames(bird11_exp10_list)
exp10_bird11_alltrials<-NULL
for(i in 1:length(bird11_exp10_files)){
  exp10_bird11_alltrials[[i]]<-readAnalyzeWL(bird11_exp10_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird11_exp10_df<-WLtrialsSummary(bird11_exp10_list,exp10_bird11_alltrials)
bird11_exp10_df$phase<- -20
bird11_exp10_df$massspecWork=bird11_exp10_df$meanWork/bird11_pecmass
bird11_exp10_df$massspecPower=bird11_exp10_df$meanPower/bird11_pecmass
bird11_exp10_df_tc<-timeCorrect(bird11_exp10_df)
exp10_bird11_summary<-data.frame(experiment="exp10",bird="bird11",
                                 exp_bird="exp10_bird11",pec_mass=bird11_pecmass,
                                 bird11_exp10_df_tc)


#experiment 14
bird11_exp14_path<-"./2015-03-18/vary freq -20% 4 pulse CP/"
bird11_exp14_list<-WLfileinfo(bird11_exp14_path)
bird11_exp14_files<-rownames(bird11_exp14_list)
exp14_bird11_alltrials<-NULL
for(i in 1:length(bird11_exp14_files)){
  exp14_bird11_alltrials[[i]]<-readAnalyzeWL(bird11_exp14_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird11_exp14_df<-WLtrialsSummary(bird11_exp14_list,exp14_bird11_alltrials)
bird11_exp14_df$phase<- -20
bird11_exp14_df$massspecWork=bird11_exp14_df$meanWork/bird11_pecmass
bird11_exp14_df$massspecPower=bird11_exp14_df$meanPower/bird11_pecmass
bird11_exp14_df_tc<-timeCorrect(bird11_exp14_df)
exp14_bird11_summary<-data.frame(experiment="exp14",bird="bird11",
                                 exp_bird="exp14_bird11",pec_mass=bird11_pecmass,
                                 bird11_exp14_df_tc)


#now correct for degradation across experiments
exp10_bird11_summary$expCorWork<-exp10_bird11_summary$timeCorWork*bird11_baseline_degradation$corFac[1]
exp10_bird11_summary$expCorPower<-exp10_bird11_summary$timeCorPower*bird11_baseline_degradation$corFac[1]
exp14_bird11_summary$expCorWork<-exp14_bird11_summary$timeCorWork*bird11_baseline_degradation$corFac[2]
exp14_bird11_summary$expCorPower<-exp14_bird11_summary$timeCorPower*bird11_baseline_degradation$corFac[2]


#combine for a bird-specific summary object
bird11_summary<-rbind(exp10_bird11_summary,exp14_bird11_summary)
#plot(bird11_summary$pulses,bird11_summary$expCorPower,col=as.factor(bird11_summary$experiment),pch=19)
#  abline(h=0)