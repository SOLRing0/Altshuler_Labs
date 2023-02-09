####bird 12#######
#Bird 12: 2015-03-19

# twitch data
bird12_peakTetanic<-0

bird12_twitch<-read.ddf.isometric("./2015-03-19/twitch003.ddf")
twitchTiming(bird12_twitch)->bird12_twitch
bird12_TTPF<-attr(bird12_twitch,"Time_to_peak_force_secs") * 1000 # convert secs to MS

#involved in experiments 8,9,10,14 only
#import relevant baselines
#WL013 baseline before 26Hz (experiment 7)
bird12_baseWL_exp08<-readAnalyzeWL("./2015-03-19/WL006.ddf")
#WL012 baseline before 30Hz (experiment 11)
bird12_baseWL_exp09<-readAnalyzeWL("./2015-03-19/WL014.ddf")
#WL010 baseline before 3pulse (experiment 12)
bird12_baseWL_exp10<-readAnalyzeWL("./2015-03-19/WL008.ddf")
#WL008 baseline before 4pulse (experiment 13)
bird12_baseWL_exp14<-readAnalyzeWL("./2015-03-19/WL009.ddf")
# assess degradation via net power output 
# ("meanPower" is mean of net power over 3 cycles)
bird12_baseWL_exp08_meanPower<-mean(ldply(bird12_baseWL_exp08$NetPower)$V1)
bird12_baseWL_exp09_meanPower<-mean(ldply(bird12_baseWL_exp09$NetPower)$V1)
bird12_baseWL_exp10_meanPower<-mean(ldply(bird12_baseWL_exp10$NetPower)$V1)
bird12_baseWL_exp14_meanPower<-mean(ldply(bird12_baseWL_exp14$NetPower)$V1)
bird12_baseline_degradation<-data.frame(bird="bird12",
                                        exp=c(8,9,10,14),
                                        meanPower=c(
                                          mean(ldply(bird12_baseWL_exp08$NetPower)$V1),
                                          mean(ldply(bird12_baseWL_exp09$NetPower)$V1),
                                          mean(ldply(bird12_baseWL_exp10$NetPower)$V1),
                                          mean(ldply(bird12_baseWL_exp14$NetPower)$V1))
)
bird12_baseline_degradation$percent<-c((100*(abs(bird12_baseWL_exp08_meanPower)/max(abs(bird12_baseline_degradation$meanPower)))),
                                       (100*(abs(bird12_baseWL_exp09_meanPower)/max(abs(bird12_baseline_degradation$meanPower)))),
                                       (100*(abs(bird12_baseWL_exp10_meanPower)/max(abs(bird12_baseline_degradation$meanPower)))),
                                       (100*(abs(bird12_baseWL_exp14_meanPower)/max(abs(bird12_baseline_degradation$meanPower))))
)
bird12_baseline_degradation$time<-c(file.info("./2015-03-19/WL006.ddf")$mtime,
                                    file.info("./2015-03-19/WL014.ddf")$mtime,
                                    file.info("./2015-03-19/WL008.ddf")$mtime,
                                    file.info("./2015-03-19/WL009.ddf")$mtime)
bird12_baseline_degradation<-bird12_baseline_degradation[order(bird12_baseline_degradation$time),]
bird12_baseline_degradation$corFac<-0.92*(max(bird12_baseline_degradation$meanPower))/bird12_baseline_degradation$meanPower
plot(bird12_baseline_degradation$time,bird12_baseline_degradation$percent,pch=19,ylim=c(0,100))
lines(bird12_baseline_degradation$time,bird12_baseline_degradation$percent)
abline(h=80,lty=2,col='red')


#####import trials#####
#Bird 12: 2015-03-19
#import experiment data
#experiment 8
bird12_exp08_path<-"./2015-03-19/vary freq -30% 3 pulse CP/"
bird12_exp08_list<-WLfileinfo(bird12_exp08_path)
bird12_exp08_files<-rownames(bird12_exp08_list)
exp08_bird12_alltrials<-NULL
for(i in 1:length(bird12_exp08_files)){
  exp08_bird12_alltrials[[i]]<-readAnalyzeWL(bird12_exp08_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird12_exp08_df<-WLtrialsSummary(bird12_exp08_list,exp08_bird12_alltrials)
bird12_exp08_df$phase<- -30
bird12_exp08_df$massspecWork=bird12_exp08_df$meanWork/bird12_pecmass
bird12_exp08_df$massspecPower=bird12_exp08_df$meanPower/bird12_pecmass
bird12_exp08_df_tc<-timeCorrect(bird12_exp08_df)
exp08_bird12_summary<-data.frame(experiment="exp08",bird="bird12",
                                 exp_bird="exp08_bird12",pec_mass=bird12_pecmass,
                                 bird12_exp08_df_tc)


#experiment 9
bird12_exp09_path<-"./2015-03-19/vary amp -20% 3 pulse/"
bird12_exp09_list<-WLfileinfo(bird12_exp09_path)
bird12_exp09_files<-rownames(bird12_exp09_list)
exp09_bird12_alltrials<-NULL
for(i in 1:length(bird12_exp09_files)){
  exp09_bird12_alltrials[[i]]<-readAnalyzeWL(bird12_exp09_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird12_exp09_df<-WLtrialsSummary(bird12_exp09_list,exp09_bird12_alltrials)
bird12_exp09_df$phase<- -20
bird12_exp09_df$massspecWork=bird12_exp09_df$meanWork/bird12_pecmass
bird12_exp09_df$massspecPower=bird12_exp09_df$meanPower/bird12_pecmass
bird12_exp09_df_tc<-timeCorrect(bird12_exp09_df)
exp09_bird12_summary<-data.frame(experiment="exp09",bird="bird12",
                                 exp_bird="exp09_bird12",pec_mass=bird12_pecmass,
                                 bird12_exp09_df_tc)


#experiment 10
bird12_exp10_path<-"./2015-03-19/vary freq -20% 3 pulse CP/"
bird12_exp10_list<-WLfileinfo(bird12_exp10_path)
bird12_exp10_files<-rownames(bird12_exp10_list)
exp10_bird12_alltrials<-NULL
for(i in 1:length(bird12_exp10_files)){
  exp10_bird12_alltrials[[i]]<-readAnalyzeWL(bird12_exp10_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird12_exp10_df<-WLtrialsSummary(bird12_exp10_list,exp10_bird12_alltrials)
bird12_exp10_df$phase<- -20
bird12_exp10_df$massspecWork=bird12_exp10_df$meanWork/bird12_pecmass
bird12_exp10_df$massspecPower=bird12_exp10_df$meanPower/bird12_pecmass
bird12_exp10_df_tc<-timeCorrect(bird12_exp10_df)
exp10_bird12_summary<-data.frame(experiment="exp10",bird="bird12",
                                 exp_bird="exp10_bird12",pec_mass=bird12_pecmass,
                                 bird12_exp10_df_tc)


#experiment 14
bird12_exp14_path<-"./2015-03-19/vary freq -20% 4 pulse CP/"
bird12_exp14_list<-WLfileinfo(bird12_exp14_path)
bird12_exp14_files<-rownames(bird12_exp14_list)
exp14_bird12_alltrials<-NULL
for(i in 1:length(bird12_exp14_files)){
  exp14_bird12_alltrials[[i]]<-readAnalyzeWL(bird12_exp14_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird12_exp14_df<-WLtrialsSummary(bird12_exp14_list,exp14_bird12_alltrials)
bird12_exp14_df$phase<- -20
bird12_exp14_df$massspecWork=bird12_exp14_df$meanWork/bird12_pecmass
bird12_exp14_df$massspecPower=bird12_exp14_df$meanPower/bird12_pecmass
bird12_exp14_df_tc<-timeCorrect(bird12_exp14_df)
exp14_bird12_summary<-data.frame(experiment="exp14",bird="bird12",
                                 exp_bird="exp14_bird12",pec_mass=bird12_pecmass,
                                 bird12_exp14_df_tc)


#now correct for degradation across experiments
exp08_bird12_summary$expCorWork<-exp08_bird12_summary$timeCorWork*bird12_baseline_degradation$corFac[1]
exp08_bird12_summary$expCorPower<-exp08_bird12_summary$timeCorPower*bird12_baseline_degradation$corFac[1]
exp09_bird12_summary$expCorWork<-exp09_bird12_summary$timeCorWork*bird12_baseline_degradation$corFac[4]
exp09_bird12_summary$expCorPower<-exp09_bird12_summary$timeCorPower*bird12_baseline_degradation$corFac[4]
exp10_bird12_summary$expCorWork<-exp10_bird12_summary$timeCorWork*bird12_baseline_degradation$corFac[2]
exp10_bird12_summary$expCorPower<-exp10_bird12_summary$timeCorPower*bird12_baseline_degradation$corFac[2]
exp14_bird12_summary$expCorWork<-exp14_bird12_summary$timeCorWork*bird12_baseline_degradation$corFac[3]
exp14_bird12_summary$expCorPower<-exp14_bird12_summary$timeCorPower*bird12_baseline_degradation$corFac[3]


#combine for a bird-specific summary object
bird12_summary<-rbind(exp08_bird12_summary,exp09_bird12_summary,exp10_bird12_summary,
                      exp14_bird12_summary)
#plot(bird12_summary$pulses,bird12_summary$expCorPower,col=as.factor(bird12_summary$experiment),pch=19)
#  abline(h=0)