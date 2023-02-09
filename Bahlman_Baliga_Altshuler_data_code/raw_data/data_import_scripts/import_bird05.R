####bird 05#######
#Bird 5: 2014-09-04

# twitch data
bird05_tetanus<-read.ddf.isometric("./2014-09-04/tetanus002.ddf")
bird05_peakTetanic<-max(bird05_tetanus$Force)*2/1000 # mult by 2 because of gear ratio then divided by 1000 to get from mN to N

bird05_twitch<-read.ddf.isometric("./2014-09-04/twitch003.ddf")
twitchTiming(bird05_twitch)->bird05_twitch
bird05_TTPF<-attr(bird05_twitch,"Time_to_peak_force_secs") * 1000 # convert secs to MS



#involved in experiments 7,9,11-13 only
#import relevant baselines
#WL013 baseline before 26Hz (experiment 7)
bird05_baseWL_exp07<-readAnalyzeWL("./2014-09-04/WL015.ddf")
#WL011 baseline before 28Hz (experiment 9)
bird05_baseWL_exp09<-readAnalyzeWL("./2014-09-04/WL016.ddf")
#WL012 baseline before 30Hz (experiment 11)
bird05_baseWL_exp11<-readAnalyzeWL("./2014-09-04/WL018.ddf")
#WL010 baseline before 3pulse (experiment 12)
bird05_baseWL_exp12<-readAnalyzeWL("./2014-09-04/WL014.ddf")
#WL008 baseline before 4pulse (experiment 13)
bird05_baseWL_exp13<-readAnalyzeWL("./2014-09-04/WL021.ddf")
# assess degradation via net power output 
# ("meanPower" is mean of net power over 3 cycles)
bird05_baseWL_exp07_meanPower<-mean(ldply(bird05_baseWL_exp07$NetPower)$V1)
bird05_baseWL_exp09_meanPower<-mean(ldply(bird05_baseWL_exp09$NetPower)$V1)
bird05_baseWL_exp11_meanPower<-mean(ldply(bird05_baseWL_exp11$NetPower)$V1)
bird05_baseWL_exp12_meanPower<-mean(ldply(bird05_baseWL_exp12$NetPower)$V1)
bird05_baseWL_exp13_meanPower<-mean(ldply(bird05_baseWL_exp13$NetPower)$V1)
bird05_baseline_degradation<-data.frame(bird="bird05",
                                        exp=c(7,9,11,12,13),
                                        meanPower=c(
                                          mean(ldply(bird05_baseWL_exp07$NetPower)$V1),
                                          mean(ldply(bird05_baseWL_exp09$NetPower)$V1),
                                          mean(ldply(bird05_baseWL_exp11$NetPower)$V1),
                                          mean(ldply(bird05_baseWL_exp12$NetPower)$V1),
                                          mean(ldply(bird05_baseWL_exp13$NetPower)$V1))
)
bird05_baseline_degradation$percent<-c((100*(abs(bird05_baseWL_exp07_meanPower)/max(abs(bird05_baseline_degradation$meanPower)))),
                                       (100*(abs(bird05_baseWL_exp09_meanPower)/max(abs(bird05_baseline_degradation$meanPower)))),
                                       (100*(abs(bird05_baseWL_exp11_meanPower)/max(abs(bird05_baseline_degradation$meanPower)))),
                                       (100*(abs(bird05_baseWL_exp12_meanPower)/max(abs(bird05_baseline_degradation$meanPower)))),
                                       (100*(abs(bird05_baseWL_exp13_meanPower)/max(abs(bird05_baseline_degradation$meanPower))))
)
bird05_baseline_degradation$time<-c(file.info("./2014-09-04/WL015.ddf")$mtime,
                                    file.info("./2014-09-04/WL016.ddf")$mtime,
                                    file.info("./2014-09-04/WL018.ddf")$mtime,
                                    file.info("./2014-09-04/WL014.ddf")$mtime,
                                    file.info("./2014-09-04/WL021.ddf")$mtime)
bird05_baseline_degradation<-bird05_baseline_degradation[order(bird05_baseline_degradation$time),]
bird05_baseline_degradation$corFac<-(max(bird05_baseline_degradation$meanPower))/bird05_baseline_degradation$meanPower
plot(bird05_baseline_degradation$time,bird05_baseline_degradation$percent,pch=19,ylim=c(0,100))
lines(bird05_baseline_degradation$time,bird05_baseline_degradation$percent)
abline(h=80,lty=2,col='red')


#####import trials#####
#Bird 5: 2014-09-04
#import experiment data
#experiment 7
bird05_exp07_path<-"./2014-09-04/vary amp -30% 3 pulse/"
bird05_exp07_list<-WLfileinfo(bird05_exp07_path)
bird05_exp07_files<-rownames(bird05_exp07_list)
exp07_bird05_alltrials<-NULL
for(i in 1:length(bird05_exp07_files)){
  exp07_bird05_alltrials[[i]]<-readAnalyzeWL(bird05_exp07_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird05_exp07_df<-WLtrialsSummary(bird05_exp07_list,exp07_bird05_alltrials)
bird05_exp07_df$phase<-c(rep(-30,nrow(bird05_exp07_df)) )
bird05_exp07_df$massspecWork=bird05_exp07_df$meanWork/bird05_pecmass
bird05_exp07_df$massspecPower=bird05_exp07_df$meanPower/bird05_pecmass
bird05_exp07_df_tc<-timeCorrect(bird05_exp07_df)
exp07_bird05_summary<-data.frame(experiment="exp07",bird="bird05",
                                 exp_bird="exp07_bird05",pec_mass=bird05_pecmass,
                                 bird05_exp07_df_tc)


#experiment 9
bird05_exp09_path<-"./2014-09-04/vary amp -20% 3 pulse/"
bird05_exp09_list<-WLfileinfo(bird05_exp09_path)
bird05_exp09_files<-rownames(bird05_exp09_list)
exp09_bird05_alltrials<-NULL
for(i in 1:length(bird05_exp09_files)){
  exp09_bird05_alltrials[[i]]<-readAnalyzeWL(bird05_exp09_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird05_exp09_df<-WLtrialsSummary(bird05_exp09_list,exp09_bird05_alltrials)
bird05_exp09_df$phase<-c(rep(-20,nrow(bird05_exp09_df)) )
bird05_exp09_df$massspecWork=bird05_exp09_df$meanWork/bird05_pecmass
bird05_exp09_df$massspecPower=bird05_exp09_df$meanPower/bird05_pecmass
bird05_exp09_df_tc<-timeCorrect(bird05_exp09_df)
exp09_bird05_summary<-data.frame(experiment="exp09",bird="bird05",
                                 exp_bird="exp09_bird05",pec_mass=bird05_pecmass,
                                 bird05_exp09_df_tc)


#experiment 11
bird05_exp11_path<-"./2014-09-04/vary amp -30% 4 pulse/"
bird05_exp11_list<-WLfileinfo(bird05_exp11_path)
bird05_exp11_files<-rownames(bird05_exp11_list)
exp11_bird05_alltrials<-NULL
for(i in 1:length(bird05_exp11_files)){
  exp11_bird05_alltrials[[i]]<-readAnalyzeWL(bird05_exp11_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird05_exp11_df<-WLtrialsSummary(bird05_exp11_list,exp11_bird05_alltrials)
bird05_exp11_df$phase<-c(rep(-30,nrow(bird05_exp11_df)) )
bird05_exp11_df$massspecWork=bird05_exp11_df$meanWork/bird05_pecmass
bird05_exp11_df$massspecPower=bird05_exp11_df$meanPower/bird05_pecmass
bird05_exp11_df_tc<-timeCorrect(bird05_exp11_df)
exp11_bird05_summary<-data.frame(experiment="exp11",bird="bird05",
                                 exp_bird="exp11_bird05",pec_mass=bird05_pecmass,
                                 bird05_exp11_df_tc)


#experiment 12
bird05_exp12_path<-"./2014-09-04/vary freq -30% 4 pulse/"
bird05_exp12_list<-WLfileinfo(bird05_exp12_path)
bird05_exp12_files<-rownames(bird05_exp12_list)
exp12_bird05_alltrials<-NULL
for(i in 1:length(bird05_exp12_files)){
  exp12_bird05_alltrials[[i]]<-readAnalyzeWL(bird05_exp12_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird05_exp12_df<-WLtrialsSummary(bird05_exp12_list,exp12_bird05_alltrials)
bird05_exp12_df$phase<-c(rep(-30,nrow(bird05_exp12_df)) )
bird05_exp12_df$massspecWork=bird05_exp12_df$meanWork/bird05_pecmass
bird05_exp12_df$massspecPower=bird05_exp12_df$meanPower/bird05_pecmass
bird05_exp12_df_tc<-timeCorrect(bird05_exp12_df)
exp12_bird05_summary<-data.frame(experiment="exp12",bird="bird05",
                                 exp_bird="exp12_bird05",pec_mass=bird05_pecmass,
                                 bird05_exp12_df_tc)


#experiment 13
bird05_exp13_path<-"./2014-09-04/vary amp -20% 4 pulse/"
bird05_exp13_list<-WLfileinfo(bird05_exp13_path)
bird05_exp13_files<-rownames(bird05_exp13_list)
exp13_bird05_alltrials<-NULL
for(i in 1:length(bird05_exp13_files)){
  exp13_bird05_alltrials[[i]]<-readAnalyzeWL(bird05_exp13_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird05_exp13_df<-WLtrialsSummary(bird05_exp13_list,exp13_bird05_alltrials)
bird05_exp13_df$phase<-c(rep(-20,nrow(bird05_exp13_df)) )
bird05_exp13_df$massspecWork=bird05_exp13_df$meanWork/bird05_pecmass
bird05_exp13_df$massspecPower=bird05_exp13_df$meanPower/bird05_pecmass
bird05_exp13_df_tc<-timeCorrect(bird05_exp13_df)
exp13_bird05_summary<-data.frame(experiment="exp13",bird="bird05",
                                 exp_bird="exp13_bird05",pec_mass=bird05_pecmass,
                                 bird05_exp13_df_tc)


#now correct for degradation across experiments
exp07_bird05_summary$expCorWork<-exp07_bird05_summary$timeCorWork*bird05_baseline_degradation$corFac[2]
exp07_bird05_summary$expCorPower<-exp07_bird05_summary$timeCorPower*bird05_baseline_degradation$corFac[2]
exp09_bird05_summary$expCorWork<-exp09_bird05_summary$timeCorWork*bird05_baseline_degradation$corFac[3]
exp09_bird05_summary$expCorPower<-exp09_bird05_summary$timeCorPower*bird05_baseline_degradation$corFac[3]
exp11_bird05_summary$expCorWork<-exp11_bird05_summary$timeCorWork*bird05_baseline_degradation$corFac[4]
exp11_bird05_summary$expCorPower<-exp11_bird05_summary$timeCorPower*bird05_baseline_degradation$corFac[4]
exp12_bird05_summary$expCorWork<-exp12_bird05_summary$timeCorWork*bird05_baseline_degradation$corFac[1]
exp12_bird05_summary$expCorPower<-exp12_bird05_summary$timeCorPower*bird05_baseline_degradation$corFac[1]
exp13_bird05_summary$expCorWork<-exp13_bird05_summary$timeCorWork*bird05_baseline_degradation$corFac[5]
exp13_bird05_summary$expCorPower<-exp13_bird05_summary$timeCorPower*bird05_baseline_degradation$corFac[5]


#combine for a bird-specific summary object
bird05_summary<-rbind(exp07_bird05_summary,exp09_bird05_summary,
                      exp11_bird05_summary,exp12_bird05_summary,exp13_bird05_summary)
#plot(bird05_summary$pulses,bird05_summary$expCorPower,col=as.factor(bird05_summary$experiment),pch=19)
#  abline(h=0)