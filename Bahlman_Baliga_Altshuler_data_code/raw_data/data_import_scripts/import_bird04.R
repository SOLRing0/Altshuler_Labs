######bird 04######
#Bird 4: 2014-08-25

# twitch data
bird04_tetanus<-read.ddf.isometric("./2014-08-25/tetanus002.ddf")
bird04_peakTetanic<-max(bird04_tetanus$Force)*2/1000 # mult by 2 because of gear ratio then divided by 1000 to get from mN to N

bird04_twitch<-read.ddf.isometric("./2014-08-25/twitch004.ddf")
twitchTiming(bird04_twitch)->bird04_twitch
bird04_TTPF<-attr(bird04_twitch,"Time_to_peak_force_secs") * 1000 # convert secs to MS


#involved in experiments 1-6 only
#import relevant baselines
#WL013 baseline before 26Hz (experiment 1)
bird04_baseWL_exp01<-readAnalyzeWL("./2014-08-25/WL013.ddf")
#WL012 baseline before 28Hz (experiment 2)
bird04_baseWL_exp02<-readAnalyzeWL("./2014-08-25/WL012.ddf")
#WL014 baseline before 30Hz (experiment 3)
bird04_baseWL_exp03<-readAnalyzeWL("./2014-08-25/WL014.ddf")
#WL015 baseline before 3pulse (experiment 4)
bird04_baseWL_exp04<-readAnalyzeWL("./2014-08-25/WL015.ddf")
#WL016 baseline before 4pulse (experiment 5)
bird04_baseWL_exp05<-readAnalyzeWL("./2014-08-25/WL016.ddf")
#WL017 baseline before 5pulse (experiment 6)
bird04_baseWL_exp06<-readAnalyzeWL("./2014-08-25/WL017.ddf")
# assess degradation via net power output 
# ("meanPower" is mean of net power over 3 cycles)
bird04_baseWL_exp01_meanPower<-mean(ldply(bird04_baseWL_exp01$NetPower)$V1)
bird04_baseWL_exp02_meanPower<-mean(ldply(bird04_baseWL_exp02$NetPower)$V1)
bird04_baseWL_exp03_meanPower<-mean(ldply(bird04_baseWL_exp03$NetPower)$V1)
bird04_baseWL_exp04_meanPower<-mean(ldply(bird04_baseWL_exp04$NetPower)$V1)
bird04_baseWL_exp05_meanPower<-mean(ldply(bird04_baseWL_exp05$NetPower)$V1)
bird04_baseWL_exp06_meanPower<-mean(ldply(bird04_baseWL_exp06$NetPower)$V1)
bird04_baseline_degradation<-data.frame(bird="bird04",
                                        exp=1:6,
                                        meanPower=c(
                                          mean(ldply(bird04_baseWL_exp01$NetPower)$V1),
                                          mean(ldply(bird04_baseWL_exp02$NetPower)$V1),
                                          mean(ldply(bird04_baseWL_exp03$NetPower)$V1),
                                          mean(ldply(bird04_baseWL_exp04$NetPower)$V1),
                                          mean(ldply(bird04_baseWL_exp05$NetPower)$V1),
                                          mean(ldply(bird04_baseWL_exp06$NetPower)$V1))
)
bird04_baseline_degradation$percent<-c((100*(abs(bird04_baseWL_exp01_meanPower)/max(abs(bird04_baseline_degradation$meanPower)))),
                                       (100*(abs(bird04_baseWL_exp02_meanPower)/max(abs(bird04_baseline_degradation$meanPower)))),
                                       (100*(abs(bird04_baseWL_exp03_meanPower)/max(abs(bird04_baseline_degradation$meanPower)))),
                                       (100*(abs(bird04_baseWL_exp04_meanPower)/max(abs(bird04_baseline_degradation$meanPower)))),
                                       (100*(abs(bird04_baseWL_exp05_meanPower)/max(abs(bird04_baseline_degradation$meanPower)))),
                                       (100*(abs(bird04_baseWL_exp06_meanPower)/max(abs(bird04_baseline_degradation$meanPower))))
)
bird04_baseline_degradation$time<-c(file.info("./2014-08-25/WL013.ddf")$mtime,
                                    file.info("./2014-08-25/WL012.ddf")$mtime,
                                    file.info("./2014-08-25/WL014.ddf")$mtime,
                                    file.info("./2014-08-25/WL015.ddf")$mtime,
                                    file.info("./2014-08-25/WL016.ddf")$mtime,
                                    file.info("./2014-08-25/WL017.ddf")$mtime)
bird04_baseline_degradation<-bird04_baseline_degradation[order(bird04_baseline_degradation$time),]
bird04_baseline_degradation$corFac<-(max(bird04_baseline_degradation$meanPower))/bird04_baseline_degradation$meanPower
plot(bird04_baseline_degradation$time,bird04_baseline_degradation$percent,pch=19,ylim=c(0,100))
lines(bird04_baseline_degradation$time,bird04_baseline_degradation$percent)
abline(h=80,lty=2,col='red')


#Bird 4: 2014-08-25
#import experiment data
#26Hz (experiment 1)
bird04_exp01_path<-"./2014-08-25/vary duration 26Hz/"
bird04_exp01_list<-WLfileinfo(bird04_exp01_path)
bird04_exp01_files<-rownames(bird04_exp01_list)
exp01_bird04_alltrials<-NULL
for(i in 1:length(bird04_exp01_files)){
  exp01_bird04_alltrials[[i]]<-readAnalyzeWL(bird04_exp01_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird04_exp01_df<-WLtrialsSummary(bird04_exp01_list,exp01_bird04_alltrials)
bird04_exp01_df$phase<- -25
bird04_exp01_df$massspecWork=bird04_exp01_df$meanWork/bird04_pecmass
bird04_exp01_df$massspecPower=bird04_exp01_df$meanPower/bird04_pecmass
bird04_exp01_df_tc<-timeCorrect(bird04_exp01_df)
exp01_bird04_summary<-data.frame(experiment="exp01",bird="bird04",
                                 exp_bird="exp01_bird04",pec_mass=bird04_pecmass,
                                 bird04_exp01_df_tc)


#28Hz (experiment 2)
bird04_exp02_path<-"./2014-08-25/vary duration 28Hz/"
bird04_exp02_list<-WLfileinfo(bird04_exp02_path)
bird04_exp02_files<-rownames(bird04_exp02_list)
exp02_bird04_alltrials<-NULL
for(i in 1:length(bird04_exp02_files)){
  exp02_bird04_alltrials[[i]]<-readAnalyzeWL(bird04_exp02_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird04_exp02_df<-WLtrialsSummary(bird04_exp02_list,exp02_bird04_alltrials)
bird04_exp02_df$phase<- -25
bird04_exp02_df$massspecWork=bird04_exp02_df$meanWork/bird04_pecmass
bird04_exp02_df$massspecPower=bird04_exp02_df$meanPower/bird04_pecmass
bird04_exp02_df_tc<-timeCorrect(bird04_exp02_df)
exp02_bird04_summary<-data.frame(experiment="exp02",bird="bird04",
                                 exp_bird="exp02_bird04",pec_mass=bird04_pecmass,
                                 bird04_exp02_df_tc)


#30Hz (experiment 3)
bird04_exp03_path<-"./2014-08-25/vary duration 30Hz/"
bird04_exp03_list<-WLfileinfo(bird04_exp03_path)
bird04_exp03_files<-rownames(bird04_exp03_list)
exp03_bird04_alltrials<-NULL
for(i in 1:length(bird04_exp03_files)){
  exp03_bird04_alltrials[[i]]<-readAnalyzeWL(bird04_exp03_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird04_exp03_df<-WLtrialsSummary(bird04_exp03_list,exp03_bird04_alltrials)
bird04_exp03_df$phase<- -25
bird04_exp03_df$massspecWork=bird04_exp03_df$meanWork/bird04_pecmass
bird04_exp03_df$massspecPower=bird04_exp03_df$meanPower/bird04_pecmass
bird04_exp03_df_tc<-timeCorrect(bird04_exp03_df)
exp03_bird04_summary<-data.frame(experiment="exp03",bird="bird04",
                                 exp_bird="exp03_bird04",pec_mass=bird04_pecmass,
                                 bird04_exp03_df_tc)


#(experiment 4)
bird04_exp04_path<-"./2014-08-25/vary phase 3pulse/"
bird04_exp04_list<-WLfileinfo(bird04_exp04_path)
bird04_exp04_files<-rownames(bird04_exp04_list)
exp04_bird04_alltrials<-NULL
for(i in 1:length(bird04_exp04_files)){
  exp04_bird04_alltrials[[i]]<-readAnalyzeWL(bird04_exp04_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird04_exp04_df<-WLtrialsSummary(bird04_exp04_list,exp04_bird04_alltrials)
bird04_exp04_df$phase<-c(-20,10,-10,-30,0,-20,-15,-25,-35,-20)
bird04_exp04_df$massspecWork=bird04_exp04_df$meanWork/bird04_pecmass
bird04_exp04_df$massspecPower=bird04_exp04_df$meanPower/bird04_pecmass
bird04_exp04_df_tc<-timeCorrect(bird04_exp04_df)
exp04_bird04_summary<-data.frame(experiment="exp04",bird="bird04",
                                 exp_bird="exp04_bird04",pec_mass=bird04_pecmass,
                                 bird04_exp04_df_tc)

#(experiment 5)
bird04_exp05_path<-"./2014-08-25/vary phase 4 pulse/"
bird04_exp05_list<-WLfileinfo(bird04_exp05_path)
bird04_exp05_files<-rownames(bird04_exp05_list)
exp05_bird04_alltrials<-NULL
for(i in 1:length(bird04_exp05_files)){
  exp05_bird04_alltrials[[i]]<-readAnalyzeWL(bird04_exp05_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird04_exp05_df<-WLtrialsSummary(bird04_exp05_list,exp05_bird04_alltrials)
bird04_exp05_df$phase<-c(-20,10,-10,-30,-20,0,-15,-25,-35,-20)
bird04_exp05_df$massspecWork=bird04_exp05_df$meanWork/bird04_pecmass
bird04_exp05_df$massspecPower=bird04_exp05_df$meanPower/bird04_pecmass
bird04_exp05_df_tc<-timeCorrect(bird04_exp05_df)
exp05_bird04_summary<-data.frame(experiment="exp05",bird="bird04",
                                 exp_bird="exp05_bird04",pec_mass=bird04_pecmass,
                                 bird04_exp05_df_tc)


#(experiment 6)
bird04_exp06_path<-"./2014-08-25/vary phase 5 pulse/"
bird04_exp06_list<-WLfileinfo(bird04_exp06_path)
bird04_exp06_files<-rownames(bird04_exp06_list)
exp06_bird04_alltrials<-NULL
for(i in 1:length(bird04_exp06_files)){
  exp06_bird04_alltrials[[i]]<-readAnalyzeWL(bird04_exp06_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird04_exp06_df<-WLtrialsSummary(bird04_exp06_list,exp06_bird04_alltrials)
bird04_exp06_df$phase<-c(-20,10,-10,-30,0,-20,-15,-25,-35,-20)
bird04_exp06_df$massspecWork=bird04_exp06_df$meanWork/bird04_pecmass
bird04_exp06_df$massspecPower=bird04_exp06_df$meanPower/bird04_pecmass
bird04_exp06_df_tc<-timeCorrect(bird04_exp06_df)
exp06_bird04_summary<-data.frame(experiment="exp06",bird="bird04",
                                 exp_bird="exp06_bird04",pec_mass=bird04_pecmass,
                                 bird04_exp06_df_tc)


#now correct for degradation across experiments
exp01_bird04_summary$expCorWork<-exp01_bird04_summary$timeCorWork*bird04_baseline_degradation$corFac[2]
exp01_bird04_summary$expCorPower<-exp01_bird04_summary$timeCorPower*bird04_baseline_degradation$corFac[2]
exp02_bird04_summary$expCorWork<-exp02_bird04_summary$timeCorWork*bird04_baseline_degradation$corFac[1]
exp02_bird04_summary$expCorPower<-exp02_bird04_summary$timeCorPower*bird04_baseline_degradation$corFac[1]
exp03_bird04_summary$expCorWork<-exp03_bird04_summary$timeCorWork*bird04_baseline_degradation$corFac[3]
exp03_bird04_summary$expCorPower<-exp03_bird04_summary$timeCorPower*bird04_baseline_degradation$corFac[3]
exp04_bird04_summary$expCorWork<-exp04_bird04_summary$timeCorWork*bird04_baseline_degradation$corFac[4]
exp04_bird04_summary$expCorPower<-exp04_bird04_summary$timeCorPower*bird04_baseline_degradation$corFac[4]
exp05_bird04_summary$expCorWork<-exp05_bird04_summary$timeCorWork*bird04_baseline_degradation$corFac[5]
exp05_bird04_summary$expCorPower<-exp05_bird04_summary$timeCorPower*bird04_baseline_degradation$corFac[5]
exp06_bird04_summary$expCorWork<-exp06_bird04_summary$timeCorWork*bird04_baseline_degradation$corFac[6]
exp06_bird04_summary$expCorPower<-exp06_bird04_summary$timeCorPower*bird04_baseline_degradation$corFac[6]


#combine for a bird-specific summary object
bird04_summary<-rbind(exp01_bird04_summary,exp02_bird04_summary,exp03_bird04_summary,
                      exp04_bird04_summary,exp05_bird04_summary,exp06_bird04_summary)
#plot(bird04_summary$pulses,bird04_summary$expCorPower,col=as.factor(bird04_summary$experiment),pch=19)
#  abline(h=0)