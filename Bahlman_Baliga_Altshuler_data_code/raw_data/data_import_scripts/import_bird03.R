####bird 03#######
#Bird 3: 2014-08-17

# twitch data
bird03_tetanus<-read.ddf.isometric("./2014-08-17/tetanus001.ddf")
bird03_peakTetanic<-max(bird03_tetanus$Force)*2/1000 # mult by 2 because of gear ratio then divided by 1000 to get from mN to N

bird03_twitch<-read.ddf.isometric("./2014-08-17/twitch004.ddf")
twitchTiming(bird03_twitch)->bird03_twitch
bird03_TTPF<-attr(bird03_twitch,"Time_to_peak_force_secs") * 1000 # convert secs to MS



#involved in experiments 4-9 and 11-13
#import relevant baselines
#WL007 baseline before 26Hz (experiment 4)
bird03_baseWL_exp04<-readAnalyzeWL("./2014-08-17/WL007.ddf")
#WL006 baseline before 28Hz (experiment 5)
bird03_baseWL_exp05<-readAnalyzeWL("./2014-08-17/WL006.ddf")
#WL008 baseline before 30Hz (experiment 6)
bird03_baseWL_exp06<-readAnalyzeWL("./2014-08-17/WL008.ddf")
#WL009 baseline before 3pulse (experiment 7)
bird03_baseWL_exp07<-readAnalyzeWL("./2014-08-17/WL009.ddf")
#WL012 baseline before 4pulse (experiment 8)
bird03_baseWL_exp08<-readAnalyzeWL("./2014-08-17/WL012.ddf")
#WL013 baseline before 3pulse (experiment 9)
bird03_baseWL_exp09<-readAnalyzeWL("./2014-08-17/WL013.ddf")
#WL013 baseline before 3pulse (experiment 11)
bird03_baseWL_exp11<-readAnalyzeWL("./2014-08-17/WL010.ddf")
#WL013 baseline before 3pulse (experiment 12)
bird03_baseWL_exp12<-readAnalyzeWL("./2014-08-17/WL014.ddf")
#WL013 baseline before 3pulse (experiment 13)
bird03_baseWL_exp13<-readAnalyzeWL("./2014-08-17/WL011.ddf")
# assess degradation via net power output 
# ("meanPower" is mean of net power over 3 cycles)
bird03_baseWL_exp04_meanPower<-mean(ldply(bird03_baseWL_exp04$NetPower)$V1)
bird03_baseWL_exp05_meanPower<-mean(ldply(bird03_baseWL_exp05$NetPower)$V1)
bird03_baseWL_exp06_meanPower<-mean(ldply(bird03_baseWL_exp06$NetPower)$V1)
bird03_baseWL_exp07_meanPower<-mean(ldply(bird03_baseWL_exp07$NetPower)$V1)
bird03_baseWL_exp08_meanPower<-mean(ldply(bird03_baseWL_exp08$NetPower)$V1)
bird03_baseWL_exp09_meanPower<-mean(ldply(bird03_baseWL_exp09$NetPower)$V1)
bird03_baseWL_exp11_meanPower<-mean(ldply(bird03_baseWL_exp11$NetPower)$V1)
bird03_baseWL_exp12_meanPower<-mean(ldply(bird03_baseWL_exp12$NetPower)$V1)
bird03_baseWL_exp13_meanPower<-mean(ldply(bird03_baseWL_exp13$NetPower)$V1)
bird03_baseline_degradation<-data.frame(bird="bird03",
                                        exp=c(4:9,11:13),
                                        meanPower=c(
                                          mean(ldply(bird03_baseWL_exp04$NetPower)$V1),
                                          mean(ldply(bird03_baseWL_exp05$NetPower)$V1),
                                          mean(ldply(bird03_baseWL_exp06$NetPower)$V1),
                                          mean(ldply(bird03_baseWL_exp07$NetPower)$V1),
                                          mean(ldply(bird03_baseWL_exp08$NetPower)$V1),
                                          mean(ldply(bird03_baseWL_exp09$NetPower)$V1),
                                          mean(ldply(bird03_baseWL_exp11$NetPower)$V1),
                                          mean(ldply(bird03_baseWL_exp12$NetPower)$V1),
                                          mean(ldply(bird03_baseWL_exp13$NetPower)$V1))
)
bird03_baseline_degradation$percent<-c((100*(abs(bird03_baseWL_exp04_meanPower)/max(abs(bird03_baseline_degradation$meanPower)))),
                                       (100*(abs(bird03_baseWL_exp05_meanPower)/max(abs(bird03_baseline_degradation$meanPower)))),
                                       (100*(abs(bird03_baseWL_exp06_meanPower)/max(abs(bird03_baseline_degradation$meanPower)))),
                                       (100*(abs(bird03_baseWL_exp07_meanPower)/max(abs(bird03_baseline_degradation$meanPower)))),
                                       (100*(abs(bird03_baseWL_exp08_meanPower)/max(abs(bird03_baseline_degradation$meanPower)))),
                                       (100*(abs(bird03_baseWL_exp09_meanPower)/max(abs(bird03_baseline_degradation$meanPower)))),
                                       (100*(abs(bird03_baseWL_exp11_meanPower)/max(abs(bird03_baseline_degradation$meanPower)))),
                                       (100*(abs(bird03_baseWL_exp12_meanPower)/max(abs(bird03_baseline_degradation$meanPower)))),
                                       (100*(abs(bird03_baseWL_exp13_meanPower)/max(abs(bird03_baseline_degradation$meanPower))))
)
bird03_baseline_degradation$time<-c(file.info("./2014-08-17/WL007.ddf")$mtime,
                                    file.info("./2014-08-17/WL006.ddf")$mtime,
                                    file.info("./2014-08-17/WL008.ddf")$mtime,
                                    file.info("./2014-08-17/WL009.ddf")$mtime,
                                    file.info("./2014-08-17/WL012.ddf")$mtime,
                                    file.info("./2014-08-17/WL013.ddf")$mtime,
                                    file.info("./2014-08-17/WL010.ddf")$mtime,
                                    file.info("./2014-08-17/WL014.ddf")$mtime,
                                    file.info("./2014-08-17/WL011.ddf")$mtime)
bird03_baseline_degradation<-bird03_baseline_degradation[order(bird03_baseline_degradation$time),]
bird03_baseline_degradation$corFac<-(max(bird03_baseline_degradation$meanPower))/bird03_baseline_degradation$meanPower
plot(bird03_baseline_degradation$time,bird03_baseline_degradation$percent,pch=19,ylim=c(0,100))
lines(bird03_baseline_degradation$time,bird03_baseline_degradation$percent)
abline(h=80,lty=2,col='red')


#####import trials#####
#Bird 3: 2014-08-17
#import experiment data
#experiment 4
bird03_exp04_path<-"./2014-08-17/vary phase 3 pulse/"
bird03_exp04_list<-WLfileinfo(bird03_exp04_path)
bird03_exp04_files<-rownames(bird03_exp04_list)
exp04_bird03_alltrials<-NULL
for(i in 1:length(bird03_exp04_files)){
  exp04_bird03_alltrials[[i]]<-readAnalyzeWL(bird03_exp04_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird03_exp04_df<-WLtrialsSummary(bird03_exp04_list,exp04_bird03_alltrials)
bird03_exp04_df$phase<-c(-20,10,-10,-30,0,-20,-15,-25,-35,-20) 
bird03_exp04_df$massspecWork=bird03_exp04_df$meanWork/bird03_pecmass
bird03_exp04_df$massspecPower=bird03_exp04_df$meanPower/bird03_pecmass
bird03_exp04_df_tc<-timeCorrect(bird03_exp04_df)
exp04_bird03_summary<-data.frame(experiment="exp04",bird="bird03",
                                 exp_bird="exp04_bird03",pec_mass=bird03_pecmass,
                                 bird03_exp04_df_tc)


#experiment 5
bird03_exp05_path<-"./2014-08-17/vary phase 4 pulse/"
bird03_exp05_list<-WLfileinfo(bird03_exp05_path)
bird03_exp05_files<-rownames(bird03_exp05_list)
exp05_bird03_alltrials<-NULL
for(i in 1:length(bird03_exp05_files)){
  exp05_bird03_alltrials[[i]]<-readAnalyzeWL(bird03_exp05_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird03_exp05_df<-WLtrialsSummary(bird03_exp05_list,exp05_bird03_alltrials)
bird03_exp05_df$phase<-c(-20,10,-10,-30,-20,0,-15,-25,-35,-20) 
bird03_exp05_df$massspecWork=bird03_exp05_df$meanWork/bird03_pecmass
bird03_exp05_df$massspecPower=bird03_exp05_df$meanPower/bird03_pecmass
bird03_exp05_df_tc<-timeCorrect(bird03_exp05_df)
exp05_bird03_summary<-data.frame(experiment="exp05",bird="bird03",
                                 exp_bird="exp05_bird03",pec_mass=bird03_pecmass,
                                 bird03_exp05_df_tc)


#experiment 6
bird03_exp06_path<-"./2014-08-17/vary phase 5 pulse/"
bird03_exp06_list<-WLfileinfo(bird03_exp06_path)
bird03_exp06_files<-rownames(bird03_exp06_list)
exp06_bird03_alltrials<-NULL
for(i in 1:length(bird03_exp06_files)){
  exp06_bird03_alltrials[[i]]<-readAnalyzeWL(bird03_exp06_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird03_exp06_df<-WLtrialsSummary(bird03_exp06_list,exp06_bird03_alltrials)
bird03_exp06_df$phase<-c(-20,10,-10,-30,0,-20,-15,-25,-35,-20) 
bird03_exp06_df$massspecWork=bird03_exp06_df$meanWork/bird03_pecmass
bird03_exp06_df$massspecPower=bird03_exp06_df$meanPower/bird03_pecmass
bird03_exp06_df_tc<-timeCorrect(bird03_exp06_df)
exp06_bird03_summary<-data.frame(experiment="exp06",bird="bird03",
                                 exp_bird="exp06_bird03",pec_mass=bird03_pecmass,
                                 bird03_exp06_df_tc)

#experiment 7
bird03_exp07_path<-"./2014-08-17/vary amp-30% 3 pulse/"
bird03_exp07_list<-WLfileinfo(bird03_exp07_path)
bird03_exp07_files<-rownames(bird03_exp07_list)
exp07_bird03_alltrials<-NULL
for(i in 1:length(bird03_exp07_files)){
  exp07_bird03_alltrials[[i]]<-readAnalyzeWL(bird03_exp07_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird03_exp07_df<-WLtrialsSummary(bird03_exp07_list,exp07_bird03_alltrials)
bird03_exp07_df$phase<-c(rep(-30,9)) 
bird03_exp07_df$massspecWork=bird03_exp07_df$meanWork/bird03_pecmass
bird03_exp07_df$massspecPower=bird03_exp07_df$meanPower/bird03_pecmass
bird03_exp07_df_tc<-timeCorrect(bird03_exp07_df)
exp07_bird03_summary<-data.frame(experiment="exp07",bird="bird03",
                                 exp_bird="exp07_bird03",pec_mass=bird03_pecmass,
                                 bird03_exp07_df_tc)


#experiment 8
bird03_exp08_path<-"./2014-08-17/vary freq -30% 3pulse/"
bird03_exp08_list<-WLfileinfo(bird03_exp08_path)
bird03_exp08_files<-rownames(bird03_exp08_list)
exp08_bird03_alltrials<-NULL
for(i in 1:length(bird03_exp08_files)){
  exp08_bird03_alltrials[[i]]<-readAnalyzeWL(bird03_exp08_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird03_exp08_df<-WLtrialsSummary(bird03_exp08_list,exp08_bird03_alltrials)
bird03_exp08_df$phase<-c(rep(-30,nrow(bird03_exp08_df))) 
bird03_exp08_df$massspecWork=bird03_exp08_df$meanWork/bird03_pecmass
bird03_exp08_df$massspecPower=bird03_exp08_df$meanPower/bird03_pecmass
bird03_exp08_df_tc<-timeCorrect(bird03_exp08_df)
exp08_bird03_summary<-data.frame(experiment="exp08",bird="bird03",
                                 exp_bird="exp08_bird03",pec_mass=bird03_pecmass,
                                 bird03_exp08_df_tc)


#experiment 9
bird03_exp09_path<-"./2014-08-17/vary amp -20% 3pulse/"
bird03_exp09_list<-WLfileinfo(bird03_exp09_path)
bird03_exp09_files<-rownames(bird03_exp09_list)
exp09_bird03_alltrials<-NULL
for(i in 1:length(bird03_exp09_files)){
  exp09_bird03_alltrials[[i]]<-readAnalyzeWL(bird03_exp09_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird03_exp09_df<-WLtrialsSummary(bird03_exp09_list,exp09_bird03_alltrials)
bird03_exp09_df$phase<-c(rep(-20,nrow(bird03_exp09_df))) 
bird03_exp09_df$massspecWork=bird03_exp09_df$meanWork/bird03_pecmass
bird03_exp09_df$massspecPower=bird03_exp09_df$meanPower/bird03_pecmass
bird03_exp09_df_tc<-timeCorrect(bird03_exp09_df)
exp09_bird03_summary<-data.frame(experiment="exp09",bird="bird03",
                                 exp_bird="exp09_bird03",pec_mass=bird03_pecmass,
                                 bird03_exp09_df_tc)


#experiment 11
bird03_exp11_path<-"./2014-08-17/vary amp -30% 4pulse/"
bird03_exp11_list<-WLfileinfo(bird03_exp11_path)
bird03_exp11_files<-rownames(bird03_exp11_list)
exp11_bird03_alltrials<-NULL
for(i in 1:length(bird03_exp11_files)){
  exp11_bird03_alltrials[[i]]<-readAnalyzeWL(bird03_exp11_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird03_exp11_df<-WLtrialsSummary(bird03_exp11_list,exp11_bird03_alltrials)
bird03_exp11_df$phase<-c(rep(-30,nrow(bird03_exp11_df))) 
bird03_exp11_df$massspecWork=bird03_exp11_df$meanWork/bird03_pecmass
bird03_exp11_df$massspecPower=bird03_exp11_df$meanPower/bird03_pecmass
bird03_exp11_df_tc<-timeCorrect(bird03_exp11_df)
exp11_bird03_summary<-data.frame(experiment="exp11",bird="bird03",
                                 exp_bird="exp11_bird03",pec_mass=bird03_pecmass,
                                 bird03_exp11_df_tc)


#experiment 12
bird03_exp12_path<-"./2014-08-17/vary freq -30% 4 pulse/"
bird03_exp12_list<-WLfileinfo(bird03_exp12_path)
bird03_exp12_files<-rownames(bird03_exp12_list)
exp12_bird03_alltrials<-NULL
for(i in 1:length(bird03_exp12_files)){
  exp12_bird03_alltrials[[i]]<-readAnalyzeWL(bird03_exp12_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird03_exp12_df<-WLtrialsSummary(bird03_exp12_list,exp12_bird03_alltrials)
bird03_exp12_df$phase<-c(rep(-30,nrow(bird03_exp12_df))) 
bird03_exp12_df$massspecWork=bird03_exp12_df$meanWork/bird03_pecmass
bird03_exp12_df$massspecPower=bird03_exp12_df$meanPower/bird03_pecmass
bird03_exp12_df_tc<-timeCorrect(bird03_exp12_df)
exp12_bird03_summary<-data.frame(experiment="exp12",bird="bird03",
                                 exp_bird="exp12_bird03",pec_mass=bird03_pecmass,
                                 bird03_exp12_df_tc)


#experiment 13
bird03_exp13_path<-"./2014-08-17/vary amp -20% 4 pulse/"
bird03_exp13_list<-WLfileinfo(bird03_exp13_path)
bird03_exp13_files<-rownames(bird03_exp13_list)
exp13_bird03_alltrials<-NULL
for(i in 1:length(bird03_exp13_files)){
  exp13_bird03_alltrials[[i]]<-readAnalyzeWL(bird03_exp13_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird03_exp13_df<-WLtrialsSummary(bird03_exp13_list,exp13_bird03_alltrials)
bird03_exp13_df$phase<-c(rep(-20,nrow(bird03_exp13_df))) 
bird03_exp13_df$massspecWork=bird03_exp13_df$meanWork/bird03_pecmass
bird03_exp13_df$massspecPower=bird03_exp13_df$meanPower/bird03_pecmass
bird03_exp13_df_tc<-timeCorrect(bird03_exp13_df)
exp13_bird03_summary<-data.frame(experiment="exp13",bird="bird03",
                                 exp_bird="exp13_bird03",pec_mass=bird03_pecmass,
                                 bird03_exp13_df_tc)


#now correct for degradation across experiments
exp04_bird03_summary$expCorWork<-exp04_bird03_summary$timeCorWork*bird03_baseline_degradation$corFac[2]
exp04_bird03_summary$expCorPower<-exp04_bird03_summary$timeCorPower*bird03_baseline_degradation$corFac[2]
exp05_bird03_summary$expCorWork<-exp05_bird03_summary$timeCorWork*bird03_baseline_degradation$corFac[1]
exp05_bird03_summary$expCorPower<-exp05_bird03_summary$timeCorPower*bird03_baseline_degradation$corFac[1]
exp06_bird03_summary$expCorWork<-exp06_bird03_summary$timeCorWork*bird03_baseline_degradation$corFac[3]
exp06_bird03_summary$expCorPower<-exp06_bird03_summary$timeCorPower*bird03_baseline_degradation$corFac[3]
exp07_bird03_summary$expCorWork<-exp07_bird03_summary$timeCorWork*bird03_baseline_degradation$corFac[4]
exp07_bird03_summary$expCorPower<-exp07_bird03_summary$timeCorPower*bird03_baseline_degradation$corFac[4]
exp08_bird03_summary$expCorWork<-exp08_bird03_summary$timeCorWork*bird03_baseline_degradation$corFac[7]
exp08_bird03_summary$expCorPower<-exp08_bird03_summary$timeCorPower*bird03_baseline_degradation$corFac[7]
exp09_bird03_summary$expCorWork<-exp09_bird03_summary$timeCorWork*bird03_baseline_degradation$corFac[8]
exp09_bird03_summary$expCorPower<-exp09_bird03_summary$timeCorPower*bird03_baseline_degradation$corFac[8]
exp11_bird03_summary$expCorWork<-exp11_bird03_summary$timeCorWork*bird03_baseline_degradation$corFac[5]
exp11_bird03_summary$expCorPower<-exp11_bird03_summary$timeCorPower*bird03_baseline_degradation$corFac[5]
exp12_bird03_summary$expCorWork<-exp12_bird03_summary$timeCorWork*bird03_baseline_degradation$corFac[9]
exp12_bird03_summary$expCorPower<-exp12_bird03_summary$timeCorPower*bird03_baseline_degradation$corFac[9]
exp13_bird03_summary$expCorWork<-exp13_bird03_summary$timeCorWork*bird03_baseline_degradation$corFac[6]
exp13_bird03_summary$expCorPower<-exp13_bird03_summary$timeCorPower*bird03_baseline_degradation$corFac[6]


#combine for a bird-specific summary object
bird03_summary<-rbind(exp04_bird03_summary,exp05_bird03_summary,exp06_bird03_summary,
                      exp07_bird03_summary,exp08_bird03_summary,exp09_bird03_summary,
                      exp11_bird03_summary,exp12_bird03_summary,exp13_bird03_summary)
#plot(bird03_summary$pulses,bird03_summary$expCorPower,col=as.factor(bird03_summary$experiment),pch=19)
#  abline(h=0)