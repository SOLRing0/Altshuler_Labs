####bird 02#######
#Bird 2: 2014-08-15

# twitch data
bird02_tetanus<-read.ddf.isometric("./2014-08-15/tetanus001.ddf")
bird02_peakTetanic<-max(bird02_tetanus$Force)*2/1000 # mult by 2 because of gear ratio then divided by 1000 to get from mN to N

bird02_twitch<-read.ddf.isometric("./2014-08-15/twitch004.ddf")
twitchTiming(bird02_twitch)->bird02_twitch
bird02_TTPF<-attr(bird02_twitch,"Time_to_peak_force_secs") * 1000 # convert secs to MS



#involved in experiments 4-8.11:13 only
#import relevant baselines
#WL013 baseline before 26Hz (experiment 4)
bird02_baseWL_exp04<-readAnalyzeWL("./2014-08-15/WL013.ddf")
#WL011 baseline before 28Hz (experiment 5)
bird02_baseWL_exp05<-readAnalyzeWL("./2014-08-15/WL011.ddf")
#WL012 baseline before 30Hz (experiment 6)
bird02_baseWL_exp06<-readAnalyzeWL("./2014-08-15/WL012.ddf")
#WL010 baseline before 3pulse (experiment 7)
bird02_baseWL_exp07<-readAnalyzeWL("./2014-08-15/WL010.ddf")
#WL008 baseline before 4pulse (experiment 8)
bird02_baseWL_exp08<-readAnalyzeWL("./2014-08-15/WL008.ddf")
#WL012 baseline before 30Hz (experiment 11)
bird02_baseWL_exp11<-readAnalyzeWL("./2014-08-15/WL009.ddf")
#WL010 baseline before 3pulse (experiment 12)
bird02_baseWL_exp12<-readAnalyzeWL("./2014-08-15/WL020.ddf") #WL021.ddf
#WL008 baseline before 4pulse (experiment 13)
bird02_baseWL_exp13<-readAnalyzeWL("./2014-08-15/WL020.ddf")
# assess degradation via net power output 
# ("meanPower" is mean of net power over 3 cycles)
bird02_baseWL_exp04_meanPower<-mean(ldply(bird02_baseWL_exp04$NetPower)$V1)
bird02_baseWL_exp05_meanPower<-mean(ldply(bird02_baseWL_exp05$NetPower)$V1)
bird02_baseWL_exp06_meanPower<-mean(ldply(bird02_baseWL_exp06$NetPower)$V1)
bird02_baseWL_exp07_meanPower<-mean(ldply(bird02_baseWL_exp07$NetPower)$V1)
bird02_baseWL_exp08_meanPower<-mean(ldply(bird02_baseWL_exp08$NetPower)$V1)
bird02_baseWL_exp11_meanPower<-mean(ldply(bird02_baseWL_exp11$NetPower)$V1)
bird02_baseWL_exp12_meanPower<-mean(ldply(bird02_baseWL_exp12$NetPower)$V1)
bird02_baseWL_exp13_meanPower<-mean(ldply(bird02_baseWL_exp13$NetPower)$V1)
bird02_baseline_degradation<-data.frame(bird="bird02",
                                        exp=c(4:8,11:13),
                                        meanPower=c(
                                          mean(ldply(bird02_baseWL_exp04$NetPower)$V1),
                                          mean(ldply(bird02_baseWL_exp05$NetPower)$V1),
                                          mean(ldply(bird02_baseWL_exp06$NetPower)$V1),
                                          mean(ldply(bird02_baseWL_exp07$NetPower)$V1),
                                          mean(ldply(bird02_baseWL_exp08$NetPower)$V1),
                                          mean(ldply(bird02_baseWL_exp11$NetPower)$V1),
                                          mean(ldply(bird02_baseWL_exp12$NetPower)$V1),
                                          mean(ldply(bird02_baseWL_exp13$NetPower)$V1))
)
bird02_baseline_degradation$percent<-c((100*(abs(bird02_baseWL_exp04_meanPower)/max(abs(bird02_baseline_degradation$meanPower)))),
                                       (100*(abs(bird02_baseWL_exp05_meanPower)/max(abs(bird02_baseline_degradation$meanPower)))),
                                       (100*(abs(bird02_baseWL_exp06_meanPower)/max(abs(bird02_baseline_degradation$meanPower)))),
                                       (100*(abs(bird02_baseWL_exp07_meanPower)/max(abs(bird02_baseline_degradation$meanPower)))),
                                       (100*(abs(bird02_baseWL_exp08_meanPower)/max(abs(bird02_baseline_degradation$meanPower)))),
                                       (100*(abs(bird02_baseWL_exp11_meanPower)/max(abs(bird02_baseline_degradation$meanPower)))),
                                       (100*(abs(bird02_baseWL_exp12_meanPower)/max(abs(bird02_baseline_degradation$meanPower)))),
                                       (100*(abs(bird02_baseWL_exp13_meanPower)/max(abs(bird02_baseline_degradation$meanPower))))
)
bird02_baseline_degradation$time<-c(file.info("./2014-08-15/WL013.ddf")$mtime,
                                    file.info("./2014-08-15/WL011.ddf")$mtime,
                                    file.info("./2014-08-15/WL012.ddf")$mtime,
                                    file.info("./2014-08-15/WL010.ddf")$mtime,
                                    file.info("./2014-08-15/WL008.ddf")$mtime,
                                    file.info("./2014-08-15/WL009.ddf")$mtime,
                                    file.info("./2014-08-15/WL021.ddf")$mtime,
                                    file.info("./2014-08-15/WL020.ddf")$mtime)
bird02_baseline_degradation<-bird02_baseline_degradation[order(bird02_baseline_degradation$time),]
bird02_baseline_degradation$corFac<-(max(bird02_baseline_degradation$meanPower))/bird02_baseline_degradation$meanPower
plot(bird02_baseline_degradation$time,bird02_baseline_degradation$percent,pch=19,ylim=c(0,100))
lines(bird02_baseline_degradation$time,bird02_baseline_degradation$percent)
abline(h=80,lty=2,col='red')


#####import trials#####
#Bird 2: 2014-08-15
#import experiment data
#experiment 4
bird02_exp04_path<-"./2014-08-15/vary phase  3 pulse/"
bird02_exp04_list<-WLfileinfo(bird02_exp04_path)
bird02_exp04_files<-rownames(bird02_exp04_list)
exp04_bird02_alltrials<-NULL
for(i in 1:length(bird02_exp04_files)){
  exp04_bird02_alltrials[[i]]<-readAnalyzeWL(bird02_exp04_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird02_exp04_df<-WLtrialsSummary(bird02_exp04_list,exp04_bird02_alltrials)
bird02_exp04_df$phase<-c(-20,10,-10,-30,0,-20,-15,-25,-35,-20) 
bird02_exp04_df$massspecWork=bird02_exp04_df$meanWork/bird02_pecmass
bird02_exp04_df$massspecPower=bird02_exp04_df$meanPower/bird02_pecmass
bird02_exp04_df_tc<-timeCorrect(bird02_exp04_df)
exp04_bird02_summary<-data.frame(experiment="exp04",bird="bird02",
                                 exp_bird="exp04_bird02",pec_mass=bird02_pecmass,
                                 bird02_exp04_df_tc)


#experiment 5
bird02_exp05_path<-"./2014-08-15/vary phase 4 pulse/"
bird02_exp05_list<-WLfileinfo(bird02_exp05_path)
bird02_exp05_files<-rownames(bird02_exp05_list)
exp05_bird02_alltrials<-NULL
for(i in 1:length(bird02_exp05_files)){
  exp05_bird02_alltrials[[i]]<-readAnalyzeWL(bird02_exp05_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird02_exp05_df<-WLtrialsSummary(bird02_exp05_list,exp05_bird02_alltrials)
bird02_exp05_df$phase<-c(-20,10,-10,-30,-20,0,-15,-25,-35,-20) 
bird02_exp05_df$massspecWork=bird02_exp05_df$meanWork/bird02_pecmass
bird02_exp05_df$massspecPower=bird02_exp05_df$meanPower/bird02_pecmass
bird02_exp05_df_tc<-timeCorrect(bird02_exp05_df)
exp05_bird02_summary<-data.frame(experiment="exp05",bird="bird02",
                                 exp_bird="exp05_bird02",pec_mass=bird02_pecmass,
                                 bird02_exp05_df_tc)

#experiment 6
bird02_exp06_path<-"./2014-08-15/vary phase 5 pulse/"
bird02_exp06_list<-WLfileinfo(bird02_exp06_path)
bird02_exp06_files<-rownames(bird02_exp06_list)
exp06_bird02_alltrials<-NULL
for(i in 1:length(bird02_exp06_files)){
  exp06_bird02_alltrials[[i]]<-readAnalyzeWL(bird02_exp06_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird02_exp06_df<-WLtrialsSummary(bird02_exp06_list,exp06_bird02_alltrials)
bird02_exp06_df$phase<-c(-20,10,-10,-30,0,-20,-15,-25,-35,-20) 
bird02_exp06_df$massspecWork=bird02_exp06_df$meanWork/bird02_pecmass
bird02_exp06_df$massspecPower=bird02_exp06_df$meanPower/bird02_pecmass
bird02_exp06_df_tc<-timeCorrect(bird02_exp06_df)
exp06_bird02_summary<-data.frame(experiment="exp06",bird="bird02",
                                 exp_bird="exp06_bird02",pec_mass=bird02_pecmass,
                                 bird02_exp06_df_tc)


#experiment 7
bird02_exp07_path<-"./2014-08-15/vary amp-3 pulse/"
bird02_exp07_list<-WLfileinfo(bird02_exp07_path)
bird02_exp07_files<-rownames(bird02_exp07_list)
exp07_bird02_alltrials<-NULL
for(i in 1:length(bird02_exp07_files)){
  exp07_bird02_alltrials[[i]]<-readAnalyzeWL(bird02_exp07_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird02_exp07_df<-WLtrialsSummary(bird02_exp07_list,exp07_bird02_alltrials)
bird02_exp07_df$phase<- -30
bird02_exp07_df$massspecWork=bird02_exp07_df$meanWork/bird02_pecmass
bird02_exp07_df$massspecPower=bird02_exp07_df$meanPower/bird02_pecmass
bird02_exp07_df_tc<-timeCorrect(bird02_exp07_df)
exp07_bird02_summary<-data.frame(experiment="exp07",bird="bird02",
                                 exp_bird="exp07_bird02",pec_mass=bird02_pecmass,
                                 bird02_exp07_df_tc)


#experiment 8
bird02_exp08_path<-"./2014-08-15/vary freq-3pulse/"
bird02_exp08_list<-WLfileinfo(bird02_exp08_path)
bird02_exp08_files<-rownames(bird02_exp08_list)
exp08_bird02_alltrials<-NULL
for(i in 1:length(bird02_exp08_files)){
  exp08_bird02_alltrials[[i]]<-readAnalyzeWL(bird02_exp08_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird02_exp08_df<-WLtrialsSummary(bird02_exp08_list,exp08_bird02_alltrials)
bird02_exp08_df$phase<- -30
bird02_exp08_df$massspecWork=bird02_exp08_df$meanWork/bird02_pecmass
bird02_exp08_df$massspecPower=bird02_exp08_df$meanPower/bird02_pecmass
bird02_exp08_df_tc<-timeCorrect(bird02_exp08_df)
exp08_bird02_summary<-data.frame(experiment="exp08",bird="bird02",
                                 exp_bird="exp08_bird02",pec_mass=bird02_pecmass,
                                 bird02_exp08_df_tc)


#experiment 11
bird02_exp11_path<-"./2014-08-15/vary amp-4pulse/"
bird02_exp11_list<-WLfileinfo(bird02_exp11_path)
bird02_exp11_files<-rownames(bird02_exp11_list)
exp11_bird02_alltrials<-NULL
for(i in 1:length(bird02_exp11_files)){
  exp11_bird02_alltrials[[i]]<-readAnalyzeWL(bird02_exp11_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird02_exp11_df<-WLtrialsSummary(bird02_exp11_list,exp11_bird02_alltrials)
bird02_exp11_df$phase<- -30
bird02_exp11_df$massspecWork=bird02_exp11_df$meanWork/bird02_pecmass
bird02_exp11_df$massspecPower=bird02_exp11_df$meanPower/bird02_pecmass
bird02_exp11_df_tc<-timeCorrect(bird02_exp11_df)
exp11_bird02_summary<-data.frame(experiment="exp11",bird="bird02",
                                 exp_bird="exp11_bird02",pec_mass=bird02_pecmass,
                                 bird02_exp11_df_tc)

#experiment 12
bird02_exp12_path<-"./2014-08-15/vary freq -30% 4pulse do-over/"
bird02_exp12_list<-WLfileinfo(bird02_exp12_path)
bird02_exp12_files<-rownames(bird02_exp12_list)
exp12_bird02_alltrials<-NULL
for(i in 1:length(bird02_exp12_files)){
  exp12_bird02_alltrials[[i]]<-readAnalyzeWL(bird02_exp12_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird02_exp12_df<-WLtrialsSummary(bird02_exp12_list,exp12_bird02_alltrials)
bird02_exp12_df$phase<- -30
bird02_exp12_df$massspecWork=bird02_exp12_df$meanWork/bird02_pecmass
bird02_exp12_df$massspecPower=bird02_exp12_df$meanPower/bird02_pecmass
bird02_exp12_df_tc<-timeCorrect(bird02_exp12_df)
exp12_bird02_summary<-data.frame(experiment="exp12",bird="bird02",
                                 exp_bird="exp12_bird02",pec_mass=bird02_pecmass,
                                 bird02_exp12_df_tc)

#experiment 13
bird02_exp13_path<-"./2014-08-15/vary amp -20% 4 pulse/"
bird02_exp13_list<-WLfileinfo(bird02_exp13_path)
bird02_exp13_files<-rownames(bird02_exp13_list)
exp13_bird02_alltrials<-NULL
for(i in 1:length(bird02_exp13_files)){
  exp13_bird02_alltrials[[i]]<-readAnalyzeWL(bird02_exp13_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird02_exp13_df<-WLtrialsSummary(bird02_exp13_list,exp13_bird02_alltrials)
bird02_exp13_df$phase<- -30
bird02_exp13_df$massspecWork=bird02_exp13_df$meanWork/bird02_pecmass
bird02_exp13_df$massspecPower=bird02_exp13_df$meanPower/bird02_pecmass
bird02_exp13_df_tc<-timeCorrect(bird02_exp13_df)
exp13_bird02_summary<-data.frame(experiment="exp13",bird="bird02",
                                 exp_bird="exp13_bird02",pec_mass=bird02_pecmass,
                                 bird02_exp13_df_tc)



#now correct for degradation across experiments
exp04_bird02_summary$expCorWork<-exp04_bird02_summary$timeCorWork*bird02_baseline_degradation$corFac[6]
exp04_bird02_summary$expCorPower<-exp04_bird02_summary$timeCorPower*bird02_baseline_degradation$corFac[6]
exp05_bird02_summary$expCorWork<-exp05_bird02_summary$timeCorWork*bird02_baseline_degradation$corFac[4]
exp05_bird02_summary$expCorPower<-exp05_bird02_summary$timeCorPower*bird02_baseline_degradation$corFac[4]
exp06_bird02_summary$expCorWork<-exp06_bird02_summary$timeCorWork*bird02_baseline_degradation$corFac[5]
exp06_bird02_summary$expCorPower<-exp06_bird02_summary$timeCorPower*bird02_baseline_degradation$corFac[5]
exp07_bird02_summary$expCorWork<-exp07_bird02_summary$timeCorWork*bird02_baseline_degradation$corFac[3]
exp07_bird02_summary$expCorPower<-exp07_bird02_summary$timeCorPower*bird02_baseline_degradation$corFac[3]
exp08_bird02_summary$expCorWork<-exp08_bird02_summary$timeCorWork*bird02_baseline_degradation$corFac[1]
exp08_bird02_summary$expCorPower<-exp08_bird02_summary$timeCorPower*bird02_baseline_degradation$corFac[1]
exp11_bird02_summary$expCorWork<-exp11_bird02_summary$timeCorWork*bird02_baseline_degradation$corFac[2]
exp11_bird02_summary$expCorPower<-exp11_bird02_summary$timeCorPower*bird02_baseline_degradation$corFac[2]
exp12_bird02_summary$expCorWork<-exp12_bird02_summary$timeCorWork*bird02_baseline_degradation$corFac[8]
exp12_bird02_summary$expCorPower<-exp12_bird02_summary$timeCorPower*bird02_baseline_degradation$corFac[8]
exp13_bird02_summary$expCorWork<-exp13_bird02_summary$timeCorWork*bird02_baseline_degradation$corFac[7]
exp13_bird02_summary$expCorPower<-exp13_bird02_summary$timeCorPower*bird02_baseline_degradation$corFac[7]


#combine for a bird-specific summary object
bird02_summary<-rbind(exp04_bird02_summary,exp05_bird02_summary,exp06_bird02_summary,
                      exp07_bird02_summary,exp08_bird02_summary,
                      exp11_bird02_summary,exp12_bird02_summary,exp13_bird02_summary)
#plot(bird02_summary$pulses,bird02_summary$expCorPower,col=as.factor(bird02_summary$experiment),pch=19)
#  abline(h=0)