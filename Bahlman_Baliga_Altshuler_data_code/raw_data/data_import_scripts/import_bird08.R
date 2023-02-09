####bird 08#######
#Bird 8: 2014-10-03

# twitch data
bird08_tetanus<-read.ddf.isometric("./2014-10-03/tetnus 01.ddf")
bird08_peakTetanic<-max(bird08_tetanus$Force)*2/1000 # mult by 2 because of gear ratio then divided by 1000 to get from mN to N

bird08_twitch<-read.ddf.isometric("./2014-10-03/twitch003.ddf")
twitchTiming(bird08_twitch)->bird08_twitch
bird08_TTPF<-attr(bird08_twitch,"Time_to_peak_force_secs") * 1000 # convert secs to MS

#involved in experiments 7-9,11-13 only
#import relevant baselines
#WL013 baseline before 26Hz (experiment 7)
bird08_baseWL_exp07<-readAnalyzeWL("./2014-10-03/WL015.ddf")
#WL011 baseline before 28Hz (experiment 9)
bird08_baseWL_exp08<-readAnalyzeWL("./2014-10-03/WL014.ddf")
#WL012 baseline before 30Hz (experiment 11)
bird08_baseWL_exp09<-readAnalyzeWL("./2014-10-03/WL016.ddf")
#WL010 baseline before 3pulse (experiment 12)
bird08_baseWL_exp11<-readAnalyzeWL("./2014-10-03/WL017.ddf")
#WL008 baseline before 4pulse (experiment 13)
bird08_baseWL_exp13<-readAnalyzeWL("./2014-10-03/WL018.ddf")
# assess degradation via net power output 
# ("meanPower" is mean of net power over 3 cycles)
bird08_baseWL_exp07_meanPower<-mean(ldply(bird08_baseWL_exp07$NetPower)$V1)
bird08_baseWL_exp08_meanPower<-mean(ldply(bird08_baseWL_exp08$NetPower)$V1)
bird08_baseWL_exp09_meanPower<-mean(ldply(bird08_baseWL_exp09$NetPower)$V1)
bird08_baseWL_exp11_meanPower<-mean(ldply(bird08_baseWL_exp11$NetPower)$V1)
bird08_baseWL_exp13_meanPower<-mean(ldply(bird08_baseWL_exp13$NetPower)$V1)
bird08_baseline_degradation<-data.frame(bird="bird08",
                                        exp=c(7:9,11,13),
                                        meanPower=c(
                                          mean(ldply(bird08_baseWL_exp07$NetPower)$V1),
                                          mean(ldply(bird08_baseWL_exp08$NetPower)$V1),
                                          mean(ldply(bird08_baseWL_exp09$NetPower)$V1),
                                          mean(ldply(bird08_baseWL_exp11$NetPower)$V1),
                                          mean(ldply(bird08_baseWL_exp13$NetPower)$V1))
)
bird08_baseline_degradation$percent<-c((100*(abs(bird08_baseWL_exp07_meanPower)/max(abs(bird08_baseline_degradation$meanPower)))),
                                       (100*(abs(bird08_baseWL_exp08_meanPower)/max(abs(bird08_baseline_degradation$meanPower)))),
                                       (100*(abs(bird08_baseWL_exp09_meanPower)/max(abs(bird08_baseline_degradation$meanPower)))),
                                       (100*(abs(bird08_baseWL_exp11_meanPower)/max(abs(bird08_baseline_degradation$meanPower)))),
                                       (100*(abs(bird08_baseWL_exp13_meanPower)/max(abs(bird08_baseline_degradation$meanPower))))
)
bird08_baseline_degradation$time<-c(file.info("./2014-10-03/WL015.ddf")$mtime,
                                    file.info("./2014-10-03/WL014.ddf")$mtime,
                                    file.info("./2014-10-03/WL016.ddf")$mtime,
                                    file.info("./2014-10-03/WL017.ddf")$mtime,
                                    file.info("./2014-10-03/WL018.ddf")$mtime)
bird08_baseline_degradation<-bird08_baseline_degradation[order(bird08_baseline_degradation$time),]
bird08_baseline_degradation$corFac<-0.937*(max(bird08_baseline_degradation$meanPower))/bird08_baseline_degradation$meanPower
plot(bird08_baseline_degradation$time,bird08_baseline_degradation$percent,pch=19,ylim=c(0,100))
lines(bird08_baseline_degradation$time,bird08_baseline_degradation$percent)
abline(h=80,lty=2,col='red')


#####import trials#####
#Bird 8: 2014-10-03
#import experiment data
#experiment 7
bird08_exp07_path<-"./2014-10-03/vary amp -30% 3 pulse/"
bird08_exp07_list<-WLfileinfo(bird08_exp07_path)
bird08_exp07_files<-rownames(bird08_exp07_list)
exp07_bird08_alltrials<-NULL
for(i in 1:length(bird08_exp07_files)){
  exp07_bird08_alltrials[[i]]<-readAnalyzeWL(bird08_exp07_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird08_exp07_df<-WLtrialsSummary(bird08_exp07_list,exp07_bird08_alltrials)
bird08_exp07_df$phase<- -30
bird08_exp07_df$massspecWork=bird08_exp07_df$meanWork/bird08_pecmass
bird08_exp07_df$massspecPower=bird08_exp07_df$meanPower/bird08_pecmass
bird08_exp07_df_tc<-timeCorrect(bird08_exp07_df)
exp07_bird08_summary<-data.frame(experiment="exp07",bird="bird08",
                                 exp_bird="exp07_bird08",pec_mass=bird08_pecmass,
                                 bird08_exp07_df_tc)


#experiment 8
bird08_exp08_path<-"./2014-10-03/vary freq -30% 3 pulse/"
bird08_exp08_list<-WLfileinfo(bird08_exp08_path)
bird08_exp08_files<-rownames(bird08_exp08_list)
exp08_bird08_alltrials<-NULL
for(i in 1:length(bird08_exp08_files)){
  exp08_bird08_alltrials[[i]]<-readAnalyzeWL(bird08_exp08_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird08_exp08_df<-WLtrialsSummary(bird08_exp08_list,exp08_bird08_alltrials)
bird08_exp08_df$phase<- -30
bird08_exp08_df$massspecWork=bird08_exp08_df$meanWork/bird08_pecmass
bird08_exp08_df$massspecPower=bird08_exp08_df$meanPower/bird08_pecmass
#can't rely on usual timeCorrection function as first trial is not same as final
#bird08_exp08_df_tc<-timeCorrect(bird08_exp08_df)
#do it manually instead. This will correct trials 4-8 only as we don't have a 
#good reference for the first three trials.
ip<-bird08_exp08_df$massspecPower[4]
fp<-bird08_exp08_df$massspecPower[8]
pd<-ip-fp
td<-bird08_exp08_df$time[8]-bird08_exp08_df$time[4]
oc<-pd/as.numeric(td)
tmp<-NULL
for (i in 1:dim(bird08_exp08_df)[1]){
  tmp[i]<-bird08_exp08_df$time[i]-bird08_exp08_df$time[1]
}
bird08_exp08_df->bird08_exp08_df_tc
bird08_exp08_df_tc$elapsed<-tmp
tmp<-tmp - tmp[4]
bird08_exp08_df_tc$timeCorWork<-(bird08_exp08_df$massspecPower + (tmp*oc))/bird08_exp08_df$freq
bird08_exp08_df_tc$timeCorPower<-bird08_exp08_df$massspecPower + (tmp*oc)
bird08_exp08_df_tc$timeCorWork[1:3]<-bird08_exp08_df_tc$massspecWork[1:3]
bird08_exp08_df_tc$timeCorPower[1:3]<-bird08_exp08_df_tc$massspecPower[1:3]
exp08_bird08_summary<-data.frame(experiment="exp08",bird="bird08",
                                 exp_bird="exp08_bird08",pec_mass=bird08_pecmass,
                                 bird08_exp08_df_tc)


#experiment 9
bird08_exp09_path<-"./2014-10-03/vary amp -20% 3 pulse/"
bird08_exp09_list<-WLfileinfo(bird08_exp09_path)
bird08_exp09_files<-rownames(bird08_exp09_list)
exp09_bird08_alltrials<-NULL
for(i in 1:length(bird08_exp09_files)){
  exp09_bird08_alltrials[[i]]<-readAnalyzeWL(bird08_exp09_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird08_exp09_df<-WLtrialsSummary(bird08_exp09_list,exp09_bird08_alltrials)
bird08_exp09_df$phase<- -20
bird08_exp09_df$massspecWork=bird08_exp09_df$meanWork/bird08_pecmass
bird08_exp09_df$massspecPower=bird08_exp09_df$meanPower/bird08_pecmass
bird08_exp09_df_tc<-timeCorrect(bird08_exp09_df)
exp09_bird08_summary<-data.frame(experiment="exp09",bird="bird08",
                                 exp_bird="exp09_bird08",pec_mass=bird08_pecmass,
                                 bird08_exp09_df_tc)


#experiment 11
bird08_exp11_path<-"./2014-10-03/vary amp -30% 4 pulse/"
bird08_exp11_list<-WLfileinfo(bird08_exp11_path)
bird08_exp11_files<-rownames(bird08_exp11_list)
exp11_bird08_alltrials<-NULL
for(i in 1:length(bird08_exp11_files)){
  exp11_bird08_alltrials[[i]]<-readAnalyzeWL(bird08_exp11_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird08_exp11_df<-WLtrialsSummary(bird08_exp11_list,exp11_bird08_alltrials)
bird08_exp11_df$phase<- -30
bird08_exp11_df$massspecWork=bird08_exp11_df$meanWork/bird08_pecmass
bird08_exp11_df$massspecPower=bird08_exp11_df$meanPower/bird08_pecmass
bird08_exp11_df_tc<-timeCorrect(bird08_exp11_df)
exp11_bird08_summary<-data.frame(experiment="exp11",bird="bird08",
                                 exp_bird="exp11_bird08",pec_mass=bird08_pecmass,
                                 bird08_exp11_df_tc)


#experiment 13
bird08_exp13_path<-"./2014-10-03/vary amp -20% 4 pulse/"
bird08_exp13_list<-WLfileinfo(bird08_exp13_path)
bird08_exp13_files<-rownames(bird08_exp13_list)
exp13_bird08_alltrials<-NULL
for(i in 1:length(bird08_exp13_files)){
  exp13_bird08_alltrials[[i]]<-readAnalyzeWL(bird08_exp13_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird08_exp13_df<-WLtrialsSummary(bird08_exp13_list,exp13_bird08_alltrials)
bird08_exp13_df$phase<- -20
bird08_exp13_df$massspecWork=bird08_exp13_df$meanWork/bird08_pecmass
bird08_exp13_df$massspecPower=bird08_exp13_df$meanPower/bird08_pecmass
bird08_exp13_df_tc<-timeCorrect(bird08_exp13_df)
exp13_bird08_summary<-data.frame(experiment="exp13",bird="bird08",
                                 exp_bird="exp13_bird08",pec_mass=bird08_pecmass,
                                 bird08_exp13_df_tc)


#now correct for degradation across experiments
exp07_bird08_summary$expCorWork<-exp07_bird08_summary$timeCorWork*bird08_baseline_degradation$corFac[2]
exp07_bird08_summary$expCorPower<-exp07_bird08_summary$timeCorPower*bird08_baseline_degradation$corFac[2]
exp08_bird08_summary$expCorWork<-exp08_bird08_summary$timeCorWork*bird08_baseline_degradation$corFac[1]
exp08_bird08_summary$expCorPower<-exp08_bird08_summary$timeCorPower*bird08_baseline_degradation$corFac[1]
exp09_bird08_summary$expCorWork<-exp09_bird08_summary$timeCorWork*bird08_baseline_degradation$corFac[3]
exp09_bird08_summary$expCorPower<-exp09_bird08_summary$timeCorPower*bird08_baseline_degradation$corFac[3]
exp11_bird08_summary$expCorWork<-exp11_bird08_summary$timeCorWork*bird08_baseline_degradation$corFac[4]
exp11_bird08_summary$expCorPower<-exp11_bird08_summary$timeCorPower*bird08_baseline_degradation$corFac[4]
exp13_bird08_summary$expCorWork<-exp13_bird08_summary$timeCorWork*bird08_baseline_degradation$corFac[5]
exp13_bird08_summary$expCorPower<-exp13_bird08_summary$timeCorPower*bird08_baseline_degradation$corFac[5]


#combine for a bird-specific summary object
bird08_summary<-rbind(exp07_bird08_summary,exp08_bird08_summary,exp09_bird08_summary,
                      exp11_bird08_summary,exp13_bird08_summary)
#plot(bird08_summary$pulses,bird08_summary$expCorPower,col=as.factor(bird08_summary$experiment),pch=19)
#  abline(h=0)