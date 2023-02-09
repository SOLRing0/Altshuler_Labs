####bird 06#######
#Bird 6: 2014-09-10

# twitch data
# tetanic not recorded
bird06_peakTetanic<-0

bird06_twitch<-read.ddf.isometric("./2014-09-10/twitch004.ddf")
twitchTiming(bird06_twitch)->bird06_twitch
bird06_TTPF<-attr(bird06_twitch,"Time_to_peak_force_secs") * 1000 # convert secs to MS


#involved in experiments 4-8 only
#import relevant baselines
#WL013 baseline before 26Hz (experiment 4)
bird06_baseWL_exp04<-readAnalyzeWL("./2014-09-10/WL006.ddf")
#WL011 baseline before 28Hz (experiment 5)
bird06_baseWL_exp05<-readAnalyzeWL("./2014-09-10/WL007.ddf")
#WL012 baseline before 30Hz (experiment 6)
bird06_baseWL_exp06<-readAnalyzeWL("./2014-09-10/WL008.ddf")
# assess degradation via net power output 
# ("meanPower" is mean of net power over 3 cycles)
bird06_baseWL_exp04_meanPower<-mean(ldply(bird06_baseWL_exp04$NetPower)$V1)
bird06_baseWL_exp05_meanPower<-mean(ldply(bird06_baseWL_exp05$NetPower)$V1)
bird06_baseWL_exp06_meanPower<-mean(ldply(bird06_baseWL_exp06$NetPower)$V1)
bird06_baseline_degradation<-data.frame(bird="bird06",
                                        exp=4:6,
                                        meanPower=c(
                                          mean(ldply(bird06_baseWL_exp04$NetPower)$V1),
                                          mean(ldply(bird06_baseWL_exp05$NetPower)$V1),
                                          mean(ldply(bird06_baseWL_exp06$NetPower)$V1))
)
bird06_baseline_degradation$percent<-c((100*(abs(bird06_baseWL_exp04_meanPower)/max(abs(bird06_baseline_degradation$meanPower)))),
                                       (100*(abs(bird06_baseWL_exp05_meanPower)/max(abs(bird06_baseline_degradation$meanPower)))),
                                       (100*(abs(bird06_baseWL_exp06_meanPower)/max(abs(bird06_baseline_degradation$meanPower))))
)
bird06_baseline_degradation$time<-c(file.info("./2014-09-10/WL006.ddf")$mtime,
                                    file.info("./2014-09-10/WL007.ddf")$mtime,
                                    file.info("./2014-09-10/WL008.ddf")$mtime)
bird06_baseline_degradation<-bird06_baseline_degradation[order(bird06_baseline_degradation$time),]
bird06_baseline_degradation$corFac<-(max(bird06_baseline_degradation$meanPower))/bird06_baseline_degradation$meanPower
plot(bird06_baseline_degradation$time,bird06_baseline_degradation$percent,pch=19,ylim=c(0,100))
lines(bird06_baseline_degradation$time,bird06_baseline_degradation$percent)
abline(h=80,lty=2,col='red')


#####import trials#####
#Bird 6: 2014-09-10
#import experiment data
#experiment 4
bird06_exp04_path<-"./2014-09-10/vary phase 3 pulse/"
bird06_exp04_list<-WLfileinfo(bird06_exp04_path)
bird06_exp04_files<-rownames(bird06_exp04_list)
exp04_bird06_alltrials<-NULL
for(i in 1:length(bird06_exp04_files)){
  exp04_bird06_alltrials[[i]]<-readAnalyzeWL(bird06_exp04_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird06_exp04_df<-WLtrialsSummary(bird06_exp04_list,exp04_bird06_alltrials)
bird06_exp04_df$phase<-c(-20,10,-10,-30,0,-20,-15,-25,-35,-20) 
bird06_exp04_df$massspecWork=bird06_exp04_df$meanWork/bird06_pecmass
bird06_exp04_df$massspecPower=bird06_exp04_df$meanPower/bird06_pecmass
bird06_exp04_df_tc<-timeCorrect(bird06_exp04_df)
exp04_bird06_summary<-data.frame(experiment="exp04",bird="bird06",
                                 exp_bird="exp04_bird06",pec_mass=bird06_pecmass,
                                 bird06_exp04_df_tc)


#experiment 5
bird06_exp05_path<-"./2014-09-10/vary phase 4 pulse/"
bird06_exp05_list<-WLfileinfo(bird06_exp05_path)
bird06_exp05_files<-rownames(bird06_exp05_list)
exp05_bird06_alltrials<-NULL
for(i in 1:length(bird06_exp05_files)){
  exp05_bird06_alltrials[[i]]<-readAnalyzeWL(bird06_exp05_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird06_exp05_df<-WLtrialsSummary(bird06_exp05_list,exp05_bird06_alltrials)
bird06_exp05_df$phase<-c(-20,10,-10,-30,-20,0,-15,-25,-35,-20) 
bird06_exp05_df$massspecWork=bird06_exp05_df$meanWork/bird06_pecmass
bird06_exp05_df$massspecPower=bird06_exp05_df$meanPower/bird06_pecmass
bird06_exp05_df_tc<-timeCorrect(bird06_exp05_df)
exp05_bird06_summary<-data.frame(experiment="exp05",bird="bird06",
                                 exp_bird="exp05_bird06",pec_mass=bird06_pecmass,
                                 bird06_exp05_df_tc)


#experiment 6
bird06_exp06_path<-"./2014-09-10/vary phase 5 pulse//"
bird06_exp06_list<-WLfileinfo(bird06_exp06_path)
bird06_exp06_files<-rownames(bird06_exp06_list)
exp06_bird06_alltrials<-NULL
for(i in 1:length(bird06_exp06_files)){
  exp06_bird06_alltrials[[i]]<-readAnalyzeWL(bird06_exp06_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird06_exp06_df<-WLtrialsSummary(bird06_exp06_list,exp06_bird06_alltrials)
bird06_exp06_df$phase<-c(-20,10,-10,-30,0,-20,-15,-25,-35,-20) 
bird06_exp06_df$massspecWork=bird06_exp06_df$meanWork/bird06_pecmass
bird06_exp06_df$massspecPower=bird06_exp06_df$meanPower/bird06_pecmass
bird06_exp06_df_tc<-timeCorrect(bird06_exp06_df)
exp06_bird06_summary<-data.frame(experiment="exp06",bird="bird06",
                                 exp_bird="exp06_bird06",pec_mass=bird06_pecmass,
                                 bird06_exp06_df_tc)


#now correct for degradation across experiments
exp04_bird06_summary$expCorWork<-exp04_bird06_summary$timeCorWork*bird06_baseline_degradation$corFac[1]
exp04_bird06_summary$expCorPower<-exp04_bird06_summary$timeCorPower*bird06_baseline_degradation$corFac[1]
exp05_bird06_summary$expCorWork<-exp05_bird06_summary$timeCorWork*bird06_baseline_degradation$corFac[2]
exp05_bird06_summary$expCorPower<-exp05_bird06_summary$timeCorPower*bird06_baseline_degradation$corFac[2]
exp06_bird06_summary$expCorWork<-exp06_bird06_summary$timeCorWork*bird06_baseline_degradation$corFac[3]
exp06_bird06_summary$expCorPower<-exp06_bird06_summary$timeCorPower*bird06_baseline_degradation$corFac[3]

#combine for a bird-specific summary object
bird06_summary<-rbind(exp04_bird06_summary,exp05_bird06_summary,exp06_bird06_summary)
#plot(bird06_summary$pulses,bird06_summary$expCorPower,col=as.factor(bird06_summary$experiment),pch=19)
#  abline(h=0)
