####bird 10#######
#Bird 10: 2015-03-13

# twitch data
bird10_peakTetanic<-0

bird10_twitch<-read.ddf.isometric("./2015-03-13/twitch001.ddf")
twitchTiming(bird10_twitch)->bird10_twitch
bird10_TTPF<-attr(bird10_twitch,"Time_to_peak_force_secs") * 1000 # convert secs to MS


#involved in experiments 10,12,14 only
#WL006 baseline before 26Hz (experiment 1)
bird10_baseWL_exp10<-readAnalyzeWL("./2015-03-13/WL008.ddf")
#WL006 baseline before 26Hz (experiment 1)
bird10_baseWL_exp12<-readAnalyzeWL("./2015-03-13/WL009.ddf")
#WL007 baseline before 28Hz (experiment 2)
bird10_baseWL_exp14<-readAnalyzeWL("./2015-03-13/WL010.ddf")
# assess degradation via net power output 
# ("meanPower" is mean of net power over 3 cycles)
bird10_baseWL_exp10_meanPower<-mean(ldply(bird10_baseWL_exp10$NetPower)$V1)
bird10_baseWL_exp12_meanPower<-mean(ldply(bird10_baseWL_exp12$NetPower)$V1)
bird10_baseWL_exp14_meanPower<-mean(ldply(bird10_baseWL_exp14$NetPower)$V1)
bird10_baseline_degradation<-data.frame(bird="bird10",
                                        exp=c(10,12,14),
                                        meanPower=c(
                                          mean(ldply(bird10_baseWL_exp10$NetPower)$V1),
                                          mean(ldply(bird10_baseWL_exp12$NetPower)$V1),
                                          mean(ldply(bird10_baseWL_exp14$NetPower)$V1))
)
bird10_baseline_degradation$percent<-c((100*(abs(bird10_baseWL_exp10_meanPower)/max(abs(bird10_baseline_degradation$meanPower)))),
                                       (100*(abs(bird10_baseWL_exp12_meanPower)/max(abs(bird10_baseline_degradation$meanPower)))),
                                       (100*(abs(bird10_baseWL_exp14_meanPower)/max(abs(bird10_baseline_degradation$meanPower))))
)
bird10_baseline_degradation$time<-c(file.info("./2015-03-13/WL008.ddf")$mtime,
                                    file.info("./2015-03-13/WL009.ddf")$mtime,
                                    file.info("./2015-03-13/WL010.ddf")$mtime)
bird10_baseline_degradation<-bird10_baseline_degradation[order(bird10_baseline_degradation$time),]
bird10_baseline_degradation$corFac<-(max(bird10_baseline_degradation$meanPower))/bird10_baseline_degradation$meanPower
plot(bird10_baseline_degradation$time,bird10_baseline_degradation$percent,pch=19,ylim=c(0,100))
lines(bird10_baseline_degradation$time,bird10_baseline_degradation$percent)
abline(h=80,lty=2,col='red')


#####import trials#####
#Bird 10: 2015-03-13
#import experiment data
#experiment 10
bird10_exp10_path<-"./2015-03-13/vary freq -20% 3 pulse CP/"
bird10_exp10_list<-WLfileinfo(bird10_exp10_path)
bird10_exp10_files<-rownames(bird10_exp10_list)
exp10_bird10_alltrials<-NULL
for(i in 1:length(bird10_exp10_files)){
  exp10_bird10_alltrials[[i]]<-readAnalyzeWL(bird10_exp10_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird10_exp10_df<-WLtrialsSummary(bird10_exp10_list,exp10_bird10_alltrials)
bird10_exp10_df$phase<- -20
bird10_exp10_df$massspecWork=bird10_exp10_df$meanWork/bird10_pecmass
bird10_exp10_df$massspecPower=bird10_exp10_df$meanPower/bird10_pecmass
bird10_exp10_df_tc<-timeCorrect(bird10_exp10_df)
exp10_bird10_summary<-data.frame(experiment="exp10",bird="bird10",
                                 exp_bird="exp10_bird10",pec_mass=bird10_pecmass,
                                 bird10_exp10_df_tc)


#experiment 12
bird10_exp12_path<-"./2015-03-13/vary freq -30% 4 pulse CP/"
bird10_exp12_list<-WLfileinfo(bird10_exp12_path)
bird10_exp12_files<-rownames(bird10_exp12_list)
exp12_bird10_alltrials<-NULL
for(i in 1:length(bird10_exp12_files)){
  exp12_bird10_alltrials[[i]]<-readAnalyzeWL(bird10_exp12_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird10_exp12_df<-WLtrialsSummary(bird10_exp12_list,exp12_bird10_alltrials)
bird10_exp12_df$phase<- -30
bird10_exp12_df$massspecWork=bird10_exp12_df$meanWork/bird10_pecmass
bird10_exp12_df$massspecPower=bird10_exp12_df$meanPower/bird10_pecmass
bird10_exp12_df_tc<-timeCorrect(bird10_exp12_df)
exp12_bird10_summary<-data.frame(experiment="exp12",bird="bird10",
                                 exp_bird="exp12_bird10",pec_mass=bird10_pecmass,
                                 bird10_exp12_df_tc)


#experiment 14
bird10_exp14_path<-"./2015-03-13/vary freq -20% 4 pulse CP/"
bird10_exp14_list<-WLfileinfo(bird10_exp14_path)
bird10_exp14_files<-rownames(bird10_exp14_list)
exp14_bird10_alltrials<-NULL
for(i in 1:length(bird10_exp14_files)){
  exp14_bird10_alltrials[[i]]<-readAnalyzeWL(bird10_exp14_files[[i]],
                                             bworth_freq=0.05,
                                             keep.cycles=3:5,GR=2)
}
#now analyze the data and produce a summary
bird10_exp14_df<-WLtrialsSummary(bird10_exp14_list,exp14_bird10_alltrials)
bird10_exp14_df$phase<- -20
bird10_exp14_df$massspecWork=bird10_exp14_df$meanWork/bird10_pecmass
bird10_exp14_df$massspecPower=bird10_exp14_df$meanPower/bird10_pecmass
#bird10_exp14_df_tc<-timeCorrect(bird10_exp14_df)
#can't rely on usual timeCorrection function as last trial is not same as first
#do it manually instead.
ip<-bird10_exp14_df$massspecPower[1]
fp<-bird10_exp14_df$massspecPower[8]
pd<-ip-fp
td<-bird10_exp14_df$time[8]-bird10_exp14_df$time[1]
oc<-pd/as.numeric(td)
tmp<-NULL
for (i in 1:dim(bird10_exp14_df)[1]){
  tmp[i]<-bird10_exp14_df$time[i]-bird10_exp14_df$time[1]
}
bird10_exp14_df->bird10_exp14_df_tc
bird10_exp14_df_tc$elapsed<-tmp
bird10_exp14_df_tc$timeCorWork<-(bird10_exp14_df$massspecPower + (tmp*oc))/bird10_exp14_df$freq
bird10_exp14_df_tc$timeCorPower<-bird10_exp14_df$massspecPower + (tmp*oc)
exp14_bird10_summary<-data.frame(experiment="exp14",bird="bird10",
                                 exp_bird="exp14_bird10",pec_mass=bird10_pecmass,
                                 bird10_exp14_df_tc)


#now correct for degradation across experiments
exp10_bird10_summary$expCorWork<-exp10_bird10_summary$timeCorWork*bird10_baseline_degradation$corFac[1]
exp10_bird10_summary$expCorPower<-exp10_bird10_summary$timeCorPower*bird10_baseline_degradation$corFac[1]
exp12_bird10_summary$expCorWork<-exp12_bird10_summary$timeCorWork*bird10_baseline_degradation$corFac[2]
exp12_bird10_summary$expCorPower<-exp12_bird10_summary$timeCorPower*bird10_baseline_degradation$corFac[2]
exp14_bird10_summary$expCorWork<-exp14_bird10_summary$timeCorWork*bird10_baseline_degradation$corFac[3]
exp14_bird10_summary$expCorPower<-exp14_bird10_summary$timeCorPower*bird10_baseline_degradation$corFac[3]


#combine for a bird-specific summary object
bird10_summary<-rbind(exp10_bird10_summary,exp12_bird10_summary,exp14_bird10_summary)
#plot(bird10_summary$pulses,bird10_summary$expCorPower,col=as.factor(bird10_summary$experiment),pch=19)
#  abline(h=0)