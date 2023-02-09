# please be sure to run 01_allFunctions and 02_dataImport prior to running
# anything in this script

################################## filter data #################################
# filter summary data to experiments 7-14 only
ampfreqdata_exp07<-allsummarydat[allsummarydat$experiment=="exp07",]
ampfreqdata_exp08<-allsummarydat[allsummarydat$experiment=="exp08",]
ampfreqdata_exp09<-allsummarydat[allsummarydat$experiment=="exp09",]
ampfreqdata_exp10<-allsummarydat[allsummarydat$experiment=="exp10",]
ampfreqdata_exp11<-allsummarydat[allsummarydat$experiment=="exp11",]
ampfreqdata_exp12<-allsummarydat[allsummarydat$experiment=="exp12",]
ampfreqdata_exp13<-allsummarydat[allsummarydat$experiment=="exp13",]
ampfreqdata_exp14<-allsummarydat[allsummarydat$experiment=="exp14",]
ampfreqdata_all<-rbind(ampfreqdata_exp07,ampfreqdata_exp08,ampfreqdata_exp08,
                       ampfreqdata_exp10,ampfreqdata_exp11,ampfreqdata_exp12,
                       ampfreqdata_exp13,ampfreqdata_exp14)


###################### mixed-model fitting and comparison ######################

# use a comparative mixed-models approach to determine how changes to 
# amplitude or frequency affect net power output

# first, we should standardize the data because their units vary widely 
# in magnitude
tmp<-scale(data.frame(expCorPower=ampfreqdata_all$expCorPower,
                      freq=ampfreqdata_all$freq,
                      amplitude=ampfreqdata_all$amplitude,
                      pulses=ampfreqdata_all$pulses,
                      phase=ampfreqdata_all$phase))
# combine with bird IDs, which will be a random effect
standdat<-data.frame(bird=ampfreqdata_all$bird,
                     tmp)
# visualize
library(lattice)
splom(standdat,pch=19)


# Bayesian mixed-models via MCMCglmm
library(MCMCglmm);library(MuMIn)
fit_allinteract<-MCMCglmm(expCorPower~freq+amplitude+pulses+phase+freq:phase+
                            freq:pulses+amplitude:phase+amplitude:pulses-1,
                          random=~bird,data=standdat,
                          nitt=1000000,thin=100,burnin=100000,verbose=FALSE,
                          pr=TRUE,pl=TRUE)
fit_freqinteract<-MCMCglmm(expCorPower~freq+amplitude+pulses+phase+freq:phase+
                             freq:pulses-1,
                           random=~bird,data=standdat,
                           nitt=1000000,thin=100,burnin=100000,verbose=FALSE,
                           pr=TRUE,pl=TRUE)
fit_ampinteract<-MCMCglmm(expCorPower~freq+amplitude+pulses+phase+
                            amplitude:phase+amplitude:pulses-1,
                           random=~bird,data=standdat,
                           nitt=1000000,thin=100,burnin=100000,verbose=FALSE,
                           pr=TRUE,pl=TRUE)
fit_nointeract<-MCMCglmm(expCorPower~freq+amplitude+pulses+phase-1,
                         random=~bird,data=standdat,
                         nitt=1000000,thin=100,burnin=100000,verbose=FALSE,
                         pr=TRUE,pl=TRUE)
fit_noamp<-MCMCglmm(expCorPower~freq+pulses+phase-1,
                         random=~bird,data=standdat,
                         nitt=1000000,thin=100,burnin=100000,verbose=FALSE,
                         pr=TRUE,pl=TRUE)
fit_null<-MCMCglmm(expCorPower~1,
                   random=~bird,data=standdat,
                   nitt=1000000,thin=100,burnin=100000,verbose=FALSE,
                   pr=TRUE,pl=TRUE)

# glance at DICs
DIC(fit_allinteract);DIC(fit_freqinteract);DIC(fit_ampinteract);
DIC(fit_nointeract);DIC(fit_noamp);DIC(fit_null);

# best model
summary(fit_freqinteract)

# how strong are the random effects (birds)?
fit_freqinteract_rands<-(fit_freqinteract$VCV[,1])/(fit_freqinteract$VCV[,1]+
                                                      fit_freqinteract$VCV[,2])
plot(fit_freqinteract_rands);mean(fit_freqinteract_rands)

# assess predictions
modslist<-lapply(c("fit_allinteract","fit_freqinteract","fit_ampinteract",
                   "fit_nointeract","fit_noamp","fit_null"),get)
modterms<-c(8,6,6,4,3,2)
mods_R2s<-NULL
DICs<-NULL
effects<-NULL
for (i in 1:length(modslist)){
  preds_tmp<-predict.MCMCglmm(modslist[[i]],marginal=NULL)
  mods_R2s[i]<-summary(lm(standdat$expCorPower~preds_tmp))$adj.r.squared
  DICs[i]<-DIC(modslist[[i]])
  effects[[i]]<-colMeans(modslist[[i]]$Sol[,1:modterms[i]])
}
modeffects<-data.frame()
for(i in seq(along=effects)) for(j in names(effects[[i]]))
  modeffects[i,j]<-round(effects[[i]][j],3)
modtable<-data.frame(mods=c("fit_allinteract","fit_freqinteract",
                            "fit_ampinteract","fit_nointeract",
                            "fit_noamp","fit_null"),
           modeffects[,1:9],
           DIC=round(DICs,2),
           deltaDIC=round(DICs-min(DICs),2),
           R2=round(mods_R2s,3))
modtable<-modtable[order(modtable$deltaDIC),]

################################# model checks #################################

# check the model via the performance package
library(performance)

# as of Fall 2019, MCMCglmm objects are not supported by the performance package
# we will make a similar model via lme4
library(lme4)
model <- lme4::lmer(expCorPower~freq+amplitude+pulses+phase+freq:phase+
                freq:pulses+(1|bird)-1, data=standdat)
# summary(model) shows that the coefficients are nearly identical to those 
# we found in fit_freqinteract
r2(model) # high r2
check_overdispersion(model) # no overdispersion
check_singularity(model) # no singularity
model_performance(model) # metrics match what we see in modtable

# power analysis
library(simr)
# use power simluation to determine both the power as well as the effect size
# for each predictor
powerSim(model,test=fixed("amplitude"),nsim=100)
powerSim(model,test=fixed("freq"),nsim=100)
powerSim(model,test=fixed("pulses"),nsim=100)
powerSim(model,test=fixed("phase"),nsim=100)
powerSim(model,test=fixed("freq:pulses"),nsim=100)
powerSim(model,test=fixed("freq:phase"),nsim=100)

# we can also use likelihood ratio tests against simpler models, e.g.
doTest(model,fcompare(~freq+amplitude+pulses+phase))
doTest(model,fcompare(~freq+pulses+phase))
doTest(model,fcompare(~freq+amplitude))
doTest(model,fcompare(~freq))
# our model favors comparably well against each of these alternates



