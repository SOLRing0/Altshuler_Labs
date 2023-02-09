#using a linear mixed effects model to account for repeated measures on the same individual 

library(nlme)
library(lme4)
library(lattice)
library(influence.ME)
library(multcompView)
library(visreg)

### ----- FIT MODEL TO AERODYNAMIC EFFICIENCY
#test to see if there is an effect due to elbow angle
test_full     <- lmer(MaxL_D ~ poly(Angle.True,3) + TI + (1|WingID), data = max_results, REML = F)
test_full_lme <- lme(MaxL_D ~ poly(Angle.True,3) + TI, random = ~1|WingID, data = max_results)
#reduced models
test_full_elbowred <- lmer(MaxL_D ~ 1 + TI + (1|WingID), data = max_results, REML = F)
test_full_TIred    <- lmer(MaxL_D ~ 1 + poly(Angle.True,3) + (1|WingID), data = max_results, REML = F)

anova(test_full,test_full_elbowred, test = "Chisq") #check for effect of elbow
anova(test_full,test_full_TIred, test = "Chisq")    #check for effect of TI

#get new predicted values
bb_ld     <- bootMer(test_full, FUN = function(x)predict(x,re.form=NA),nsim=100)
median_ld <-apply(bb_ld$t,2,quantile,0.5)
lower_ld  <-apply(bb_ld$t,2,quantile,0.025)
upper_ld  <-apply(bb_ld$t,2,quantile,0.975)

### ------ Fit model for CLMAX - quadratic
test_full_clmax <- lmer(CLmax ~ poly(Angle.True,1) + TI + (1|WingID), data = max_results, REML = F)
test_full_lme_clmax <- lme(CLmax ~ poly(Angle.True,1) + TI, random = ~1|WingID, data = max_results)
#reduced models
test_full_elbowred_clmax <- lmer(CLmax ~ 1 + TI + (1|WingID), data = max_results, REML = F)
test_full_TIred_clmax    <- lmer(CLmax ~ 1 + poly(Angle.True,1) + (1|WingID), data = max_results, REML = F)

anova(test_full_clmax,test_full_elbowred_clmax, test = "Chisq") #check for effect of elbow
anova(test_full_clmax,test_full_TIred_clmax, test = "Chisq")    #check for effect of TI

#get new predicted values
bb_clmax <- bootMer(test_full_clmax, FUN = function(x)predict(x,re.form=NA),nsim=100)
median_clmax <-apply(bb_clmax$t,2,quantile,0.5)
lower_clmax <-apply(bb_clmax$t,2,quantile,0.025)
upper_clmax <-apply(bb_clmax$t,2,quantile,0.975)

### ------ Fit model for CDMIN 
test_full_cdmin     <- lmer(CDmin ~ poly(Angle.True,2) + TI +(1|WingID), data = max_results, REML = F)
test_full_lme_cdmin <- lme(CDmin ~ poly(Angle.True,2) +TI, random = ~1|WingID, data = max_results)
#reduced models
test_full_elbowred_cdmin <- lmer(CDmin ~ 1 + TI+(1|WingID), data = max_results, REML = F)
test_full_TIred_cdmin    <- lmer(CDmin ~ 1 + poly(Angle.True,2) + (1|WingID), data = max_results, REML = F)

anova(test_full_cdmin,test_full_elbowred_cdmin, test = "Chisq") #check for effect of elbow
anova(test_full_cdmin,test_full_TIred_cdmin, test = "Chisq")    #check for effect of TI

#get new predicted values
bb_cdmin <- bootMer(test_full_cdmin, FUN = function(x)predict(x,re.form=NA),nsim=100)
median_cdmin <-apply(bb_cdmin$t,2,quantile,0.5)
lower_cdmin <-apply(bb_cdmin$t,2,quantile,0.025)
upper_cdmin <-apply(bb_cdmin$t,2,quantile,0.975)

### ------ Fit model for PITCHING SLOPE - quartic
test_full_cmcl     <- lmer(slopemean ~ poly(Angle.True,3) + TI + (1|WingID), data = fit_results, REML = F)
test_full_cmcl_lme <- lme(slopemean ~ poly(Angle.True,3) + TI, random = ~1|WingID, data = fit_results)
#reduced models
test_full_elbowred_cmcl <- lmer(slopemean ~ 1 +  TI + (1|WingID), data = fit_results, REML = F)
test_full_TIred_cmcl    <- lmer(slopemean ~ 1 + poly(Angle.True,3) + (1|WingID), data = fit_results, REML = F)

anova(test_full_cmcl,test_full_elbowred_cmcl, test = "Chisq") #check for effect of elbow
anova(test_full_cmcl,test_full_TIred_cmcl, test = "Chisq")    #check for effect of TI

#get new predicted values
bb_cmcl<- bootMer(test_full_cmcl, FUN = function(x)predict(x,re.form=NA),nsim=100)
median_cmcl <-apply(bb_cmcl$t,2,quantile,0.5)
lower_cmcl <-apply(bb_cmcl$t,2,quantile,0.025)
upper_cmcl <-apply(bb_cmcl$t,2,quantile,0.975)

### ------ Fit model for pitching intercept - cubic
#test to see if there is an effect due to elbow angle
test_full_cm0 <- lmer(Cm_AC ~ poly(Angle.True,1) + TI + (1|WingID), data = fit_results, REML = F)
test_full_cm0_lme <- lme(Cm_AC ~ poly(Angle.True,1) + TI, random = ~1|WingID, data = fit_results)
#reduced models
test_full_elbowred_cm0 <- lmer(Cm_AC ~ 1 + TI + (1|WingID), data = fit_results, REML = F)
test_full_TIred_cm0<- lmer(Cm_AC ~ poly(Angle.True,1) + 1 + (1|WingID), data = fit_results, REML = F)

anova(test_full_cm0,test_full_elbowred_cm0, test = "Chisq") #check for effect of elbow
anova(test_full_cm0,test_full_TIred_cm0, test = "Chisq")    #check for effect of TI

#get new predicted values
bb_cm0<- bootMer(test_full_cm0, FUN = function(x)predict(x,re.form=NA),nsim=100)
median_cm0 <-apply(bb_cm0$t,2,quantile,0.5)
lower_cm0 <-apply(bb_cm0$t,2,quantile,0.025)
upper_cm0 <-apply(bb_cm0$t,2,quantile,0.975)
