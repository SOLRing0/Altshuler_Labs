
x <- data.frame(matrix(nrow = 5, ncol = 12))

 test_full_clmax3 <- lmer(CLmax ~ poly(Angle.True,3) + TI + (1|WingID), data = max_results, REML = F)
 test_full_clmax1 <- lmer(CLmax ~ poly(Angle.True,1) + TI + (1|WingID), data = max_results, REML = F)
 test_full_clmax <- lmer(CLmax ~ poly(Angle.True,2) + TI + (1|WingID), data = max_results, REML = F)
x[1,] <- c(AIC(test_full_clmax), AIC(test_full_clmax1), AIC(test_full_clmax3)) # stores AIC values in a vector
### ---- BEST RESULT FOR CLmax -> LINEAR MODEL
##-----------------------------------------------
test_full1 <- lmer(MaxL_D ~ poly(Angle.True,1) + TI + (1|WingID), data = max_results, REML = F)
test_full <- lmer(MaxL_D ~ poly(Angle.True,2) + TI + (1|WingID), data = max_results, REML = F)
test_full3 <- lmer(MaxL_D ~ poly(Angle.True,3) + TI + (1|WingID), data = max_results, REML = F)
x[2,] <- c(AIC(test_full), AIC(test_full1), AIC(test_full3)) # stores AIC values in a vector
### ---- BEST RESULT FOR L/D -> CUBIC MODEL
##-----------------------------------------------
test_full_cdmin1 <- lmer(CDmin ~ poly(Angle.True,1) + TI +(1|WingID), data = max_results, REML = F)
test_full_cdmin <- lmer(CDmin ~ poly(Angle.True,2) + TI +(1|WingID), data = max_results, REML = F)
test_full_cdmin3 <- lmer(CDmin ~ poly(Angle.True,3) + TI +(1|WingID), data = max_results, REML = F)
x[3,] <- c(AIC(test_full_cdmin1), AIC(test_full_cdmin), AIC(test_full_cdmin3)) # stores AIC values in a vector
### ---- BEST RESULT FOR CDmin -> QUADRATIC MODEL
##-----------------------------------------------
test_full_cmcl1 <- lmer(slopemean ~ poly(Angle.True,1) + TI + (1|WingID), data = fit_results, REML = F)
test_full_cmcl2 <- lmer(slopemean ~ poly(Angle.True,2) + TI + (1|WingID), data = fit_results, REML = F)
test_full_cmcl <- lmer(slopemean ~ poly(Angle.True,3) + TI + (1|WingID), data = fit_results, REML = F)
x[4,] <- c(AIC(test_full_cmcl1), AIC(test_full_cmcl2), AIC(test_full_cmcl)) # stores AIC values in a vector
### ---- BEST RESULT FOR dCm/dCL -> CUBIC MODEL
##-----------------------------------------------
test_full_cm0 <- lmer(Cm_AC ~ poly(Angle.True,1) + TI + (1|WingID), data = fit_results, REML = F)
test_full_cm02 <- lmer(Cm_AC ~ poly(Angle.True,2) + TI + (1|WingID), data = fit_results, REML = F)
test_full_cm03 <- lmer(Cm_AC ~ poly(Angle.True,3) + TI + (1|WingID), data = fit_results, REML = F)
x[5,] <- c(AIC(test_full_cm0), AIC(test_full_cm02), AIC(test_full_cm03)) # stores AIC values in a vector
### ---- BEST RESULT FOR dCm/dCL -> LINEAR MODEL

for (i in 1:5){
x[i,4:6] <- x[i,1:3] - min(x[i,1:3])   # AIC differences
x[i,7:9]   <- exp(-0.5 * x[i,4:6])       # relative likelihoods of models
x[i,10:12] <- x[i,7:9]/sum(x[i,7:9])}                   # Akaike weights 

