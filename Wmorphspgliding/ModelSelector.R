## ---- MODEL SELECTOR -----
####caret random forests#####
library(caret)
library(randomForest)
library(tidyr)
setwd(working_directory) 
source('distancecalc.R')

##Run all files in Main Script First.

#Predefine lists
nomodels = 21
predict.elb        <- vector("list",nomodels)
predict.elb_wt     <- vector("list",nomodels)
predict.invivo     <- vector("list",nomodels)
error              <- vector("list",nomodels)
error_wt           <- vector("list",nomodels)
error_wtvivo       <- vector("list",nomodels)
sqerror            <- vector("list",nomodels)
sqerror_wt         <- vector("list",nomodels)
sqerror_wtvivo     <- vector("list",nomodels)
windmodel          <- vector("list",nomodels)
gustmodel          <- vector("list",nomodels)
predict.elb        <- vector("list",nomodels)
modelerr           <- data.frame(matrix(nrow = nomodels, ncol = 8))
elbow.range        <- data.frame(matrix(nrow = len_invivo, ncol = nomodels))
names(modelerr)    <- c("ID","modelname","meanerr_wt","meanerr_lab","meanerr_wtvivo",
                        "meansqerr_wt","meansqerr_lab","meansqerr_wtvivo")
modelerr$ID        <- letters[1:nomodels]
modelerr$modeltype <- "Linear - 2 PCs additive"


ctrl      <-trainControl(method = "cv",number = 10) #set up methods for cross-validation

## ------- All tested models -------------
#----- Linear models -----
modelerr$modeltype[1:6] <- "Linear"

#2 PCs
model.lm.2PC.add     <- train(elbow.angle~PC1+PC2, data = pc.labman, method = "lm", trControl = ctrl)
predict.elb[[1]]     <- predict(model.lm.2PC.add, newdata = pc_results_labman)
predict.elb_wt[[1]]  <- predict(model.lm.2PC.add, newdata = pc_results_wt)
predict.invivo[[1]]  <- predict(model.lm.2PC.add, newdata = pc_results_invivo[,7:16])
modelerr$modelname[1]<- "Linear - 2 PCs additive"

model.lm.2PC.int     <- train(elbow.angle~PC1*PC2, data = pc.labman, method = "lm", trControl = ctrl)
predict.elb[[2]]     <- predict(model.lm.2PC.int, newdata = pc_results_labman)
predict.elb_wt[[2]]  <- predict(model.lm.2PC.int, newdata = pc_results_wt)
predict.invivo[[2]]  <- predict(model.lm.2PC.int, newdata = pc_results_invivo[,7:16])
modelerr$modelname[2]<- "Linear - 2 PCs interative"

#4 PCs
model.lm.4PC.add     <- train(elbow.angle~PC1+PC2+PC3+PC4, data = pc.labman, method = "lm",trControl = ctrl)
predict.elb[[3]]     <- predict(model.lm.4PC.add, newdata = pc_results_labman)
predict.elb_wt[[3]]  <- predict(model.lm.4PC.add, newdata = pc_results_wt)
predict.invivo[[3]]  <- predict(model.lm.4PC.add,    newdata = pc_results_invivo[,7:16])
modelerr$modelname[3]<- "Linear - 4 PCs additive"

model.lm.4PC.int     <- train(elbow.angle~PC1*PC2*PC3*PC4, data = pc.labman, method = "lm",trControl = ctrl)
predict.elb[[4]]     <- predict(model.lm.4PC.int, newdata = pc_results_labman)
predict.elb_wt[[4]]  <- predict(model.lm.4PC.int, newdata = pc_results_wt)
predict.invivo[[4]]  <- predict(model.lm.4PC.int, newdata = pc_results_invivo[,7:16])
modelerr$modelname[4]<- "Linear - 4 PCs interative"

#10 PCs
model.lm.10PC.add    <- train(elbow.angle~PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = pc.labman, method = "lm",trControl = ctrl)
predict.elb[[5]]     <- predict(model.lm.10PC.add, newdata = pc_results_labman)
predict.elb_wt[[5]]  <- predict(model.lm.10PC.add, newdata = pc_results_wt)
predict.invivo[[5]]  <- predict(model.lm.10PC.add, newdata = pc_results_invivo[,7:16])
modelerr$modelname[5]<- "Linear - 10 PCs additive"

model.lm.10PC.int    <- train(elbow.angle~PC1 * PC2 * PC3 * PC4 * PC5 * PC6 * PC7 * PC8 * PC9 * PC10, data = pc.labman, method = "lm",trControl = ctrl)
predict.elb[[6]]     <- predict(model.lm.10PC.int, newdata = pc_results_labman)
predict.elb_wt[[6]]  <- predict(model.lm.10PC.int, newdata = pc_results_wt)
predict.invivo[[6]]  <- predict(model.lm.10PC.int, newdata = pc_results_invivo[,7:16])
modelerr$modelname[6]<- "Linear - 10 PCs interative"

#----- Loess models
modelerr$modeltype[7:8] <- "Loess"
loessGrid<-expand.grid(span = seq(0.1,0.9,0.1), degree = 1)

model.loess.2PC.add  <- train(elbow.angle~PC1+PC2, data = pc.labman, method = "gamLoess",trControl = ctrl, tuneGrid=loessGrid)
predict.elb[[7]]     <- predict(model.loess.2PC.add, newdata = pc_results_labman)
predict.elb_wt[[7]]  <- predict(model.loess.2PC.add, newdata = pc_results_wt)
predict.invivo[[7]]  <- predict(model.loess.2PC.add, newdata = pc_results_invivo[,7:16])
modelerr$modelname[7]<- "Loess - 2 PCs additive"

model.loess.4PC.add  <- train(elbow.angle~PC1+PC2+PC3+PC4, data = pc.labman, method = "gamLoess",trControl = ctrl, tuneGrid=loessGrid)
predict.elb[[8]]     <- predict(model.loess.4PC.add, newdata = pc_results_labman)
predict.elb_wt[[8]]  <- predict(model.loess.4PC.add, newdata = pc_results_wt)
predict.invivo[[8]]  <- predict(model.loess.4PC.add, newdata = pc_results_invivo[,7:16])
modelerr$modelname[8]<- "Loess - 4 PCs additive"

# ---- Linear Polynomial models
modelerr$modeltype[9:12] <- "Linear"

model.poly2            <- lm(elbow.angle ~ poly(PC2,2) + poly(PC1,2), data = pc_results_labman)
predict.elb[[9]]      <- predict(model.poly2, newdata = pc_results_labman)
predict.elb_wt[[9]]   <- predict(model.poly2, newdata = pc_results_wt)
predict.invivo[[9]]   <- predict(model.poly2, newdata = pc_results_invivo[,7:16])
modelerr$modelname[9] <- "Linear Polynomial PC1^2 + PC2^2"

model.poly3            <- lm(elbow.angle ~ poly(PC2,3) + poly(PC1,3), data = pc_results_labman)
predict.elb[[10]]      <- predict(model.poly3, newdata = pc_results_labman)
predict.elb_wt[[10]]   <- predict(model.poly3, newdata = pc_results_wt)
predict.invivo[[10]]   <- predict(model.poly3, newdata = pc_results_invivo[,7:16])
modelerr$modelname[10] <- "Linear Polynomial PC1^3 + PC2^3"

model.cos              <- lm(elbow.angle ~ cos(PC2) + cos(PC1), data = pc_results_labman)
predict.elb[[11]]      <- predict(model.cos, newdata = pc_results_labman)
predict.elb_wt[[11]]   <- predict(model.cos, newdata = pc_results_wt)
predict.invivo[[11]]   <- predict(model.cos, newdata = pc_results_invivo[,7:16])
modelerr$modelname[11] <- "Linear Cosine"

model.sin              <- lm(elbow.angle ~ sin(PC2) + sin(PC1), data = pc_results_labman)
predict.elb[[12]]      <- predict(model.sin, newdata = pc_results_labman)
predict.elb_wt[[12]]   <- predict(model.sin, newdata = pc_results_wt)
predict.invivo[[12]]   <- predict(model.sin, newdata = pc_results_invivo[,7:16])
modelerr$modelname[12] <- "Linear Sinusoidal"


modelerr$modelname[13]    <- "Closest Neighbor - 2PCs"
modelerr$modelname[14]    <- "Avg of 3 Closest Neighbors - 2PCs"
modelerr$modelname[15]    <- "Closest Neighbor - 6PCs"
modelerr$modelname[16]    <- "Avg of 3 Closest Neighbors - 6PCs"
modelerr$modelname[17]    <- "Closest Neighbor - 10PCs"
modelerr$modelname[18]    <- "Avg of 3 Closest Neighbors - 10PCs"

modelerr$modeltype[13:18] <- "Closest Neighbor"

###----- Build distance calculator to predict elbow angle based on closest
for(i in 1:36){
  check <- distancecalc(pc_results_wt[i,7],pc_results_wt[i,8],pc_results_wt[i,9],pc_results_wt[i,10],pc_results_wt[i,11],
                        pc_results_wt[i,12],pc_results_wt[i,13],pc_results_wt[i,14],pc_results_wt[i,15],pc_results_wt[i,16])
  
  predict.elb_wt[[13]][i]      <- check$closest_first2
  predict.elb_wt[[14]][i]      <- check$closest_first2_avg3
  predict.elb_wt[[15]][i]      <- check$closest_first6
  predict.elb_wt[[16]][i]      <- check$closest_first6_avg3
  predict.elb_wt[[17]][i]      <- check$closest_all10
  predict.elb_wt[[18]][i]      <- check$closest_all10_avg3

}

for (i in 1:len_labman){
check_full <- distancecalc(pc_results_labman[i,7],pc_results_labman[i,8],pc_results_labman[i,9],pc_results_labman[i,10],pc_results_labman[i,11],
                           pc_results_labman[i,12],pc_results_labman[i,13],pc_results_labman[i,14],pc_results_labman[i,15],pc_results_labman[i,16])

predict.elb[[13]][i]         <- check_full$closest_first2
predict.elb[[14]][i]         <- check_full$closest_first2_avg3
predict.elb[[15]][i]         <- check_full$closest_first6
predict.elb[[16]][i]         <- check_full$closest_first6_avg3
predict.elb[[17]][i]         <- check_full$closest_all10
predict.elb[[18]][i]         <- check_full$closest_all10_avg3
}


for (i in 1:len_invivo){
  check_invivo <- distancecalc(pc_results_invivo[i,7],pc_results_invivo[i,8],pc_results_invivo[i,9],pc_results_invivo[i,10],pc_results_invivo[i,11],
                             pc_results_invivo[i,12],pc_results_invivo[i,13],pc_results_invivo[i,14],pc_results_invivo[i,15],pc_results_invivo[i,16])
  
  predict.invivo[[13]][i]         <- check_invivo$closest_first2
  predict.invivo[[14]][i]         <- check_invivo$closest_first2_avg3
  predict.invivo[[15]][i]         <- check_invivo$closest_first6
  predict.invivo[[16]][i]         <- check_invivo$closest_first6_avg3
  predict.invivo[[17]][i]         <- check_invivo$closest_all10
  predict.invivo[[18]][i]         <- check_invivo$closest_all10_avg3
}


## ----- Random Forest Models ---------
modelerr$modeltype[19:21] <- "Random Forest"

pc.labman <-pc_results_labman[,c(4,5,7:16)]


#rf.full will contain the random forests model that uses PCs 1-10 as predictors of elbow angle
rf.10pc                <- train(elbow.angle~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=pc.labman,method='rf',trControl=ctrl)
predict.elb[[19]]      <- predict(rf.10pc, newdata=pc.labman)
predict.elb_wt[[19]]   <- predict(rf.10pc, newdata=pc_results_wt)
predict.invivo[[19]]   <- predict(rf.10pc, newdata = pc_results_invivo[,7:16])
modelerr$modelname[19] <- "Random Forest 10PCs"

rf.4pc                  <- train(elbow.angle~PC1+PC2+PC3+PC4,data=pc.labman,method='rf',trControl=ctrl)
predict.elb[[20]]       <- predict(rf.4pc, newdata=pc.labman)
predict.elb_wt[[20]]    <- predict(rf.4pc, newdata=pc_results_wt)
predict.invivo[[20]]    <- predict(rf.4pc, newdata = pc_results_invivo[,7:16])
modelerr$modelname[20]  <- "Random Forest 4PCs"

rf.2pc                  <- train(elbow.angle~PC1+PC2,data=pc.labman,method='rf',trControl=ctrl)
predict.elb[[21]]       <- predict(rf.2pc, newdata=pc.labman)
predict.elb_wt[[21]]    <- predict(rf.2pc, newdata=pc_results_wt)
predict.invivo[[21]]    <- predict(rf.2pc, newdata = pc_results_invivo[,7:16])
modelerr$modelname[21]  <- "Random Forest 2PCs"

### ---- Calculate Mean Square Error ------
for (i in 1:nomodels){
  #Absolute Error
  error[[i]]        <- predict.elb[[i]]    - pc_results_labman$elbow.angle
  error_wt[[i]]     <- predict.elb_wt[[i]] - pc_results_wt$elbow.angle
  error_wtvivo[[i]] <- error_wt[[i]][which(pc_results_wt$elbow.angle>80)]
  #Square Error
  sqerror[[i]]        <- (predict.elb[[i]]    - pc_results_labman$elbow.angle)^2
  sqerror_wt[[i]]     <- (predict.elb_wt[[i]] - pc_results_wt$elbow.angle)^2
  sqerror_wtvivo[[i]] <- (error_wt[[i]][which(pc_results_wt$elbow.angle>80)])^2
  #Mean Absolute Error
  modelerr$meanerr_lab[i]        <- mean(abs(error[[i]]))
  modelerr$meanerr_wt[i]     <- mean(abs(error_wt[[i]]))
  modelerr$meanerr_wtvivo[i] <- mean(abs(error_wtvivo[[i]]))
  #Mean Square Error
  modelerr$meansqerr_lab[i]        <- mean(sqerror[[i]])
  modelerr$meansqerr_wt[i]     <- mean(sqerror_wt[[i]])
  modelerr$meansqerr_wtvivo[i] <- mean(sqerror_wtvivo[[i]])
  #Root mean Square Error
  modelerr$rmserr_lab[i]        <- sqrt(mean(sqerror[[i]]))
  modelerr$rmserr_wt[i]     <- sqrt(mean(sqerror_wt[[i]]))
  modelerr$rmserr_wtvivo[i] <- sqrt(mean(sqerror_wtvivo[[i]]))
  #Elbow range per model prediction
  elbow.range[,i]           <- predict.invivo[[i]]
  #Wind Model Goodness of fit
  windmodel[[i]]               <- lm(predict.invivo[[i]]~pc_results_invivo$Wind.Speed)
  wind                         <- summary(windmodel[[i]])
  modelerr$coeff.wind[i]       <- wind$coefficients[2]
  modelerr$confint2.5.wind[i]  <- confint(windmodel[[i]])[2]
  modelerr$confint97.5.wind[i] <- confint(windmodel[[i]])[4]
  
  #Gust Model Goodness of fit
  gustmodel[[i]]               <- lm(predict.invivo[[i]]~pc_results_invivo$Max.Wind.Gust)
  gust                         <- summary(gustmodel[[i]])
  modelerr$coeff.gust[i]       <- gust$coefficients[2]
  modelerr$confint2.5.gust[i]  <- confint(gustmodel[[i]])[2]
  modelerr$confint97.5.gust[i] <- confint(gustmodel[[i]])[4]
}
#convert elbow range to long form for plotting
names(elbow.range)         <- letters[1:nomodels]
elbow.range.long           <- gather(elbow.range,ID,elbow.angle,a:letters[nomodels])
elbow.range.long$modeltype <- modelerr$modeltype[match(elbow.range.long$ID, modelerr$ID)]

modelerr$score = sqrt(modelerr$rmserr_lab^2 + (modelerr$rmserr_wt)^2 + (modelerr$rmserr_wtvivo)^2)

##-----------------------------------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------------------------------

rf.4pc.man        <-train(manus.angle~PC1+PC2+PC3+PC4,data=pc.labman,method='rf',trControl=ctrl)
predict.man       <-predict(rf.4pc.man,newdata=pc.labman)
predict.man_wt    <-predict(rf.4pc.man,newdata=pc_results_wt)

sqerror.man      <- (predict.man - pc_results_labman$manus.angle)^2
meansqerr.man    <- mean(sqerror.man)/length(pc_results_labman)
sqerror_wt.man   <- (predict.man_wt - pc_results_wt$manus.angle)^2



  