###----- This code checks the sensitivity of the photos to position of wings
setwd(working_directory)
p<-5 # number of landmarks 
k<-2 # number of dimensions

sensitivity            <- read.table('sensitivity_landmarks.txt', header = TRUE, stringsAsFactors = FALSE)
names.sens.CH          <- sensitivity[,1]
sens.CH                <- arrayspecs(sensitivity[,5:ncol(sensitivity)], 7, k)
dimnames(sens.CH)[[3]] <- names.sens.CH 
len_sens               <- length(sens.CH[1,1,])
sens.CH                <- sens.CH[-c(1,2),,] #removes lm 1 and 2
sens.CH[,2,]           <- -sens.CH[,2,]  #flips about y axis

##### Combine the datasets#####
#NOTE: it is important for the error analysis that rawdata.three.d always is input first
combined.landmarks.sens <- abind(rawdata.three.d,InVivo.CH,WT.CH,sens.CH) 

###------ PCA ANALYSIS
##### Superimposition #####
combined.super.sens<-gpagen(combined.landmarks.sens) # Procrustes superimposition #plot(combined.super)
combined.size.sens <-combined.super.sens$Csize       # Centroid size of each specimen

##### 2D and 3D arrays of coordinates #####
combined.3D.sens  <-combined.super.sens$coords       #3D array of proc. superimposed coordinates
combined.2D.sens  <-two.d.array(combined.3D.sens)    #create 2D array

##### Data Splitting #####
labman.2D.sens    <- combined.2D.sens[1:len_labman, ]  # the remaining 91 rows, all lab manipulations
InVivo.2D.sens    <- combined.2D.sens[(len_labman+1):(len_labman+len_invivo), ] 
WT.2D.sens        <- combined.2D.sens[(len_labman+len_invivo+1):(len_labman+len_invivo+len_WT), ] 
sens.2D.sens      <- combined.2D.sens[(len_labman+len_invivo+len_WT+1):(len_labman+len_invivo+len_WT+len_sens), ]

labman.3D.sens    <- arrayspecs(labman.2D.sens,p,k) # 3D array

##### PCA on data subset #####
labman.PCA.sens   <- plotTangentSpace(labman.3D.sens, tol = 0) # Should look very similar to alldata.PCA

##------- TRANSFORM WT AND INVIVO DATA
##### Calculate PC scores of unretained data #####
C.sens              <- labman.3D.sens
ref.unret.sens      <- mshape(C.sens)
x.unret.sens        <- c(ref.unret.sens[1,1],ref.unret.sens[1,2],ref.unret.sens[2,1],ref.unret.sens[2,2],ref.unret.sens[3,1],ref.unret.sens[3,2],ref.unret.sens[4,1],ref.unret.sens[4,2],ref.unret.sens[5,1],ref.unret.sens[5,2])
InVivo.Updated.sens <- InVivo.2D.sens
WT.Updated.sens     <- WT.2D.sens
sens.Updated.sens   <- sens.2D.sens

#--- Remove mean shape from data
for (j in 1:nrow(InVivo.2D.sens)){
  InVivo.Updated.sens[j,] = InVivo.2D.sens[j,] - x.unret.sens
}
for (j in 1:nrow(WT.2D.sens)){
  WT.Updated.sens[j,] = WT.2D.sens[j,] - x.unret.sens
}
for (j in 1:nrow(sens.2D.sens)){
  sens.Updated.sens[j,] = sens.2D.sens[j,] - x.unret.sens
}

#####  Multiply by Eigenvectors ##### 
InVivo.points.sens <- (InVivo.Updated.sens) %*% labman.PCA.sens$rotation 
WT.points.sens     <- (WT.Updated.sens) %*% labman.PCA.sens$rotation 
sens.points.sens   <- (sens.Updated.sens) %*% labman.PCA.sens$rotation 
labman.points.sens <- labman.PCA.sens$pc.scores

##### --------  Create Master DataFrame ----------- ##### 
tmp               <- rbind(labman.points.sens,InVivo.points.sens)
tmp               <- rbind(tmp,WT.points.sens)
tmp               <- rbind(tmp,sens.points.sens)
pc_results.sens        <- data.frame(matrix(nrow = nrow(tmp), ncol = 8))
names(pc_results.sens) <- c("ID","dataset","GullID","elbow.angle","manus.angle","digit.angle","PC1","PC2")

#PC RESULTS
pc_results.sens$PC1 <- tmp[,1]
pc_results.sens$PC2 <- tmp[,2]
pc_results.sens$PC3 <- tmp[,3]
pc_results.sens$PC4 <- tmp[,4]
pc_results.sens$PC5 <- tmp[,5]
pc_results.sens$PC6 <- tmp[,6]
pc_results.sens$PC7 <- tmp[,7]
pc_results.sens$PC8 <- tmp[,8]
pc_results.sens$PC9 <- tmp[,9]
pc_results.sens$PC10 <- tmp[,10]

# ID 
pc_results.sens$ID[1:len_labman]                                              <- names
pc_results.sens$ID[(len_labman+1):(len_labman+len_invivo)]                    <- names.InVivo.CH
pc_results.sens$ID[(len_labman+len_invivo+1):(len_labman+len_invivo+len_WT)]  <- names.WT.CH
pc_results.sens$ID[(len_labman+len_invivo+len_WT+1):(len_labman+len_invivo+len_WT+len_sens)]  <- names.sens.CH

# DATASET
pc_results.sens$dataset[1:len_labman]                                                                <- "labman"
pc_results.sens$dataset[(len_labman+1):(len_labman+len_invivo)]                                      <- "invivo"
pc_results.sens$dataset[(len_labman+len_invivo+1):(len_labman+len_invivo+len_WT)]                    <- "windtunnel"
pc_results.sens$dataset[(len_labman+len_invivo+len_WT+1):(len_labman+len_invivo+len_WT+len_sens)]    <- "sensitivity"
# GULLID
pc_results.sens$GullID[1:len_labman]                                               <- rawdata$GullID
pc_results.sens$GullID[(len_labman+len_invivo+len_WT+1):(len_labman+len_invivo+len_WT+len_sens)] <- sensitivity$GullID
# ELBOW ANGLE
pc_results.sens$elbow.angle[1:len_labman]                            <- rawdata$elbow.angle # prosillica
pc_results.sens$elbow.angle[(len_labman+1):(len_labman+len_invivo)]  <- NA #invivo
tmp <- wtangles$elbow.angle[match(pc_results.sens$ID[(len_labman+len_invivo+1):(len_labman+len_invivo+len_WT)], wtangles$WingID)]
pc_results.sens$elbow.angle[(len_labman+len_invivo+1):(len_labman+len_invivo+len_WT)]  <- tmp
pc_results.sens$elbow.angle <- as.numeric(pc_results.sens$elbow.angle)

# MANUS ANGLE
pc_results.sens$manus.angle[1:len_labman]                            <- rawdata$manus.angle # prosillica
pc_results.sens$manus.angle[(len_labman+1):(len_labman+len_invivo)]  <- NA  #invivo
tmp <- wtangles$manus.angle[match(pc_results.sens$ID[(len_labman+len_invivo+1):(len_labman+len_invivo+len_WT)], wtangles$WingID)]
pc_results.sens$manus.angle[(len_labman+len_invivo+1):(len_labman+len_invivo+len_WT)]  <- tmp
pc_results.sens$manus.angle <- as.numeric(pc_results.sens$manus.angle)

# DIGIT ANGLE
pc_results.sens$digit.angle[1:len_labman]  <- rawdata$digit.angle
pc_results.sens$orientation                <- NA
pc_results.sens$orientation[(len_labman+len_invivo+1):(len_labman+len_invivo+len_WT)] <- position.WT

#- subset data
pc_results_sens <- subset(pc_results.sens, dataset %in% "sensitivity" )
pc_results_lab_sens <- subset(pc_results.sens, dataset %in% "labman" )
#### ---------- Fit Model Selected Prior in "ModelSelector.R"
pc.labman.sens <-pc.labman<-pc_results_lab_sens[,c(4,5,7:16)]
rf.4pc.sens             <-train(elbow.angle~PC1+PC2+PC3+PC4,data=pc.labman.sens,method='rf',trControl=ctrl)
rf.4pc.man.sens         <-train(manus.angle~PC1+PC2+PC3+PC4,data=pc.labman.sens,method='rf',trControl=ctrl)

pc_results_sens$elbow.angle <- predict(rf.4pc.sens,    newdata = pc_results_sens[,7:10])
pc_results_sens$manus.angle <- predict(rf.4pc.man.sens,    newdata = pc_results_sens[,7:10])

pc_results_sens$true.elbow.angle[pc_results_sens$GullID == "Gull3"] <- wtangles$elbow.angle[wtangles$WingID == "Gull3"] 
pc_results_sens$true.elbow.angle[pc_results_sens$GullID == "17_0243"] <- wtangles$elbow.angle[wtangles$WingID == "17_0243"] 
#pc_results_sens$elbow.angle <- predict(model.4pc,    newdata = pc_results_sens[,7:10])

pc_results_sens$true.manus.angle[pc_results_sens$GullID == "Gull3"] <- wtangles$manus.angle[wtangles$WingID == "Gull3"] 
pc_results_sens$true.manus.angle[pc_results_sens$GullID == "17_0243"] <- wtangles$manus.angle[wtangles$WingID == "17_0243"] 
#pc_results_sens$manus.angle <- predict(model.4pc.man,    newdata = pc_results_sens[,7:10])

pc_results_sens$rot.x   <- sensitivity$rot.x
pc_results_sens$rot.y   <- sensitivity$rot.y
pc_results_sens$err.elb.abs <- abs(pc_results_sens$elbow.angle - pc_results_sens$true.elbow.angle)
pc_results_sens$err.man.abs <- abs(pc_results_sens$manus.angle - pc_results_sens$true.manus.angle)
pc_results_sens$err.elb <- abs(pc_results_sens$elbow.angle - pc_results_sens$true.elbow.angle)/pc_results_sens$true.elbow.angle
pc_results_sens$err.man <- abs(pc_results_sens$manus.angle - pc_results_sens$true.manus.angle)/pc_results_sens$true.manus.angle

# create the dataframes that will be used for plotting
plotdf.y <- pc_results_sens[pc_results_sens$rot.y == 0,]
tmp <- pc_results_sens[(pc_results_sens$GullID == "Gull3" & pc_results_sens$rot.x == 0),]
tmp1 <- pc_results_sens[(pc_results_sens$GullID == "17_0243" & pc_results_sens$rot.x == -20),]
plotdf.x <- rbind(tmp,tmp1)
