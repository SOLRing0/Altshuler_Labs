########## Geometric Morphometric Projection-izer - Planform ##########
### Created: Vikram Baliga 2017-06-26 
### Revised: Christina Harvey 2017-06-28
### Using geomorph v. 3.0.3

### This code:
### 1) imports landmark data from lab manipulations, invivo and wind tunnel(camber projection or planform)
### 2) performs Procrustes superimposition with all data and creates wireframe
### 3) seperates the data and performs PCA on only the labmanipulation
### 4) using the eigenvectors projects the invivo and windtunnel data onto the labmanipulation PCA 

#### Define Libraries #####

library(geomorph)  #geometric morphometrics
library(abind)     #We will use this later to combine arrays
library(gridBase)  #plotting
library(ggplot2)   #plotting
library(gridExtra) #plotting
library(reshape2)  #melt for loes fit
library(stringr)
library(caret)
library(randomForest)
#### Read in Previous Info #####
setwd(working_directory)
wtangles    <- read.csv('windtunnelangles.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
photodat    <- read.csv('photometadata.csv', stringsAsFactors = FALSE,strip.white = TRUE, na.strings = c("") )
photodat$ID <- paste(photodat$Date,"_",photodat$Location,"_",photodat$PhotoID, sep = "")

#### Initialize any variables #####
p<-5 # number of landmarks 
k<-2 # number of dimensions

##### Importing Data #####
#pull in lab manipulation - Prosillica 
rawdata                        <- read.csv(file.choose(),header = TRUE, stringsAsFactors=FALSE,strip.white = TRUE, na.strings = c("") ) # landmarks.txt file from Christina
rawdata                        <- rawdata[,-c(8:11)]
names                          <- rawdata$ID
testid                         <- rawdata$TestID
rawdata.three.d                <- arrayspecs(rawdata[,8:17], 5, k) # create 3D array from data
dimnames(rawdata.three.d)[[3]] <- paste(names,testid,sep="")
len_labman                     <- length(rawdata.three.d[1,1,])
rawdata.three.d[,2,]           <- -rawdata.three.d[,2,]  #flips about y axis

#pull in Christina's invivo digitized data
tmp                      <- read.table(file.choose(),header = TRUE, stringsAsFactors = FALSE )
names.InVivo.CH          <- tmp[,1]
InVivo.CH                <- arrayspecs(tmp[,2:ncol(tmp)], 7, k)
dimnames(InVivo.CH)[[3]] <- names.InVivo.CH 
len_invivo               <- length(InVivo.CH[1,1,])
InVivo.CH                <- InVivo.CH[-c(1,2),,] #removes lm 1 and 2

#pull in wind tunnel outlines
tmp                       <- read.table(file.choose(),header = TRUE, stringsAsFactors = FALSE )
names.WT.CH               <- tmp[,1]
position.WT               <- tmp[,16]
WT.CH                     <- arrayspecs(tmp[,2:15], 7, k)
dimnames(WT.CH)[[3]]      <- paste(names.WT.CH,position.WT,sep = "")
len_WT                    <- length(WT.CH[1,1,])
WT.CH                     <- WT.CH[-c(1,2),,] #removes lm 1 and 2

##### Combine the datasets#####
#NOTE: it is important for the error analysis that rawdata.three.d always is input first
combined.landmarks <- abind(rawdata.three.d,InVivo.CH,WT.CH) 

###------ PCA ANALYSIS
##### Superimposition #####
combined.super<-gpagen(combined.landmarks) # Procrustes superimposition #plot(combined.super)
combined.size <-combined.super$Csize       # Centroid size of each specimen

##### 2D and 3D arrays of coordinates #####
combined.3D  <-combined.super$coords       #3D array of proc. superimposed coordinates
combined.2D  <-two.d.array(combined.3D)    #create 2D array

##### Data Splitting #####
labman.2D    <- combined.2D[1:len_labman, ]  
InVivo.2D    <- combined.2D[(len_labman+1):(len_labman+len_invivo), ] 
WT.2D        <- combined.2D[(len_labman+len_invivo+1):(len_labman+len_invivo+len_WT), ] 

labman.3D    <- arrayspecs(labman.2D,p,k) # 3D array

##### PCA on data subset #####
labman.PCA   <- plotTangentSpace(labman.3D, tol = 0) 

##------- TRANSFORM WT AND INVIVO DATA
##### Calculate PC scores of unretained data #####
C              <- labman.3D
ref.unret      <- mshape(C)
x.unret        <- c(ref.unret[1,1],ref.unret[1,2],ref.unret[2,1],ref.unret[2,2],ref.unret[3,1],ref.unret[3,2],ref.unret[4,1],ref.unret[4,2],ref.unret[5,1],ref.unret[5,2])
InVivo.Updated <- InVivo.2D
WT.Updated = WT.2D
#--- Remove mean shape from data
for (j in 1:nrow(InVivo.2D)){
  InVivo.Updated[j,] = InVivo.2D[j,] - x.unret
}
for (j in 1:nrow(WT.2D)){
  WT.Updated[j,] = WT.2D[j,] - x.unret
}
#####  Multiply by Eigenvectors ##### 
InVivo.points <- (InVivo.Updated) %*% labman.PCA$rotation 
WT.points     <- (WT.Updated) %*% labman.PCA$rotation 
labman.points <- labman.PCA$pc.scores

##### --------  Create Master DataFrame ----------- ##### 
tmp               <- rbind(labman.points,InVivo.points)
tmp               <- rbind(tmp,WT.points)
pc_results        <- data.frame(matrix(nrow = nrow(tmp), ncol = 8))
names(pc_results) <- c("ID","dataset","GullID","elbow.angle","manus.angle","digit.angle","PC1","PC2")

#PC RESULTS
pc_results$PC1 <- - tmp[,1] #NEGATIVE ADDED JUST TO FLIP THE GRAPH For ease of plotting
pc_results$PC2 <- tmp[,2]
pc_results$PC3 <- tmp[,3]
pc_results$PC4 <- tmp[,4]
pc_results$PC5 <- tmp[,5]
pc_results$PC6 <- tmp[,6]
pc_results$PC7 <- tmp[,7]
pc_results$PC8 <- tmp[,8]
pc_results$PC9 <- tmp[,9]
pc_results$PC10 <- tmp[,10]

# ID 
pc_results$ID[1:len_labman]                                              <- names
pc_results$ID[(len_labman+1):(len_labman+len_invivo)]                    <- names.InVivo.CH
pc_results$ID[(len_labman+len_invivo+1):(len_labman+len_invivo+len_WT)]  <- names.WT.CH

# DATASET
pc_results$dataset[1:len_labman]                                              <- "labman"
pc_results$dataset[(len_labman+1):(len_labman+len_invivo)]                    <- "invivo"
pc_results$dataset[(len_labman+len_invivo+1):(len_labman+len_invivo+len_WT)]  <- "windtunnel"
# GULLID
pc_results$GullID[1:len_labman]                                               <- rawdata$GullID

# ELBOW ANGLE
pc_results$elbow.angle[1:len_labman]                            <- rawdata$elbow.angle # prosillica
pc_results$elbow.angle[(len_labman+1):(len_labman+len_invivo)]  <- NA #invivo
tmp <- wtangles$elbow.angle[match(pc_results$ID[(len_labman+len_invivo+1):(len_labman+len_invivo+len_WT)], wtangles$WingID)]
pc_results$elbow.angle[(len_labman+len_invivo+1):(len_labman+len_invivo+len_WT)]  <- tmp
pc_results$elbow.angle <- as.numeric(pc_results$elbow.angle)

# MANUS ANGLE
pc_results$manus.angle[1:len_labman]                            <- rawdata$manus.angle # prosillica
pc_results$manus.angle[(len_labman+1):(len_labman+len_invivo)]  <- NA  #invivo
tmp <- wtangles$manus.angle[match(pc_results$ID[(len_labman+len_invivo+1):(len_labman+len_invivo+len_WT)], wtangles$WingID)]
pc_results$manus.angle[(len_labman+len_invivo+1):(len_labman+len_invivo+len_WT)]  <- tmp
pc_results$manus.angle <- as.numeric(pc_results$manus.angle)

# DIGIT ANGLE
pc_results$digit.angle[1:len_labman]  <- rawdata$digit.angle
pc_results$orientation[(len_labman+len_invivo+1):(len_labman+len_invivo+len_WT)] <- position.WT

#- subset data
pc_results_invivo <- subset(pc_results, dataset %in% "invivo" )
pc_results_wt     <- subset(pc_results,dataset %in% "windtunnel")
pc_results_labman <- subset(pc_results,dataset %in% "labman")

#### ---------- Fit Model Selected Prior in "ModelSelector.R"
pc.labman <-pc_results_labman[,c(4,5,7:16)]
ctrl      <-trainControl(method = "cv",number = 10) #set up methods for cross-validation

rf.4pc            <- train(elbow.angle~PC1+PC2+PC3+PC4,data=pc.labman,method='rf',trControl=ctrl)
rf.4pc.man        <- train(manus.angle~PC1+PC2+PC3+PC4,data=pc.labman,method='rf',trControl=ctrl)
pc_results_invivo$elbow.angle <- predict(rf.4pc,     newdata = pc_results_invivo[,7:10])
pc_results_invivo$manus.angle <- predict(rf.4pc.man, newdata = pc_results_invivo[,7:10])

pc_results_invivo <- merge(pc_results_invivo, photodat[,c("ID","Date","Location","Time","Temperature","Wind.Speed","Max.Wind.Gust","Glide.Style")], by = c("ID"), all.x = TRUE, all.y = FALSE)

#update units from km/h to m/s
pc_results_invivo$Wind.Speed    = pc_results_invivo$Wind.Speed/3.6 
pc_results_invivo$Max.Wind.Gust = pc_results_invivo$Max.Wind.Gust/3.6 


