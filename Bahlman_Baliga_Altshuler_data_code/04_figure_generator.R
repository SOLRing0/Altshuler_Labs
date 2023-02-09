# please be sure to run 01_allFunctions and 02_dataImport prior to running
# anything in this script

################################ figure prep work ##############################

# set up color pallete
#colz<-colorRampPalette(c("black", "red"))(20) ## (n)
colz<-viridis(20,option="D")
#plot(1:20,c(rep(0,20)),pch=15,cex=5,col=colz) # just to visualize
reddd<-viridis(20,option="A")[14] #"#A52E13"
# create color vectors
threecolz<-viridis(20,option="D")[c(1,10,15)]
#plot(1:3,c(rep(0,3)),pch=15,cex=5,col=threecolz)
sevencolz<-viridis(8,option="A")
sevencolz<-sevencolz[c(4,1,2,3,5,6,7)]
#plot(1:7,c(rep(0,7)),pch=15,cex=5,col=sevencolz)

################################### figure 1 ###################################
# read in the data file used to generate figure 1; from 2014-10-01
fig1dat<-read.ddf.WL("./raw_data/2014-10-01/many loops 20.ddf")
# now select all cycles
fig1cycles<-selectCycles_p2p(fig1dat,keep.cycles=1:20)
# and analyze to compute work & power output
fig1analyzed<-analyzeWorkLoop(fig1cycles,GR=2)


#### fig 1 panel F ####
par(mar=c(5,4,4,6)+0.3)
plot(fig1dat$Time[1:8000],fig1dat$Position[1:8000]/2,
     axes=FALSE,ylim=c(-2,3),tck=0.02,bty="n",ylab="",xlab="",
     pch=19,col="white")
lines(fig1dat$Time[1:8000],fig1dat$Position[1:8000]/2,
      lwd=1,ylab="",xlab="",lty=2)
tmp<-fig1dat[fig1dat$Stim==1,]
tmp$Position<-tmp$Position/2;tmp$Force<-tmp$Force*2/1000
a<-1:161;b<-a[seq(1,length(a),8)];c<-a[seq(1,length(a),4)]
for (i in b){
  rect(tmp$Time[i],-2,tmp$Time[i+7],
       2.5,col=rgb(0,0,0,alpha=0.1),border=NA)
}
axis(2,ylim=c(-2,3),tck=0.02,las=1)
mtext("position (mm)",side=2,line=2.5)
par(new=TRUE)
plot(fig1dat$Time[1:8000],(fig1dat$Force[1:8000]*2)/1000,
     type="n",
     ylim=c(0,3),axes=FALSE,ylab="",xlab="")
lines(fig1dat$Time[1:8000],(fig1dat$Force[1:8000]*2)/1000,
      lwd=2,col=reddd,ylab="",xlab="",type="l",lty=1) #lty =5
for (i in 1:20){
  rect(fig1analyzed$Time[i][[1]][1],2.75,
       fig1analyzed$Time[i][[1]][length(fig1analyzed$Time[i][[1]])],3,
       col=colz[i])
}
axis(side=4,c(0,1.5,3),tck=0.02,las=1,col=reddd,lwd=2)
mtext("force (N)",side=4,line=2.5)
axis(side=1,tck=0.02)
mtext("time (sec)",side=1,line=2.5)
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1),mgp=c(3,1,0),las=0) # reset


#### fig 1 panel G ####
plot((fig1analyzed$Position[20][[1]]),(fig1analyzed$Force[20][[1]]),
     ylim=c(500,2700),xlim=c(-1.5,2),type="n",
     tck=0.02,bty="n",ylab="force (mN)",xlab="position (mm)")
for(i in 1:length(fig1analyzed$Position)){
  lines((fig1analyzed$Position[i][[1]]),(fig1analyzed$Force[i][[1]]),
        col=colz[i],lwd=2)
}


#### fig 1 panel H ####
plot(1:20,ldply(fig1analyzed$Work)$V1,pch=21,col="black",bg=colz,
     xlim=c(0,20),ylim=c(0,2.5),cex=2,
     ylab="work (mJ)",xlab="cycle number",
     tck=0.02,bty="n")
rect(2.5,0,5.5,2.5,col=rgb(0,0,0,alpha=0.1),border=NA)

#### fig 1 panel I ####

## work loops are the top row and 
## instantaneous power are the second

par(mfrow=c(2,3))
for (i in 1:3){
  plot(fig1analyzed$Position[[i]],
       fig1analyzed$Force[[i]],
       xlim=c(-1.5,2),ylim=c(0,4000),
       type='n',
       xlab="position",
       ylab="force (mN)",
       tck=0.02,bty='n',las=1)
  lines(fig1analyzed$Position[[i]],
        fig1analyzed$Force[[i]],
        col=colz[i+2],lty=1,lwd=2)
}
for (i in 1:3){
  plot(seq(0,1,1/(length(fig1analyzed$InstantPower[[i]])-1)),
       fig1analyzed$InstantPower[[i]]/bird07_pecmass,
       ylim=c(-450,900),axes=FALSE,
       xlab="percent cycle",ylab="instantaneous power (W/kg)",
       type="n",tck=0.02,bty='n',las=1)
  axis(2,at=c(-450,0,450,900),tck=0.02,las=1)
  axis(1,at=c(0,0.5,1),tck=0.02,las=1)
  abline(h=0)
  rect(0.5,-450,
       1,900,
       col=rgb(0,0,0,alpha=0.05),border=NA) 
  lines(seq(0,1,1/(length(fig1analyzed$InstantPower[[i]])-1)),
        fig1analyzed$InstantPower[[i]]/bird07_pecmass,lwd=2,
        col=colz[i+2]) 

}
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1),mgp=c(3,1,0),las=0) # reset


###### fig 1 panels C-D ######
## Elements make by these plots were later used in Illustrator to compose 
## panels C and D
tmp<-exp11_bird03_alltrials[[1]]$Time[[3]][exp11_bird03_alltrials[[1]]$Stim[[3]]==1]
plot(exp11_bird03_alltrials[[1]]$Time[[3]],
     exp11_bird03_alltrials[[1]]$Position[[3]],
     type='n',
     xlab="time",
     ylab="position",
     tck=0.02,bty='n',las=1)
lines(exp11_bird03_alltrials[[1]]$Time[[3]],
      exp11_bird03_alltrials[[1]]$Position[[3]],
      lwd=2,lty=2)
rect(tmp[1],-1.5,tmp[8],2,
     col=rgb(0,0,0,alpha=0.1),border=NA)
## This grey rectangle is actually the stimulus time, not upstroke vs downstroke
plot(exp11_bird03_alltrials[[1]]$Time[[3]],
     exp11_bird03_alltrials[[1]]$Force[[3]],
     type='n',
     xlab="time",
     ylab="force",
     tck=0.02,bty='n',las=1)
lines(exp11_bird03_alltrials[[1]]$Time[[3]],
      exp11_bird03_alltrials[[1]]$Force[[3]],
      lwd=2)
rect(tmp[1],500,tmp[8],6000,
     col=rgb(0,0,0,alpha=0.1),border=NA)
## This grey rectangle is actually the stimulus time, not upstroke vs downstroke

plot(exp11_bird03_alltrials[[1]]$Position[[3]][1:179],
     exp11_bird03_alltrials[[1]]$Force[[3]][1:179],
     ylim=c(0,6000),xlim=c(-2,2),
     type='n',
     xlab="position",
     ylab="force",
     tck=0.02,bty='n',las=1)
lines(exp11_bird03_alltrials[[1]]$Position[[3]][1:179],
      exp11_bird03_alltrials[[1]]$Force[[3]][1:179],
      lwd=2)
plot(exp11_bird03_alltrials[[1]]$Position[[3]][180:358],
     exp11_bird03_alltrials[[1]]$Force[[3]][180:358],
     ylim=c(0,6000),xlim=c(-2,2),
     type='n',
     xlab="position",
     ylab="force",
     tck=0.02,bty='n',las=1)
lines(exp11_bird03_alltrials[[1]]$Position[[3]][180:358],
      exp11_bird03_alltrials[[1]]$Force[[3]][180:358],
      lwd=2)
plot(exp11_bird03_alltrials[[1]]$Position[[3]],
     exp11_bird03_alltrials[[1]]$Force[[3]],
     ylim=c(0,6000),xlim=c(-2,2),
     type='n',
     xlab="position",
     ylab="force",
     tck=0.02,bty='n',las=1)
lines(exp11_bird03_alltrials[[1]]$Position[[3]],
      exp11_bird03_alltrials[[1]]$Force[[3]],
      lwd=2)


################################### figure s1 ##################################
library(lubridate)

baselines<-list(bird01_baseline_degradation,bird02_baseline_degradation,
                bird03_baseline_degradation,bird04_baseline_degradation,
                bird05_baseline_degradation,bird06_baseline_degradation,
                bird07_baseline_degradation,bird08_baseline_degradation,
                bird09_baseline_degradation,bird10_baseline_degradation,
                bird11_baseline_degradation,bird12_baseline_degradation)

for (i in 1:length(baselines)){
  for (j in 1:nrow(baselines[[i]])){
    baselines[[i]]$elapsed[j]<-as.numeric(as.duration(baselines[[i]]$time[1] %--% baselines[[i]]$time[j]),"minutes")
  }
}

ldply(baselines)->baselinesdf

colzz<-c("#a6cee3",
         "#1f78b4",
         "#b2df8a",
         "#33a02c",
         "#fb9a99",
         "#e31a1c",
         "#fdbf6f",
         "#ff7f00",
         "#cab2d6",
         "#6a3d9a",
         "#ffff99",
         "#b15928")

plot(baselinesdf$elapsed,
     baselinesdf$percent,
     type="n",ylim=c(0,100),
     ylab="%",xlab="time elapsed (mins)",
     tck=0.02,bty='n',las=1)
for (i in 1:length(baselines)){
  points(baselines[[i]]$elapsed,baselines[[i]]$percent,
         col=colzz[i],pch=19)
  lines(baselines[[i]]$elapsed,baselines[[i]]$percent,
        col=colzz[i],lwd=2)
}


## not used, but you may also enjoy:
# par(mfrow=c(4,3))
# for (i in 1:length(baselines)){
#   plot(baselines[[i]]$elapsed,
#        baselines[[i]]$percent,
#        pch=19,ylim=c(0,100),
#        ylab="%",xlab="time elapsed (mins)",
#        tck=0.02,bty='n',las=1)
#   lines(baselines[[i]]$elapsed,baselines[[i]]$percent)
#   abline(h=80,lty=2,col='red')
# }
# par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1),mgp=c(3,1,0),las=0) # reset



################################### figure 2 ###################################
# panels A, B, C use data from experiments 1-3
t<-dplyr::filter(allsummarydat,
                 experiment=="exp01"|experiment=="exp02"|experiment=="exp03")
attach(t)
t<-t[order(pulses,bird),]
detach(t)

# assign line types
linezz<-c(1:3) # bird 1 is solid; bird 4 is dashed; bird 7 is dotted

#### fig 3 panel A ####
plot(t$pulses[t$exp_bird=="exp01_bird01"],
     t$expCorPower[t$exp_bird=="exp01_bird01"],
     xlim=c(1,7),ylim=c(-200,150),
     pch=19,col=threecolz[1],
     xlab="stimulus duration (# pulses at 300 Hz)",
     ylab="net power (W/kg)",
     tck=0.02,bty='n',las=1)
abline(h=0)
lines(t$pulses[t$exp_bird=="exp01_bird01"],
      t$expCorPower[t$exp_bird=="exp01_bird01"],
      lty=linezz[1],col=threecolz[1])
points(t$pulses[t$exp_bird=="exp02_bird01"],
       t$expCorPower[t$exp_bird=="exp02_bird01"],
       pch=19,col=threecolz[2])
lines(t$pulses[t$exp_bird=="exp02_bird01"],
      t$expCorPower[t$exp_bird=="exp02_bird01"],
      lty=linezz[1],col=threecolz[2])
points(t$pulses[t$exp_bird=="exp03_bird01"],
       t$expCorPower[t$exp_bird=="exp03_bird01"],
       pch=19,col=threecolz[3])
lines(t$pulses[t$exp_bird=="exp03_bird01"],
      t$expCorPower[t$exp_bird=="exp03_bird01"],
      lty=linezz[1],col=threecolz[3])
points(t$pulses[t$exp_bird=="exp01_bird04"],
       t$expCorPower[t$exp_bird=="exp01_bird04"],
       pch=19,col=threecolz[1])
lines(t$pulses[t$exp_bird=="exp01_bird04"],
      t$expCorPower[t$exp_bird=="exp01_bird04"],
      lty=linezz[2],col=threecolz[1])
points(t$pulses[t$exp_bird=="exp02_bird04"],
       t$expCorPower[t$exp_bird=="exp02_bird04"],
       pch=19,col=threecolz[2])
lines(t$pulses[t$exp_bird=="exp02_bird04"],
      t$expCorPower[t$exp_bird=="exp02_bird04"],
      lty=linezz[2],col=threecolz[2])
points(t$pulses[t$exp_bird=="exp03_bird04"],
       t$expCorPower[t$exp_bird=="exp03_bird04"],
       pch=19,col=threecolz[3])
lines(t$pulses[t$exp_bird=="exp03_bird04"],
      t$expCorPower[t$exp_bird=="exp03_bird04"],
      lty=linezz[2],col=threecolz[3])
points(t$pulses[t$exp_bird=="exp01_bird07"],
       t$expCorPower[t$exp_bird=="exp01_bird07"],
       pch=19,col=threecolz[1])
lines(t$pulses[t$exp_bird=="exp01_bird07"],
      t$expCorPower[t$exp_bird=="exp01_bird07"],
      lty=linezz[3],col=threecolz[1])
points(t$pulses[t$exp_bird=="exp02_bird07"],
       t$expCorPower[t$exp_bird=="exp02_bird07"],
       pch=19,col=threecolz[2])
lines(t$pulses[t$exp_bird=="exp02_bird07"],
      t$expCorPower[t$exp_bird=="exp02_bird07"],
      lty=linezz[3],col=threecolz[2])
points(t$pulses[t$exp_bird=="exp03_bird07"],
       t$expCorPower[t$exp_bird=="exp03_bird07"],
       pch=19,col=threecolz[3])
lines(t$pulses[t$exp_bird=="exp03_bird07"],
      t$expCorPower[t$exp_bird=="exp03_bird07"],
      lty=linezz[3],col=threecolz[3])


#### fig 2 panel B ####
par(mfrow=c(1,7))
for (j in c(2,3,4,1,5,6,7)){
  plot(exp02_bird04_alltrials[[j]]$Position[[1]],
       exp02_bird04_alltrials[[j]]$Force[[1]],
       xlim=c(-1.5,2),ylim=c(0,4000),
       type='n',
       xlab="position",
       ylab="force",
       tck=0.02,bty='n',las=1)
  for (i in 1:3){
    lines(exp02_bird04_alltrials[[j]]$Position[[i]],
          exp02_bird04_alltrials[[j]]$Force[[i]],
          col=threecolz[2],lty=1,lwd=2)
  }
}
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1),mgp=c(3,1,0),las=0) # reset


#### fig 2 panel C ####
mult<-seq(1,1.3,0.05) # vector for placement of stimulations
mult<-(mult[c(4,1,2,3,5,6,7)])
plot(exp02_bird04_alltrials[[7]]$Time[[1]][-length(exp02_bird04_alltrials[[7]]$Time[[1]])],
     exp02_bird04_alltrials[[7]]$InstantPower[[1]]/bird04_pecmass,
     ylim=c(-1200,850),axes=FALSE,
     xlab="time (sec)",ylab="instantaneous power (W/kg)",
     type="n",tck=0.02,bty='n',las=1)
axis(2,at=c(-850,-425,0,425,850),tck=0.02,las=1)
axis(1,xlim=c(0.110,0.145),tck=0.02,las=1)
abline(h=0)
for (i in c(2,3,4,1,5,6,7)){
  lines(exp02_bird04_alltrials[[i]]$Time[[1]][-length(exp02_bird04_alltrials[[i]]$Time[[1]])],
        exp02_bird04_alltrials[[i]]$InstantPower[[1]]/bird04_pecmass,
        col=sevencolz[i],lwd=2)
  tmp<-filter(data.frame(time=exp02_bird04_alltrials[[i]]$Time[[1]],
                         stim=exp02_bird04_alltrials[[i]]$Stim[[1]]),stim>0)
  points(tmp$time,
         (-900*(mult[i]))*tmp$stim,
         pch=15,col=sevencolz[i])
}
rect((max(exp02_bird04_alltrials[[7]]$Time[[1]])+min(exp02_bird04_alltrials[[7]]$Time[[1]]))/2,-850,
     (max(exp02_bird04_alltrials[[7]]$Time[[1]])),850,
     col=rgb(0,0,0,alpha=0.05),border=NA)


# panels D, E, F use data from experiments 4-6
# panel D
t<-dplyr::filter(allsummarydat,
                 experiment=="exp04"|experiment=="exp05"|experiment=="exp06")
attach(t)
t<-t[order(phase,bird),]
detach(t)

# create color vectors
threecolz<-viridis(20,option="D")[c(1,10,15)]
#plot(1:3,c(rep(0,3)),pch=15,cex=5,col=threecolz)
tencolz<-viridis(10,option="A")[1:8]
#plot(1:8,c(rep(0,8)),pch=15,cex=5,col=tencolz)

# assign line types
# bird 2 is solid; bird 3 is dashed; bird 4 is dotted; bird 6 is dot-dashed
linezzz<-c(1:4) 

#### fig 2 panel D ####
plot(t$phase[t$exp_bird=="exp04_bird02"],
     t$expCorPower[t$exp_bird=="exp04_bird02"],
     xlim=c(-40,10),ylim=c(-200,150),
     pch=19,col=threecolz[1],
     xlab="onset (% cycle before peak length)",
     ylab="net power (W/kg)",
     tck=0.02,bty='n',las=1)
abline(h=0)
lines(t$phase[t$exp_bird=="exp04_bird02"],
      t$expCorPower[t$exp_bird=="exp04_bird02"],
      lty=linezzz[1],col=threecolz[1])
points(t$phase[t$exp_bird=="exp05_bird02"],
       t$expCorPower[t$exp_bird=="exp05_bird02"],
       pch=19,col=threecolz[2])
lines(t$phase[t$exp_bird=="exp05_bird02"],
      t$expCorPower[t$exp_bird=="exp05_bird02"],
      lty=linezzz[1],col=threecolz[2])
points(t$phase[t$exp_bird=="exp06_bird02"],
       t$expCorPower[t$exp_bird=="exp06_bird02"],
       pch=19,col=threecolz[3])
lines(t$phase[t$exp_bird=="exp06_bird02"],
      t$expCorPower[t$exp_bird=="exp06_bird02"],
      lty=linezzz[1],col=threecolz[3])
points(t$phase[t$exp_bird=="exp04_bird03"],
       t$expCorPower[t$exp_bird=="exp04_bird03"],
       pch=19,col=threecolz[1])
lines(t$phase[t$exp_bird=="exp04_bird03"],
      t$expCorPower[t$exp_bird=="exp04_bird03"],
      lty=linezzz[2],col=threecolz[1])
points(t$phase[t$exp_bird=="exp05_bird03"],
       t$expCorPower[t$exp_bird=="exp05_bird03"],
       pch=19,col=threecolz[2])
lines(t$phase[t$exp_bird=="exp05_bird03"],
      t$expCorPower[t$exp_bird=="exp05_bird03"],
      lty=linezzz[2],col=threecolz[2])
points(t$phase[t$exp_bird=="exp06_bird03"],
       t$expCorPower[t$exp_bird=="exp06_bird03"],
       pch=19,col=threecolz[3])
lines(t$phase[t$exp_bird=="exp06_bird03"],
      t$expCorPower[t$exp_bird=="exp06_bird03"],
      lty=linezzz[2],col=threecolz[3])
points(t$phase[t$exp_bird=="exp04_bird04"],
       t$expCorPower[t$exp_bird=="exp04_bird04"],
       pch=19,col=threecolz[1])
lines(t$phase[t$exp_bird=="exp04_bird04"],
      t$expCorPower[t$exp_bird=="exp04_bird04"],
      lty=linezzz[3],col=threecolz[1])
points(t$phase[t$exp_bird=="exp05_bird04"],
       t$expCorPower[t$exp_bird=="exp05_bird04"],
       pch=19,col=threecolz[2])
lines(t$phase[t$exp_bird=="exp05_bird04"],
      t$expCorPower[t$exp_bird=="exp05_bird04"],
      lty=linezzz[3],col=threecolz[2])
points(t$phase[t$exp_bird=="exp06_bird04"],
       t$expCorPower[t$exp_bird=="exp06_bird04"],
       pch=19,col=threecolz[3])
lines(t$phase[t$exp_bird=="exp06_bird04"],
      t$expCorPower[t$exp_bird=="exp06_bird04"],
      lty=linezzz[3],col=threecolz[3])
points(t$phase[t$exp_bird=="exp04_bird06"],
       t$expCorPower[t$exp_bird=="exp04_bird06"],
       pch=19,col=threecolz[1])
lines(t$phase[t$exp_bird=="exp04_bird06"],
      t$expCorPower[t$exp_bird=="exp04_bird06"],
      lty=linezzz[4],col=threecolz[1])
points(t$phase[t$exp_bird=="exp05_bird06"],
       t$expCorPower[t$exp_bird=="exp05_bird06"],
       pch=19,col=threecolz[2])
lines(t$phase[t$exp_bird=="exp05_bird06"],
      t$expCorPower[t$exp_bird=="exp05_bird06"],
      lty=linezzz[4],col=threecolz[2])
points(t$phase[t$exp_bird=="exp06_bird06"],
       t$expCorPower[t$exp_bird=="exp06_bird06"],
       pch=19,col=threecolz[3])
lines(t$phase[t$exp_bird=="exp06_bird06"],
      t$expCorPower[t$exp_bird=="exp06_bird06"],
      lty=linezzz[4],col=threecolz[3])


#### fig 2 panel E ####
par(mfrow=c(1,8))
for (j in c(9,4,8,1,7,3,6,2)){
  plot(exp05_bird03_alltrials[[j]]$Position[[1]],
       exp05_bird03_alltrials[[j]]$Force[[1]],
       xlim=c(-1.5,2),ylim=c(0,7000),
       type='n',
       xlab="position",
       ylab="force",
       tck=0.02,bty='n',las=1)
  for (i in 1:3){
    lines(exp05_bird03_alltrials[[j]]$Position[[i]],
          exp05_bird03_alltrials[[j]]$Force[[i]],
          col=threecolz[2],lty=1,lwd=2)
  }
}
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1),mgp=c(3,1,0),las=0) # reset



#### fig 2 panel F ####
mult<-rev(c(seq(1,1.2,0.05),seq(1.2,1.35,0.05))) # vector for placement of stimulations
plot(exp05_bird02_alltrials[[7]]$Time[[3]][-length(exp05_bird02_alltrials[[7]]$Time[[3]])],
     exp05_bird02_alltrials[[7]]$InstantPower[[3]]/bird02_pecmass,
     ylim=c(-1050,700),axes=FALSE,
     xlab="time (sec)",ylab="instantaneous power (W/kg)",
     type="n",tck=0.02,bty='n',las=1)
axis(2,at=c(-700,-350,0,350,700),tck=0.02,las=1)
axis(1,xlim=c(0.110,0.145),tck=0.02,las=1)
abline(h=0)
lines(exp05_bird02_alltrials[[9]]$Time[[3]][-length(exp05_bird02_alltrials[[9]]$Time[[3]])],
      exp05_bird02_alltrials[[9]]$InstantPower[[3]]/bird02_pecmass,
      col=tencolz[1],lwd=2)
tmp<-filter(data.frame(time=exp05_bird02_alltrials[[9]]$Time[[3]],
                       stim=exp05_bird02_alltrials[[9]]$Stim[[3]]),stim>0)
points(tmp$time,(-750*(mult[9]))*tmp$stim,pch=15,col=tencolz[1])
lines(exp05_bird02_alltrials[[4]]$Time[[3]][-length(exp05_bird02_alltrials[[4]]$Time[[3]])],
      exp05_bird02_alltrials[[4]]$InstantPower[[3]]/bird02_pecmass,
      col=tencolz[2],lwd=2)
tmp<-filter(data.frame(time=exp05_bird02_alltrials[[4]]$Time[[3]],
                       stim=exp05_bird02_alltrials[[4]]$Stim[[3]]),stim>0)
points(tmp$time,(-750*(mult[8]))*tmp$stim,pch=15,col=tencolz[2])
lines(exp05_bird02_alltrials[[8]]$Time[[3]][-length(exp05_bird02_alltrials[[8]]$Time[[3]])],
      exp05_bird02_alltrials[[8]]$InstantPower[[3]]/bird02_pecmass,
      col=tencolz[3],lwd=2)
tmp<-filter(data.frame(time=exp05_bird02_alltrials[[8]]$Time[[3]],
                       stim=exp05_bird02_alltrials[[8]]$Stim[[3]]),stim>0)
points(tmp$time,(-750*(mult[7]))*tmp$stim,pch=15,col=tencolz[3])
lines(exp05_bird02_alltrials[[1]]$Time[[3]][-length(exp05_bird02_alltrials[[1]]$Time[[3]])],
      exp05_bird02_alltrials[[1]]$InstantPower[[3]]/bird02_pecmass,
      col=tencolz[4],lwd=2)
tmp<-filter(data.frame(time=exp05_bird02_alltrials[[1]]$Time[[3]],
                       stim=exp05_bird02_alltrials[[1]]$Stim[[3]]),stim>0)
points(tmp$time,(-750*(mult[6]))*tmp$stim,pch=15,col=tencolz[4])
lines(exp05_bird02_alltrials[[7]]$Time[[3]][-length(exp05_bird02_alltrials[[7]]$Time[[3]])],
      exp05_bird02_alltrials[[7]]$InstantPower[[3]]/bird02_pecmass,
      col=tencolz[5],lwd=2)
tmp<-filter(data.frame(time=exp05_bird02_alltrials[[7]]$Time[[3]],
                       stim=exp05_bird02_alltrials[[7]]$Stim[[3]]),stim>0)
points(tmp$time,(-750*(mult[4]))*tmp$stim,pch=15,col=tencolz[5])
lines(exp05_bird02_alltrials[[3]]$Time[[3]][-length(exp05_bird02_alltrials[[3]]$Time[[3]])],
      exp05_bird02_alltrials[[3]]$InstantPower[[3]]/bird02_pecmass,
      col=tencolz[6],lwd=2)
tmp<-filter(data.frame(time=exp05_bird02_alltrials[[3]]$Time[[3]],
                       stim=exp05_bird02_alltrials[[3]]$Stim[[3]]),stim>0)
points(tmp$time,(-750*(mult[3]))*tmp$stim,pch=15,col=tencolz[6])
lines(exp05_bird02_alltrials[[6]]$Time[[3]][-length(exp05_bird02_alltrials[[6]]$Time[[3]])],
      exp05_bird02_alltrials[[6]]$InstantPower[[3]]/bird02_pecmass,
      col=tencolz[7],lwd=2)
tmp<-filter(data.frame(time=exp05_bird02_alltrials[[6]]$Time[[3]],
                       stim=exp05_bird02_alltrials[[6]]$Stim[[3]]),stim>0)
points(tmp$time,(-750*(mult[2]))*tmp$stim,pch=15,col=tencolz[7])
lines(exp05_bird02_alltrials[[2]]$Time[[3]][-length(exp05_bird02_alltrials[[2]]$Time[[3]])],
      exp05_bird02_alltrials[[2]]$InstantPower[[3]]/bird02_pecmass,
      col=tencolz[8],lwd=2)
tmp<-filter(data.frame(time=exp05_bird02_alltrials[[2]]$Time[[3]],
                       stim=exp05_bird02_alltrials[[2]]$Stim[[3]]),stim>0)
points(tmp$time,(-750*(mult[1]))*tmp$stim,pch=15,col=tencolz[8])
rect((max(exp05_bird02_alltrials[[7]]$Time[[3]])+min(exp05_bird02_alltrials[[7]]$Time[[3]]))/2,-700,
     (max(exp05_bird02_alltrials[[7]]$Time[[3]])),700,
     col=rgb(0,0,0,alpha=0.05),border=NA)



################################### figure 3 ###################################
library(ggplot2);library(ggridges)
library(tidybayes);library(tidyr);library(forcats)

# set up colors and data first
# panels A, B, C use data from experiments 13-14
cols_abcramp<-colorRampPalette(c("#000000","#A020F0"))
cols_abcdos<-c("#662D91","#A020F0")
lines_abc<-1:4
ca<-dplyr::filter(allsummarydat,experiment=="exp13") # panel C, amplitude
cf<-dplyr::filter(allsummarydat,experiment=="exp14") # panel C, frequency

# panels D, E, F use data from experiments 11-12
cols_deframp<-colorRampPalette(c("#000000","red"))
cols_defdos<-c("#C1272D","#FF3333")
lines_def<-1:4
fa<-dplyr::filter(allsummarydat,experiment=="exp11") # panel F, amplitude
ff<-dplyr::filter(allsummarydat,experiment=="exp12") # panel F, frequency

# panels G, H, I use data from experiments 9-10
cols_ghiramp<-colorRampPalette(c("#000000","green"))
cols_ghidos<-c("#009245","#8CC63F")
lines_ghi<-1:4
ha<-dplyr::filter(allsummarydat,experiment=="exp09") # panel I, amplitude
hf<-dplyr::filter(allsummarydat,experiment=="exp10") # panel I, frequency
# exclusions:
# exclude all 15.5% trials b/c stimulation was accidentally carried out longer
ha<-ha[-(as.numeric(row.names(ha[ha$amplitude==3.267,]))),]

# panels J, K, L use data from experiments 7-8
cols_jklramp<-colorRampPalette(c("#000000","blue"))
cols_jkldos<-c("#2A51E0","#85CFEC")
lines_jkl<-1:4
la<-dplyr::filter(allsummarydat,experiment=="exp07") # panel L, amplitude
lf<-dplyr::filter(allsummarydat,experiment=="exp08") # panel L, frequency


#### fig 3 panel A ####
tidyr::gather(data.frame(fit_freqinteract$Sol[,1:6]))->df
df%>%mutate(var=fct_reorder(key,value,.fun=median))->m
ggplot(m,aes(x=value,y=var))+
  geom_halfeyeh(.width=0.95) +
  coord_cartesian(xlim=c(-1,1)) +
  geom_vline(xintercept=0) + theme_ridges()


#### fig 3 panel C ####
plot(ca$shortening_velocity,ca$expCorPower,
     xlim=c(150,200),ylim=c(0,150),
     type='n',
     xlab="mean muscle shortening velocity (mm/sec)",
     ylab="net power (W/kg)",
     tck=0.02,bty='n',las=1)
fig3c_birdordera<-c("bird02","bird03","bird05","bird08")
fig3c_birdorderf<-c("bird09","bird10","bird11","bird12")
for (i in 1:4){
  tmp<-ca[ca$bird==fig3c_birdordera[i],]
  tmp<-tmp[order(tmp$shortening_velocity),]
  points(tmp$shortening_velocity,tmp$expCorPower,
         pch=19,col=cols_abcdos[1])
  lines(tmp$shortening_velocity,tmp$expCorPower,
        col=cols_abcdos[1],lty=lines_abc[i],lwd=2)
}
for (i in 1:4){
  tmp<-cf[cf$bird==fig3c_birdorderf[i],]
  tmp<-tmp[order(tmp$shortening_velocity),]
  points(tmp$shortening_velocity,tmp$expCorPower,
         pch=1,col=cols_abcdos[2])
  lines(tmp$shortening_velocity,tmp$expCorPower,
        col=cols_abcdos[2],lty=lines_abc[i],lwd=2)
}


#### fig 3 panel D ####
plot(fa$shortening_velocity,fa$expCorPower,
     xlim=c(150,200),ylim=c(0,150),
     type='n',
     xlab="mean muscle shortening velocity (mm/sec)",
     ylab="net power (W/kg)",
     tck=0.02,bty='n',las=1)
fig3f_boa<-c("bird02","bird03","bird05","bird08")
fig3f_bof<-c("bird02","bird03","bird05","bird10")
for (i in 1:4){
  tmp<-fa[fa$bird==fig3f_boa[i],]
  tmp<-tmp[order(tmp$shortening_velocity),]
  points(tmp$shortening_velocity,tmp$expCorPower,
         pch=19,col=cols_defdos[1])
  lines(tmp$shortening_velocity,tmp$expCorPower,
        col=cols_defdos[1],lty=lines_def[i],lwd=2)
}
for (i in 1:4){
  tmp<-ff[ff$bird==fig3f_bof[i],]
  tmp<-tmp[order(tmp$shortening_velocity),]
  points(tmp$shortening_velocity,tmp$expCorPower,
         pch=1,col=cols_defdos[2])
  lines(tmp$shortening_velocity,tmp$expCorPower,
        col=cols_defdos[2],lty=lines_def[i],lwd=2)
}


#### fig 3 panel E ####
plot(ha$shortening_velocity,ha$expCorPower,
     xlim=c(150,200),ylim=c(0,150),
     type='n',
     xlab="mean muscle shortening velocity (mm/sec)",
     ylab="net power (W/kg)",
     tck=0.02,bty='n',las=1)
fig3i_boa<-c("bird03","bird05","bird08","bird12")
fig3i_bof<-c("bird09","bird10","bird11","bird12")
for (i in 1:4){
  tmp<-ha[ha$bird==fig3i_boa[i],]
  tmp<-tmp[order(tmp$shortening_velocity),]
  points(tmp$shortening_velocity,tmp$expCorPower,
         pch=19,col=cols_ghidos[1])
  lines(tmp$shortening_velocity,tmp$expCorPower,
        col=cols_ghidos[1],lty=lines_ghi[i],lwd=2)
}
for (i in 1:4){
  tmp<-hf[hf$bird==fig3i_bof[i],]
  tmp<-tmp[order(tmp$shortening_velocity),]
  points(tmp$shortening_velocity,tmp$expCorPower,
         pch=1,col=cols_ghidos[2])
  lines(tmp$shortening_velocity,tmp$expCorPower,
        col=cols_ghidos[2],lty=lines_ghi[i],lwd=2)
}


#### fig 3 panel F ####
plot(la$shortening_velocity,la$expCorPower,
     xlim=c(150,200),ylim=c(0,150),
     type='n',
     xlab="mean muscle shortening velocity (mm/sec)",
     ylab="net power (W/kg)",
     tck=0.02,bty='n',las=1)
fig3l_boa<-c("bird02","bird03","bird05","bird08")
fig3l_bof<-c("bird02","bird03","bird08","bird12") 
for (i in 1:4){
  tmp<-la[la$bird==fig3l_boa[i],]
  tmp<-tmp[order(tmp$shortening_velocity),]
  points(tmp$shortening_velocity,tmp$expCorPower,
         pch=19,col=cols_jkldos[1])
  lines(tmp$shortening_velocity,tmp$expCorPower,
        col=cols_jkldos[1],lty=lines_jkl[i],lwd=2)
}
for (i in 1:4){
  tmp<-lf[lf$bird==fig3l_bof[i],]
  tmp<-tmp[order(tmp$shortening_velocity),]
  points(tmp$shortening_velocity,tmp$expCorPower,
         pch=1,col=cols_jkldos[2])
  lines(tmp$shortening_velocity,tmp$expCorPower,
        col=cols_jkldos[2],lty=lines_jkl[i],lwd=2)
}


#### fig 3 panels C, D, E, F together #### 
par(mfrow=c(1,4))
plot(ca$shortening_velocity,ca$expCorPower,
     xlim=c(150,200),ylim=c(0,150),
     type='n',
     xlab="mean muscle shortening velocity (mm/sec)",
     ylab="net power (W/kg)",
     tck=0.02,bty='n',las=1)
fig3c_birdordera<-c("bird02","bird03","bird05","bird08")
fig3c_birdorderf<-c("bird09","bird10","bird11","bird12")
for (i in 1:4){
  tmp<-ca[ca$bird==fig3c_birdordera[i],]
  tmp<-tmp[order(tmp$shortening_velocity),]
  points(tmp$shortening_velocity,tmp$expCorPower,
         pch=19,col=cols_abcdos[1])
  lines(tmp$shortening_velocity,tmp$expCorPower,
        col=cols_abcdos[1],lty=lines_abc[i],lwd=2)
}
for (i in 1:4){
  tmp<-cf[cf$bird==fig3c_birdorderf[i],]
  tmp<-tmp[order(tmp$shortening_velocity),]
  points(tmp$shortening_velocity,tmp$expCorPower,
         pch=1,col=cols_abcdos[2])
  lines(tmp$shortening_velocity,tmp$expCorPower,
        col=cols_abcdos[2],lty=lines_abc[i],lwd=2)
}
plot(fa$shortening_velocity,fa$expCorPower,
     xlim=c(150,200),ylim=c(0,150),
     type='n',
     xlab="mean muscle shortening velocity (mm/sec)",
     ylab="net power (W/kg)",
     tck=0.02,bty='n',las=1)
fig3f_boa<-c("bird02","bird03","bird05","bird08")
fig3f_bof<-c("bird02","bird03","bird05","bird10")
for (i in 1:4){
  tmp<-fa[fa$bird==fig3f_boa[i],]
  tmp<-tmp[order(tmp$shortening_velocity),]
  points(tmp$shortening_velocity,tmp$expCorPower,
         pch=19,col=cols_defdos[1])
  lines(tmp$shortening_velocity,tmp$expCorPower,
        col=cols_defdos[1],lty=lines_def[i],lwd=2)
}
for (i in 1:4){
  tmp<-ff[ff$bird==fig3f_bof[i],]
  tmp<-tmp[order(tmp$shortening_velocity),]
  points(tmp$shortening_velocity,tmp$expCorPower,
         pch=1,col=cols_defdos[2])
  lines(tmp$shortening_velocity,tmp$expCorPower,
        col=cols_defdos[2],lty=lines_def[i],lwd=2)
}
plot(ha$shortening_velocity,ha$expCorPower,
     xlim=c(150,200),ylim=c(0,150),
     type='n',
     xlab="mean muscle shortening velocity (mm/sec)",
     ylab="net power (W/kg)",
     tck=0.02,bty='n',las=1)
fig3i_boa<-c("bird03","bird05","bird08","bird12")
fig3i_bof<-c("bird09","bird10","bird11","bird12")
for (i in 1:4){
  tmp<-ha[ha$bird==fig3i_boa[i],]
  tmp<-tmp[order(tmp$shortening_velocity),]
  points(tmp$shortening_velocity,tmp$expCorPower,
         pch=19,col=cols_ghidos[1])
  lines(tmp$shortening_velocity,tmp$expCorPower,
        col=cols_ghidos[1],lty=lines_ghi[i],lwd=2)
}
for (i in 1:4){
  tmp<-hf[hf$bird==fig3i_bof[i],]
  tmp<-tmp[order(tmp$shortening_velocity),]
  points(tmp$shortening_velocity,tmp$expCorPower,
         pch=1,col=cols_ghidos[2])
  lines(tmp$shortening_velocity,tmp$expCorPower,
        col=cols_ghidos[2],lty=lines_ghi[i],lwd=2)
}
plot(la$shortening_velocity,la$expCorPower,
     xlim=c(150,200),ylim=c(0,150),
     type='n',
     xlab="mean muscle shortening velocity (mm/sec)",
     ylab="net power (W/kg)",
     tck=0.02,bty='n',las=1)
fig3l_boa<-c("bird02","bird03","bird05","bird08")
fig3l_bof<-c("bird02","bird03","bird08","bird12")
for (i in 1:4){
  tmp<-la[la$bird==fig3l_boa[i],]
  tmp<-tmp[order(tmp$shortening_velocity),]
  points(tmp$shortening_velocity,tmp$expCorPower,
         pch=19,col=cols_jkldos[1])
  lines(tmp$shortening_velocity,tmp$expCorPower,
        col=cols_jkldos[1],lty=lines_jkl[i],lwd=2)
}
for (i in 1:4){
  tmp<-lf[lf$bird==fig3l_bof[i],]
  tmp<-tmp[order(tmp$shortening_velocity),]
  points(tmp$shortening_velocity,tmp$expCorPower,
         pch=1,col=cols_jkldos[2])
  lines(tmp$shortening_velocity,tmp$expCorPower,
        col=cols_jkldos[2],lty=lines_jkl[i],lwd=2)
}
par(mfrow=c(1,1))


#### fig 3 panel G ####
plot(seq(0,1,1/(length(exp13_bird02_alltrials[[7]]$InstantPower[[3]])-1)),
     exp13_bird02_alltrials[[7]]$InstantPower[[3]]/bird02_pecmass,
     ylim=c(-450,900),axes=FALSE,
     xlab="percent cycle",ylab="instantaneous power (W/kg)",
     type="n",tck=0.02,bty='n',las=1)
axis(2,at=c(-450,0,450,900),tck=0.02,las=1)
axis(1,at=c(0,0.5,1),tck=0.02,las=1)
abline(h=0)
fig3a_to<-c(2,4,5,7,1,8,3)
for (i in 1:7){
  lines(seq(0,1,1/(length(exp13_bird02_alltrials[[fig3a_to[i]]]$InstantPower[[3]])-1)),
        exp13_bird02_alltrials[[fig3a_to[i]]]$InstantPower[[3]]/bird02_pecmass,
        col=cols_abcramp(7)[i],lwd=2)
}
rect(0.5,-450,
     1,900,
     col=rgb(0,0,0,alpha=0.05),border=NA)


#### fig 3 panel H ####
plot(seq(0,1,1/(length(exp11_bird02_alltrials[[7]]$InstantPower[[3]])-1)),
     exp11_bird02_alltrials[[7]]$InstantPower[[3]]/bird02_pecmass,
     ylim=c(-450,900),axes=FALSE,
     xlab="percent cycle",ylab="instantaneous power (W/kg)",
     type="n",tck=0.02,bty='n',las=1)
axis(2,at=c(-450,0,450,900),tck=0.02,las=1)
axis(1,at=c(0,0.5,1),tck=0.02,las=1)
abline(h=0)
fig3d_to<-c(2,4,5,7,1,6,3)
for (i in 1:7){
  lines(seq(0,1,1/(length(exp11_bird02_alltrials[[fig3d_to[i]]]$InstantPower[[3]])-1)),
        exp11_bird02_alltrials[[fig3d_to[i]]]$InstantPower[[3]]/bird02_pecmass,
        col=cols_deframp(7)[i],lwd=2)
}
rect(0.5,-450,
     1,900,
     col=rgb(0,0,0,alpha=0.05),border=NA)


#### fig 3 panel I ####
plot(seq(0,1,1/(length(exp09_bird03_alltrials[[7]]$InstantPower[[3]])-1)),
     exp09_bird03_alltrials[[7]]$InstantPower[[3]]/bird03_pecmass,
     ylim=c(-450,900),axes=FALSE,
     xlab="percent cycle",ylab="instantaneous power (W/kg)",
     type="n",tck=0.02,bty='n',las=1)
axis(2,at=c(-450,0,450,900),tck=0.02,las=1)
axis(1,at=c(0,0.5,1),tck=0.02,las=1)
abline(h=0)
fig3g_to<-c(2,4,5,7,1,8,3)
for (i in 1:7){
  lines(seq(0,1,1/(length(exp09_bird03_alltrials[[fig3g_to[i]]]$InstantPower[[3]])-1)),
        exp09_bird03_alltrials[[fig3g_to[i]]]$InstantPower[[3]]/bird03_pecmass,
        col=cols_ghiramp(7)[i],lwd=2)
}
rect(0.5,-450,
     1,900,
     col=rgb(0,0,0,alpha=0.05),border=NA)


#### fig 4 panel J ####
plot(seq(0,1,1/(length(exp07_bird02_alltrials[[7]]$InstantPower[[3]])-1)),
     exp07_bird02_alltrials[[7]]$InstantPower[[3]]/bird02_pecmass,
     ylim=c(-450,900),axes=FALSE,
     xlab="percent cycle",ylab="instantaneous power (W/kg)",
     type="n",tck=0.02,bty='n',las=1)
axis(2,at=c(-450,0,450,900),tck=0.02,las=1)
axis(1,at=c(0,0.5,1),tck=0.02,las=1)
abline(h=0)
fig3j_to<-c(2,4,5,7,1,6,3)
for (i in 1:7){
  lines(seq(0,1,1/(length(exp07_bird02_alltrials[[fig3j_to[i]]]$InstantPower[[3]])-1)),
        exp07_bird02_alltrials[[fig3j_to[i]]]$InstantPower[[3]]/bird02_pecmass,
        col=cols_jklramp(7)[i],lwd=2)
}
rect(0.5,-450,
     1,900,
     col=rgb(0,0,0,alpha=0.05),border=NA)


#### fig 3 panel K ####
plot(seq(0,1,1/(length(exp14_bird09_alltrials[[7]]$InstantPower[[3]])-1)),
     exp14_bird09_alltrials[[7]]$InstantPower[[3]]/bird09_pecmass,
     ylim=c(-450,900),axes=FALSE,
     xlab="percent cycle",ylab="instantaneous power (W/kg)",
     type="n",tck=0.02,bty='n',las=1)
axis(2,at=c(-450,0,450,900),tck=0.02,las=1)
axis(1,at=c(0,0.5,1),tck=0.02,las=1)
abline(h=0)
fig3b_to<-c(2,6,1,7,5,3)
for (i in 1:6){
  lines(seq(0,1,1/(length(exp14_bird09_alltrials[[fig3b_to[i]]]$InstantPower[[3]])-1)),
        exp14_bird09_alltrials[[fig3b_to[i]]]$InstantPower[[3]]/bird09_pecmass,
        col=cols_abcramp(6)[i],lwd=2)
}
rect(0.5,-450,
     1,900,
     col=rgb(0,0,0,alpha=0.05),border=NA)


#### fig 4 panel L ####
plot(seq(0,1,1/(length(exp12_bird02_alltrials[[7]]$InstantPower[[3]])-1)),
     exp12_bird02_alltrials[[7]]$InstantPower[[3]]/bird09_pecmass,
     ylim=c(-450,900),axes=FALSE,
     xlab="percent cycle",ylab="instantaneous power (W/kg)",
     type="n",tck=0.02,bty='n',las=1)
axis(2,at=c(-450,0,450,900),tck=0.02,las=1)
axis(1,at=c(0,0.5,1),tck=0.02,las=1)
abline(h=0)
fig3e_to<-c(2,4,6,1,7,5,3)
for (i in 1:7){
  lines(seq(0,1,1/(length(exp12_bird02_alltrials[[fig3e_to[i]]]$InstantPower[[3]])-1)),
        exp12_bird02_alltrials[[fig3e_to[i]]]$InstantPower[[3]]/bird02_pecmass,
        col=cols_deframp(7)[i],lwd=2)
}
rect(0.5,-450,
     1,900,
     col=rgb(0,0,0,alpha=0.05),border=NA)


#### fig 3 panel M ####
plot(seq(0,1,1/(length(exp10_bird09_alltrials[[7]]$InstantPower[[3]])-1)),
     exp10_bird09_alltrials[[7]]$InstantPower[[3]]/bird09_pecmass,
     ylim=c(-450,900),axes=FALSE,
     xlab="percent cycle",ylab="instantaneous power (W/kg)",
     type="n",tck=0.02,bty='n',las=1)
axis(2,at=c(-450,0,450,900),tck=0.02,las=1)
axis(1,at=c(0,0.5,1),tck=0.02,las=1)
abline(h=0)
fig3h_to<-c(2,4,6,1,7,5,3)
for (i in 1:7){
  lines(seq(0,1,1/(length(exp10_bird09_alltrials[[fig3h_to[i]]]$InstantPower[[3]])-1)),
        exp10_bird09_alltrials[[fig3h_to[i]]]$InstantPower[[3]]/bird09_pecmass,
        col=cols_ghiramp(7)[i],lwd=2)
}
rect(0.5,-450,
     1,900,
     col=rgb(0,0,0,alpha=0.05),border=NA)


#### fig 3 panel N ####
plot(seq(0,1,1/(length(exp08_bird02_alltrials[[7]]$InstantPower[[3]])-1)),
     exp08_bird02_alltrials[[7]]$InstantPower[[3]]/bird02_pecmass,
     ylim=c(-450,900),axes=FALSE,
     xlab="percent cycle",ylab="instantaneous power (W/kg)",
     type="n",tck=0.02,bty='n',las=1)
axis(2,at=c(-450,0,450,900),tck=0.02,las=1)
axis(1,at=c(0,0.5,1),tck=0.02,las=1)
abline(h=0)
fig3k_to<-c(2,4,6,1,7,5,3)
for (i in 1:7){
  lines(seq(0,1,1/(length(exp08_bird02_alltrials[[fig3k_to[i]]]$InstantPower[[3]])-1)),
        exp08_bird02_alltrials[[fig3k_to[i]]]$InstantPower[[3]]/bird02_pecmass,
        col=cols_jklramp(7)[i],lwd=2)
}
rect(0.5,-450,
     1,900,
     col=rgb(0,0,0,alpha=0.05),border=NA)


################################### figure 4 ###################################

###### fig 4 panel A ######

## Work loops of 4 pulse -20 phase (in vivo)
## Amplitude
par(mfrow=c(1,8))
for (j in c(2,4,5,7,1,6,8,3)){
  plot(exp13_bird02_alltrials[[j]]$Position[[3]],
       exp13_bird02_alltrials[[j]]$Force[[3]],
       xlim=c(0,4),ylim=c(0,4000),
       type='n',
       xlab="position",
       ylab="force (mN)",
       tck=0.02,bty='n',las=1)
  for (i in 3:3){
    lines(exp13_bird02_alltrials[[j]]$Position[[i]],
          exp13_bird02_alltrials[[j]]$Force[[i]],
          col="#662D91",lty=1,lwd=2)
  }
}
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1),mgp=c(3,1,0),las=0) # reset
## Frequency
par(mfrow=c(1,8))
for (j in c(2,6,1,7,5,3)){
  plot(exp14_bird09_alltrials[[j]]$Position[[3]]-1,
       exp14_bird09_alltrials[[j]]$Force[[3]],
       xlim=c(0,4),ylim=c(0,4000),
       type='n',
       xlab="position",
       ylab="force (mN)",
       tck=0.02,bty='n',las=1)
  for (i in 3:3){
    lines(exp14_bird09_alltrials[[j]]$Position[[i]]-1,
          exp14_bird09_alltrials[[j]]$Force[[i]],
          col="#662D91",lty=1,lwd=2)
  }
}
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1),mgp=c(3,1,0),las=0) # reset

###### fig 4 panel B ######
## Work loops of 3 pulse -30 phase (least extreme)
## Amplitude
par(mfrow=c(1,8))
for (j in c(2,4,5,7,1,6,8,3)){
  plot(exp07_bird03_alltrials[[j]]$Position[[3]],
       exp07_bird03_alltrials[[j]]$Force[[3]],
       xlim=c(-2,2),ylim=c(0,6000),
       type='n',
       xlab="position",
       ylab="force",
       tck=0.02,bty='n',las=1)
  for (i in 3:3){
    lines(exp07_bird03_alltrials[[j]]$Position[[i]],
          exp07_bird03_alltrials[[j]]$Force[[i]],
          col="#5C9BDC",lty=1,lwd=2)
  }
}
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1),mgp=c(3,1,0),las=0) # reset
## Frequency
par(mfrow=c(1,8))
for (j in c(2,4,6,1,7,5,3)){
  plot(exp08_bird03_alltrials[[j]]$Position[[3]],
       exp08_bird03_alltrials[[j]]$Force[[3]],
       xlim=c(-2,2),ylim=c(0,6000),
       type='n',
       xlab="position",
       ylab="force",
       tck=0.02,bty='n',las=1)
  for (i in 3:3){
    lines(exp08_bird03_alltrials[[j]]$Position[[i]],
          exp08_bird03_alltrials[[j]]$Force[[i]],
          col="#5C9BDC",lty=1,lwd=2)
  }
}
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1),mgp=c(3,1,0),las=0) # reset



###### fig 4 panel C & D ######
## Time course of 4 pulse -20 phase (in vivo)
## Amplitude
par(mfrow=c(1,8))
for (j in c(2,4,5,7,1,6,8,3)){
  plot(seq(0,1,1/(length(exp13_bird02_alltrials[[j]]$Time[[3]])-1)),
       exp13_bird02_alltrials[[j]]$Force[[3]],
       ylim=c(0,4000),axes=FALSE,
       type='n',
       xlab="percent cycle",
       ylab="force",
       tck=0.02,bty='n',las=1)
  axis(2,at=c(0,2000,4000),tck=0.02,las=1)
  axis(1,at=c(0,0.5,1),tck=0.02,las=1)
  for (i in 3:3){
    lines(seq(0,1,1/(length(exp13_bird02_alltrials[[j]]$Time[[3]])-1)),
          exp13_bird02_alltrials[[j]]$Force[[i]],
          col="#662D91",lty=1,lwd=2)
  }
}
par(mfrow=c(1,8))
for (j in c(2,4,5,7,1,6,8,3)){
  plot(seq(0,1,1/(length(exp13_bird02_alltrials[[j]]$InstantPower[[3]])-1)),
       exp13_bird02_alltrials[[j]]$InstantPower[[3]]/bird02_pecmass,
       ylim=c(-450,900),axes=FALSE,
       xlab="percent cycle",
       ylab="instantaneous power (W/kg)",
       type="n",tck=0.02,bty='n',las=1)
  axis(2,at=c(-450,0,450,900),tck=0.02,las=1)
  axis(1,at=c(0,0.5,1),tck=0.02,las=1)
  for (i in 3:3){
    lines(seq(0,1,1/(length(exp13_bird02_alltrials[[j]]$InstantPower[[3]])-1)),
          exp13_bird02_alltrials[[j]]$InstantPower[[i]]/bird02_pecmass,
          col="black",lty=1,lwd=2)
  }
}
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1),mgp=c(3,1,0),las=0) # reset

## Frequency
par(mfrow=c(1,8))
for (j in c(2,6,1,7,5,3)){
  plot(seq(0,1,1/(length(exp14_bird09_alltrials[[j]]$Time[[3]])-1)),
       exp14_bird09_alltrials[[j]]$Force[[3]],
       ylim=c(0,4000),axes=FALSE,
       type='n',
       xlab="percent cycle",
       ylab="force",
       tck=0.02,bty='n',las=1)
  axis(2,at=c(0,2000,4000),tck=0.02,las=1)
  axis(1,at=c(0,0.5,1),tck=0.02,las=1)
  for (i in 3:3){
    lines(seq(0,1,1/(length(exp14_bird09_alltrials[[j]]$Time[[3]])-1)),
          exp14_bird09_alltrials[[j]]$Force[[i]],
          col="#662D91",lty=1,lwd=2)
  }
}
par(mfrow=c(1,8))
for (j in c(2,6,1,7,5,3)){
  plot(seq(0,1,1/(length(exp14_bird09_alltrials[[j]]$InstantPower[[3]])-1)),
       exp14_bird09_alltrials[[j]]$InstantPower[[3]]/bird09_pecmass,
       ylim=c(-450,900),axes=FALSE,
       xlab="percent cycle",
       ylab="instantaneous power (W/kg)",
       type="n",tck=0.02,bty='n',las=1)
  axis(2,at=c(-450,0,450,900),tck=0.02,las=1)
  axis(1,at=c(0,0.5,1),tck=0.02,las=1)
  for (i in 3:3){
    lines(seq(0,1,1/(length(exp14_bird09_alltrials[[j]]$InstantPower[[3]])-1)),
          exp14_bird09_alltrials[[j]]$InstantPower[[i]]/bird09_pecmass,
          col="black",lty=1,lwd=2)
  }
}
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1),mgp=c(3,1,0),las=0) # reset

## Time course of 3 pulse -30 phase (least extreme)
## Amplitude
par(mfrow=c(1,8))
for (j in c(2,4,5,7,1,6,8,3)){
  plot(seq(0,1,1/(length(exp07_bird02_alltrials[[j]]$Time[[3]])-1)),
       exp07_bird02_alltrials[[j]]$Force[[3]],
       ylim=c(0,4000),axes=FALSE,
       type='n',
       xlab="percent cycle",
       ylab="force",
       tck=0.02,bty='n',las=1)
  axis(2,at=c(0,2000,4000),tck=0.02,las=1)
  axis(1,at=c(0,0.5,1),tck=0.02,las=1)
  for (i in 3:3){
    lines(seq(0,1,1/(length(exp13_bird02_alltrials[[j]]$Time[[3]])-1)),
          exp13_bird02_alltrials[[j]]$Force[[i]],
          col="#662D91",lty=1,lwd=2)
  }
}
par(mfrow=c(1,8))
for (j in c(2,4,5,7,1,6,8,3)){
  plot(seq(0,1,1/(length(exp07_bird02_alltrials[[j]]$InstantPower[[3]])-1)),
       exp07_bird02_alltrials[[j]]$InstantPower[[3]]/bird02_pecmass,
       ylim=c(-450,900),axes=FALSE,
       xlab="percent cycle",
       ylab="instantaneous power (W/kg)",
       type="n",tck=0.02,bty='n',las=1)
  axis(2,at=c(-450,0,450,900),tck=0.02,las=1)
  axis(1,at=c(0,0.5,1),tck=0.02,las=1)
  for (i in 3:3){
    lines(seq(0,1,1/(length(exp13_bird02_alltrials[[j]]$InstantPower[[3]])-1)),
          exp13_bird02_alltrials[[j]]$InstantPower[[i]]/bird02_pecmass,
          col="black",lty=1,lwd=2)
  }
}
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1),mgp=c(3,1,0),las=0) # reset

## Frequency
par(mfrow=c(1,8))
for (j in c(2,4,6,1,7,5,3)){
  plot(seq(0,1,1/(length(exp08_bird02_alltrials[[j]]$Time[[3]])-1)),
       exp08_bird02_alltrials[[j]]$Force[[3]],
       ylim=c(0,4000),axes=FALSE,
       type='n',
       xlab="percent cycle",
       ylab="force",
       tck=0.02,bty='n',las=1)
  axis(2,at=c(0,2000,4000),tck=0.02,las=1)
  axis(1,at=c(0,0.5,1),tck=0.02,las=1)
  for (i in 3:3){
    lines(seq(0,1,1/(length(exp08_bird02_alltrials[[j]]$Time[[3]])-1)),
          exp08_bird02_alltrials[[j]]$Force[[i]],
          col="#662D91",lty=1,lwd=2)
  }
}
par(mfrow=c(1,8))
for (j in c(2,4,6,1,7,5,3)){
  plot(seq(0,1,1/(length(exp08_bird02_alltrials[[j]]$InstantPower[[3]])-1)),
       exp08_bird02_alltrials[[j]]$InstantPower[[3]]/bird02_pecmass,
       ylim=c(-450,900),axes=FALSE,
       xlab="percent cycle",
       ylab="instantaneous power (W/kg)",
       type="n",tck=0.02,bty='n',las=1)
  axis(2,at=c(-450,0,450,900),tck=0.02,las=1)
  axis(1,at=c(0,0.5,1),tck=0.02,las=1)
  for (i in 3:3){
    lines(seq(0,1,1/(length(exp08_bird02_alltrials[[j]]$InstantPower[[3]])-1)),
          exp08_bird02_alltrials[[j]]$InstantPower[[i]]/bird02_pecmass,
          col="black",lty=1,lwd=2)
  }
}
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1),mgp=c(3,1,0),las=0) # reset



################################### figure 5 ###################################
allsummarydat[allsummarydat$amplitude>0.79,] -> dat
colorz<-c("#999999","#FD9567","#721F81")
plot(dat$duty_cycle,dat$expCorPower,
     pch=19,tck=0.02,bty="n",xlim=c(0,1),
     col=c(rep(colorz[1],170),
           rep(colorz[2],35), #exp7
           rep(colorz[3],33), #exp8
           rep(colorz[2],35), #exp11
           rep(colorz[3],33), #exp12
           rep(colorz[2],36), #exp13
           rep(colorz[2],36), #exp9
           rep(colorz[3],35), #exp10
           rep(colorz[3],35)  #exp14
     ),
     ylab="net power (J)",xlab="duty cycle")


################################### figure 6 ###################################
## This figure was made entirely in Adobe Illustrator, based on data from 
## Figure 3 (see above)

