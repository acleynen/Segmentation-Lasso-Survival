library(plyr)
library(shape)
library(glmnet)
library(survival)

library(car)
library(ggplot2)
library(ggtext)
library(ggpubr)
library(survminer)
library(zoo)
library(pROC)
library(Rmisc)
library(pdfCluster)
library(Segmentor3IsBack)
library(dplyr)
library(parallel)
library(MCMCpack)


source("scripts/Methods.R")
source("scripts/Criteria.R")
source("scripts/RunMethods.R")
set.seed(1234)


#####################
load("Data/Test30000.RData")
load("Data/Train30000.RData")

N=30000
n=ncol(Xcov[[1]])
D=nrow(Regions)
RealR=rep(0,N)
BetaClass=rep(0,N)

for (i in 1:D)
{
	RealR[Regions$start[i]:(Regions$start[i]+Regions$length[i]-1)]<-1
	BetaClass[Regions$start[i]:(Regions$start[i]+Regions$length[i]-1)]<-i
}


RunAllMethods<-function(i)
{
	set.seed(i)
	print(i)

	XcovPat=Xcov[[i]]
	XcovPatTest=XcovTest[[i]]
	
	SurvPat=SurvivalPatients[[i]]
	SurvPatTest=SurvivalPatientsTest[[i]]
	
	PatientSelect<-sample(1:n,300,replace=F) # splits the data in a 300-Select and 200-Fit sets
	PatientFit<-(1:n)[!is.element(1:n,PatientSelect)]

	XcovSelect<-XcovPat[,PatientSelect]
	XcovFit<-XcovPat[,PatientFit]

	SurvivalSelect<-SurvPat[PatientSelect,]
	SurvivalFit<-SurvPat[PatientFit,]

	
	################################
	
	

	Gold=RunGold(XcovSelect,SurvivalSelect,XcovFit,SurvivalFit)
	Lasso=RunLasso(XcovSelect,SurvivalSelect,XcovFit,SurvivalFit)
	SegPval=RunCoxUni(XcovSelect,SurvivalSelect,XcovFit,SurvivalFit)
	AdjPv=RunMultipleTesting(XcovSelect,SurvivalSelect,XcovFit,SurvivalFit)
	
 
  ######################################"

	AUCGold=AllAUC(Gold$Ragreg,Gold$Surv,SurvPatTest,XcovPatTest)
	LseAUC=AllAUC(Lasso$RagregL2,Lasso$Surv_L2,SurvPatTest,XcovPatTest)
	LsminAUC=AllAUC(Lasso$RagregL1,Lasso$Surv_L1,SurvPatTest,XcovPatTest)
	PvalAUC=AllAUC(SegPval$Ragreg,SegPval$Surv,SurvPatTest,XcovPatTest)
	AUCAdjPv=AllAUC(AdjPv$Ragreg,AdjPv$Surv,SurvPatTest,XcovPatTest)	
	
	save(Gold,AUCGold,Lasso,LseAUC,LsminAUC,SegPval,PvalAUC,AdjPv,AUCAdjPv, file=paste("SimulationLong",i,".RData",sep=""))
	######################################"
	
	if (i<21)
	{	
		SegOdds=RunSegOdds(XcovSelect,SurvivalSelect,XcovFit,SurvivalFit)
		SegOddsLasso=RunSegLasso(XcovSelect,SurvivalSelect,XcovFit,SurvivalFit,SegOdds$MC)
		
	 
		######################################"
		
		AUCSegLassL1=AllAUC(SegOddsLasso$RagregL1,SegOddsLasso$SurvL1,SurvPatTest,XcovPatTest)
		AUCSegLassL2=AllAUC(SegOddsLasso$RagregL2,SegOddsLasso$SurvL2,SurvPatTest,XcovPatTest)
		SFKAUC=AllAUC(SegOdds$Ragreg,SegOdds$Surv,SurvPatTest,XcovPatTest)
		
		save(Gold,AUCGold,Lasso,LseAUC,LsminAUC,SegOdds,SFKAUC,SegPval,PvalAUC,AdjPv,AUCAdjPv,SegOddsLasso,AUCSegLassL1,AUCSegLassL2, file=paste("SimulationLong",i,".RData",sep=""))
	}
}



mclapply(1:100, FUN = RunAllMethods, mc.cores=3)


