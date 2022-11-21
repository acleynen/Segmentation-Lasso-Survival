#################################################"
# gold standard for Cox-Model

RunGold<-function(XcovS,SurvivalS,XcovF,SurvivalF)
{

	Ragreg<-data.frame(start=Regions$start,	end=Regions$start+Regions$length-1)
	XGold<-Xagreg(Ragreg,XcovF,SurvivalF)
	Surv_coxGold <- coxph(Surv(OS, Event) ~ ., XGold)
	
	GrGold<-CreateParamVector(Ragreg,N)
	RGold<-Criteria(Surv_coxGold,GrGold,BetaClass)

	RagregGold=Ragreg
	return(list(R2=RGold$R2,ARI=RGold$ARI,Ncoeff=RGold$Ncoeff,Id=RGold$Id,Ragreg=RagregGold,Surv=Surv_coxGold))
}

#################################################"
# Lasso-penalized cox
RunLasso<-function(XcovS,SurvivalS,XcovF,SurvivalF)
{

	x=t(XcovS)


	fit<-glmnet(x,Surv(time=SurvivalS$time,event=SurvivalS$status),family="cox",alpha=1)
	cv.fit<-cv.glmnet(x,Surv(time=SurvivalS$time,event=SurvivalS$status),family="cox",alpha=1)

	Coefficients <- coef(fit, s = cv.fit$lambda.min)
	 Active.Index <- which(Coefficients != 0)
	Coefficients2 <- coef(fit, s = cv.fit$lambda.1se)
	 Active.Index2 <- which(Coefficients2 != 0)


	#########################
	# Active2 : selection with 1se penalty
	########################"

	da=diff(Active.Index2)
	Ragreg<-data.frame(start=Active.Index2[c(1,which(da>5)+1)],	end=Active.Index2[c(which(da>5),length(Active.Index2))])
	print(Ragreg)	
	
	XLA2<-Xagreg(Ragreg,XcovF,SurvivalF)

	Nb=apply(XLA2,2,sum)
	rm=which(Nb<2)
	if (length(rm)>0)	Xtemp=XLA2[,-rm] else Xtemp=XLA2
	Surv_L2 <- coxph(Surv(OS, Event) ~ ., Xtemp)
	if (sum(is.na(Surv_L2$coefficients))>0)
	{
		rm=which(is.na(Surv_L2$coefficients))
		Xtemp=Xtemp[,-(rm+2)]
		Surv_L2 <- coxph(Surv(OS, Event) ~ ., Xtemp)
	}	
	pv=summary(Surv_L2)$coefficients[,5]
	a=sum(pv>0.05)

	while (a>0)
	{
		b=which.max(pv)
		Xtemp=Xtemp[,-(b+2)]
		Surv_L2 <- coxph(Surv(OS, Event) ~ ., Xtemp)
		pv=summary(Surv_L2)$coefficients[,5]
		a=sum(pv>0.05)		
	}

	Ragreg<-Ragreg[is.element(paste("Region",Ragreg[,1],sep=""),colnames(Xtemp)),]
	GrL2<-CreateParamVector(Ragreg,N)
	RL2<-Criteria(Surv_L2,GrL2,BetaClass)
	RagregL2=Ragreg
	ColL2<-rep(0,N)
	for (i in 1:nrow(Ragreg))
		ColL2[Ragreg[i,1]:Ragreg[i,2]]<-1		

	
	#########################
	# Active1: selection with smin penalty
	########################"

	da=diff(Active.Index)
	Ragreg<-data.frame(start=Active.Index[c(1,which(da>5)+1)],	end=Active.Index[c(which(da>5),length(Active.Index))])
	XLA1<-Xagreg(Ragreg,XcovF,SurvivalF)

	Nb=apply(XLA1,2,sum)
	rm=which(Nb<2)
	if (length(rm)>0)	Xtemp=XLA1[,-rm] else Xtemp=XLA1
	Surv_L1 <- coxph(Surv(OS, Event) ~ ., Xtemp)
	if (sum(is.na(Surv_L1$coefficients))>0)
	{
		rm=which(is.na(Surv_L1$coefficients))
		Xtemp=Xtemp[,-(rm+2)]
		Surv_L1 <- coxph(Surv(OS, Event) ~ ., Xtemp)
	}	
	pv=summary(Surv_L1)$coefficients[,5]
	a=sum(pv>0.05)

	while (a>0)
	{
		b=which.max(pv)
		Xtemp=Xtemp[,-(b+2)]
		Surv_L1 <- coxph(Surv(OS, Event) ~ ., Xtemp)
		pv=summary(Surv_L1)$coefficients[,5]
		a=sum(pv>0.05)		
	}

	Ragreg<-Ragreg[is.element(paste("Region",Ragreg[,1],sep=""),colnames(Xtemp)),]

	GrL1<-CreateParamVector(Ragreg,N)
	RL1<-Criteria(Surv_L1,GrL1,BetaClass)
	RagregL1=Ragreg
	ColL1<-rep(0,N)
	for (i in 1:nrow(Ragreg))
		ColL1[Ragreg[i,1]:Ragreg[i,2]]<-1

	
	return(list(R2L1=RL1$R2,ARIL1=RL1$ARI,NcoeffL1=RL1$Ncoeff,RagregL1=RagregL1, IdL1=RL1$Id, IdL2=RL2$Id, ColL1=ColL1,R2L2=RL2$R2,ARIL2=RL2$ARI,NcoeffL2=RL2$Ncoeff,RagregL2=RagregL2, ColL2=ColL2, Surv_L1=Surv_L1, Surv_L2=Surv_L2))	
}
##################################
# Multiple testing correction
##################################

RunMultipleTesting<-function(XcovS,SurvivalS,XcovF,SurvivalF)
{
	xsecond=cbind("OS"=SurvivalS$time,"Event"=SurvivalS$status,t(XcovS))
	xsecond=data.frame(xsecond)

	CoxU=CoxUni(xsecond)
	
	AdjPv=p.adjust(CoxU$coxpval, method ="BH")
	ActAdjPv=which(AdjPv<0.05) 

	da=diff(ActAdjPv)
	Ragreg<-data.frame(start=ActAdjPv[c(1,which(da>5)+1)],	end=ActAdjPv[c(which(da>5),length(ActAdjPv))])
	XApv<-Xagreg(Ragreg,XcovF,SurvivalF)

	Nb=apply(XApv,2,sum)
	rm=which(Nb<2)
	if (length(rm)>0)	Xtemp=XApv[,-rm] else Xtemp=XApv
	Surv_XApv <- coxph(Surv(OS, Event) ~ ., Xtemp)
	if (sum(is.na(Surv_XApv$coefficients))>0)
	{
		rm=which(is.na(Surv_XApv$coefficients))
		Xtemp=Xtemp[,-(rm+2)]
		Surv_XApv <- coxph(Surv(OS, Event) ~ ., Xtemp)
	}	


	pv=summary(Surv_XApv)$coefficients[,5]
	a=sum(pv>0.05)

	while (a>0)
	{
		b=which.max(pv)
		Xtemp=Xtemp[,-(b+2)]
		Surv_XApv <- coxph(Surv(OS, Event) ~ ., Xtemp)
		pv=summary(Surv_XApv)$coefficients[,5]
		a=sum(pv>0.05)		
	}

	Ragreg<-Ragreg[is.element(paste("Region",Ragreg[,1],sep=""),colnames(Xtemp)),]


	GrPv<-CreateParamVector(Ragreg,N)
	RLPv<-Criteria(Surv_XApv,GrPv,BetaClass)

	BrPv<-sort(unique(c(Ragreg[,1],Ragreg[,2]+1)))

	RagregPv=Ragreg
	ColPv<-rep(0,N)
	for (i in 1:nrow(Ragreg))
		ColPv[Ragreg[i,1]:Ragreg[i,2]]<-1
		
	return(list(R2=RLPv$R2,ARI=RLPv$ARI,Ncoeff=RLPv$Ncoeff, Id=RLPv$Id, Ragreg=RagregPv,Col=ColPv,Surv=Surv_XApv))
	
}



##################################
# Segmentation of Cox-univariate p-values
##################################

RunCoxUni<-function(XcovS,SurvivalS,XcovF,SurvivalF)
{
	xsecond=cbind("OS"=SurvivalS$time,"Event"=SurvivalS$status,t(XcovS))
	xsecond=data.frame(xsecond)

	CoxU=CoxUni(xsecond)

	Pval=log(CoxU$coxpval/(1-CoxU$coxpval))
	Pval[is.na(Pval)]<-3
	SP=Segmentor(Pval,model=2,Kmax=150)
	Ksp=SelectModel(SP)
	BSP=c(0,getBreaks(SP)[Ksp,1:Ksp])
	ParamGChr=getParameters(SP)[Ksp,1:Ksp]
	e=which(ParamGChr<log(0.01/0.99))
	Ragreg<-data.frame(start=BSP[e]+1,end=BSP[e+1])
	XCUA<-Xagreg(Ragreg,XcovF,SurvivalF)

	Nb=apply(XCUA,2,sum)
	rm=which(Nb<2)
	if (length(rm)>0)	Xtemp=XCUA[,-rm] else Xtemp=XCUA
	Surv_XCUA <- coxph(Surv(OS, Event) ~ ., Xtemp)
	
	if (sum(is.na(Surv_XCUA$coefficients))>0)
	{
		rm=which(is.na(Surv_XCUA$coefficients))
		Xtemp=Xtemp[,-(rm+2)]
		Surv_XCUA <- coxph(Surv(OS, Event) ~ ., Xtemp)
	}	

	
	
	pv=summary(Surv_XCUA)$coefficients[,5]
	a=sum(pv>0.05)

	while (a>0)
	{
		b=which.max(pv)
		Xtemp=Xtemp[,-(b+2)]
		Surv_XCUA <- coxph(Surv(OS, Event) ~ ., Xtemp)
		pv=summary(Surv_XCUA)$coefficients[,5]
		a=sum(pv>0.05)		
	}

	Ragreg<-Ragreg[is.element(paste("Region",Ragreg[,1],sep=""),colnames(Xtemp)),]


	GrPv<-CreateParamVector(Ragreg,N)
	RLPv<-Criteria(Surv_XCUA,GrPv,BetaClass)

	BrPv<-sort(unique(c(Ragreg[,1],Ragreg[,2]+1)))

	RagregPv=Ragreg
	ColPv<-rep(0,N)
	for (i in 1:nrow(Ragreg))
		ColPv[Ragreg[i,1]:Ragreg[i,2]]<-1
		
	return(list(R2=RLPv$R2,ARI=RLPv$ARI,Ncoeff=RLPv$Ncoeff, Id=RLPv$Id, Ragreg=RagregPv,Col=ColPv,Surv=Surv_XCUA))
	
}




##################################
# Segmentation of Fisher parameter
##################################

RunSegOdds<-function(XcovS,SurvivalS,XcovF,SurvivalF) # watch out, super long!
{
	Sos=Surv(time=SurvivalS$time,event=SurvivalS$status)
	Data=data.frame(Sos)
	alpha=c(seq(0.01,2,by=0.01),seq(2.1,3,by=0.1),seq(3.5,5,by=0.5))
	N=nrow(XcovS)
	TableSeg=lapply(1:N,CreateTable,Data=Data,Xcov=XcovS)
	A=t(sapply(alpha,Vecx,TableSeg=TableSeg))
	MC=MatCost(A,alpha)	
	SegCO<-SegMatCost(MC$CostMat,40)
	OptS=OptSeg(40,MatBreaks=SegCO$BestB,MatParam=MC$Param,MatCost=SegCO$BestCost)


	## K selected with slope heuristic
	e<-which(abs(log(OptS$Param))>log2(1.25))
	if (is.element(1,e))
	{
		Ragreg<-data.frame(start=c(1,OptS$Segmentation[e-1]),end=c(OptS$Segmentation,N+1)[e]-1)
	} else Ragreg<-data.frame(start=c(OptS$Segmentation[e-1]),end=c(OptS$Segmentation,N+1)[e]-1)

	XOptS<-Xagreg(Ragreg,XcovF,SurvivalF)
	Nb=apply(XOptS,2,sum)
	rm=which(Nb<2)
	if (length(rm)>0)	Xtemp=XOptS[,-rm] else Xtemp=XOptS
	Surv_XOptS <- coxph(Surv(OS, Event) ~ ., Xtemp)
	if (sum(is.na(Surv_XOptS$coefficients))>0)
	{
		rm=which(is.na(Surv_XOptS$coefficients))
		Xtemp=Xtemp[,-(rm+2)]
		Surv_XOptS <- coxph(Surv(OS, Event) ~ ., Xtemp)
	}	
	
	pv=summary(Surv_XOptS)$coefficients[,5]
	a=sum(pv>0.05)

	while (a>0)
	{
		b=which.max(pv)
		Xtemp=Xtemp[,-(b+2)]
		Surv_XOptS <- coxph(Surv(OS, Event) ~ ., Xtemp)
		pv=summary(Surv_XOptS)$coefficients[,5]
		a=sum(pv>0.05)		
	}

	Ragreg<-Ragreg[is.element(paste("Region",Ragreg[,1],sep=""),colnames(Xtemp)),]
	GrSegF<-CreateParamVector(Ragreg,N)
	RLSF<-Criteria(Surv_XOptS,GrSegF,BetaClass)
	RagregSF=Ragreg
	ColSF<-rep(0,N)
	for (i in 1:nrow(Ragreg))
		ColSF[Ragreg[i,1]:Ragreg[i,2]]<-1
	
	return(list(R2=RLSF$R2,ARI=RLSF$ARI,Ncoeff=RLSF$Ncoeff,Id=RLSF$Id, Ragreg=RagregSF,Col=ColSF,Surv=Surv_XOptS,MC=MC))
}


##################################
# Segmentation odds + Lasso 
##################################

RunSegLasso<-function(XcovS,SurvivalS,XcovF,SurvivalF,MC) # watch out, super long!
{

	N=nrow(XcovS)
	SegCO<-SegMatCost(MC$CostMat,300)
	OptS=OptSeg(300,MatBreaks=SegCO$BestB,MatParam=MC$Param,MatCost=SegCO$BestCost)
	Ragreg<-data.frame(start=c(1,OptS$Segmentation),end=c(OptS$Segmentation,N+1)-1)
	XS<-Xagreg(Ragreg,XcovS,SurvivalS)
	RagregSeg=Ragreg
	
	x=as.matrix(XS[,-c(1,2)])
	fit<-glmnet(x,Surv(time=SurvivalS$time,event=SurvivalS$status),family="cox",alpha=1)
	cv.fit<-cv.glmnet(x,Surv(time=SurvivalS$time,event=SurvivalS$status),family="cox",alpha=1)
	Coefficients <- coef(fit, s = cv.fit$lambda.min)
	Active.Index <- which(Coefficients != 0)
	Coefficients2 <- coef(fit, s = cv.fit$lambda.1se)
	 Active.Index2 <- which(Coefficients2 != 0)
	#smin
	da=diff(Active.Index)
	Ragreg<-data.frame(start=RagregSeg$start[Active.Index[c(1,which(da>0)+1)]],	end=RagregSeg$end[Active.Index[c(which(da>0),length(Active.Index))]])
	XLA1<-Xagreg(Ragreg,XcovF,SurvivalF)

	Nb=apply(XLA1,2,sum)
	rm=which(Nb<2)
	if (length(rm)>0)	Xtemp=XLA1[,-rm] else Xtemp=XLA1
	Surv_L1 <- coxph(Surv(OS, Event) ~ ., Xtemp)
	if (sum(is.na(Surv_L1$coefficients))>0)
	{
		rm=which(is.na(Surv_L1$coefficients))
		Xtemp=Xtemp[,-(rm+2)]
		Surv_L1 <- coxph(Surv(OS, Event) ~ ., Xtemp)
	}		
	pv=summary(Surv_L1)$coefficients[,5]
	a=sum(pv>0.05)

	while (a>0)
	{
		b=which.max(pv)
		Xtemp=Xtemp[,-(b+2)]
		Surv_L1 <- coxph(Surv(OS, Event) ~ ., Xtemp)
		pv=summary(Surv_L1)$coefficients[,5]
		a=sum(pv>0.05)		
	}

	Ragreg<-Ragreg[is.element(paste("Region",Ragreg[,1],sep=""),colnames(Xtemp)),]

	GrL1<-CreateParamVector(Ragreg,N)
	RL1<-Criteria(Surv_L1,GrL1,BetaClass)
	RagregL1=Ragreg
	ColL1<-rep(0,N)
	for (i in 1:nrow(Ragreg))
		ColL1[Ragreg[i,1]:Ragreg[i,2]]<-1

	#1se
	da=diff(Active.Index2)
	Ragreg<-data.frame(start=RagregSeg$start[Active.Index2[c(1,which(da>0)+1)]],	end=RagregSeg$end[Active.Index2[c(which(da>0),length(Active.Index2))]])
	XLA2<-Xagreg(Ragreg,XcovF,SurvivalF)

	Nb=apply(XLA2,2,sum)
	rm=which(Nb<2)
	if (length(rm)>0)	Xtemp=XLA2[,-rm] else Xtemp=XLA1
	Surv_L2 <- coxph(Surv(OS, Event) ~ ., Xtemp)
	if (sum(is.na(Surv_L2$coefficients))>0)
	{
		rm=which(is.na(Surv_L2$coefficients))
		Xtemp=Xtemp[,-(rm+2)]
		Surv_L2 <- coxph(Surv(OS, Event) ~ ., Xtemp)
	}		
	pv=summary(Surv_L2)$coefficients[,5]
	a=sum(pv>0.05)

	while (a>0)
	{
		b=which.max(pv)
		Xtemp=Xtemp[,-(b+2)]
		Surv_L2 <- coxph(Surv(OS, Event) ~ ., Xtemp)
		pv=summary(Surv_L2)$coefficients[,5]
		a=sum(pv>0.05)		
	}

	Ragreg<-Ragreg[is.element(paste("Region",Ragreg[,1],sep=""),colnames(Xtemp)),]

	GrL2<-CreateParamVector(Ragreg,N)
	RL2<-Criteria(Surv_L2,GrL2,BetaClass)
	RagregL2=Ragreg
	ColL2<-rep(0,N)
	for (i in 1:nrow(Ragreg))
		ColL2[Ragreg[i,1]:Ragreg[i,2]]<-1
		
	return(list(R2L1=RL1$R2,ARIL1=RL1$ARI,NcoeffL1=RL1$Ncoeff,RagregL1=RagregL1, IdL1=RL1$Id, ColL1=ColL1, SurvL1=Surv_L1, R2L2=RL2$R2,ARIL2=RL2$ARI,NcoeffL2=RL2$Ncoeff,RagregL2=RagregL2, IdL2=RL2$Id, ColL2=ColL2, SurvL2=Surv_L2))	
}

################################################
# Evaluate on Test data
################################################

AllAUC<-function(Ragreg,Surv_Model,SurvivalPatientsTest,XcovTest)
{
		EvTime=sapply(seq(50,2000,by=25),EventAUC,OS=SurvivalPatientsTest$time,Event=SurvivalPatientsTest$status)
		if (sum(is.na(Surv_Model$coefficients))>0)
		{
			Ragreg<-Ragreg[-which(is.na(Surv_Model$coefficients)),]
			Surv_Model$coefficients<-Surv_Model$coefficients[-which(is.na(Surv_Model$coefficients))]
		}
		if (nrow(Ragreg)>0)
		{
			XTest<-XagregTest(Ragreg,XcovTest,SurvivalPatientsTest)
			if (nrow(Ragreg)==1) Risk<-((XTest*Surv_Model$coefficients)[,3:ncol(XTest)]) else 		Risk<-apply((XTest*Surv_Model$coefficients)[,3:ncol(XTest)],1,sum)
			scaledRisk=(Risk-min(Risk))/(max(Risk)-min(Risk))
			AucMeth<-AUC(scaledRisk,EvTime)
		} else AucMeth<-rep(0,length(	seq(50,2000,by=25)))
		return(AucMeth)
}


