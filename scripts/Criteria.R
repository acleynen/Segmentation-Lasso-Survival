Criteria<-function(Surv_model,ParamVector,Beta)
{
	Ncoeff=length(Surv_model$coefficients)
	R2=2*diff(Surv_model$loglik)/length(Surv_model$coefficients)
	ARI<-adj.rand.index(Beta, ParamVector)
	regions=unique(Beta)
	Id=rep(FALSE,length(regions)-1)
	for (i in 1:(length(regions)-1))
		if (sum((ParamVector>0)*(Beta==i))>0)
			Id[i]<-TRUE
	return(Criteria<-list(R2=R2,ARI=ARI,Id=Id,Ncoeff=Ncoeff,Surv_model=Surv_model,ParamVector=ParamVector))
}

CreateParamVector<-function(Ragreg,N)
{
	ParamVector<-rep(0,N)
	nr=nrow(Ragreg)
	for (x in 1:nr)
		ParamVector[Ragreg[x,1]:Ragreg[x,2]]<-x
	return(ParamVector)
}

Agreg<-function(p,Xcov,Ragreg)
{
	pb=rep(0,nrow(Ragreg))
	for (i in 1:nrow(Ragreg))
	{
		ind=Ragreg$start[i]:Ragreg$end[i]
		if (sum(Xcov[ind,p])>0)
			pb[i]<-1
	}	
	return(pb)
}

AgregTest<-function(p,Xcov,Ragreg)
{
	pb=rep(0,nrow(Ragreg))
	for (i in 1:nrow(Ragreg))
	{
		ind=Ragreg$start[i]:Ragreg$end[i]
		if (sum(Xcov[ind,p])>0)
			pb[i]<-1
	}	
	return(pb)
}


Xagreg<-function(Ragreg,Xcov,SurvivalPatients)
{
	n=ncol(Xcov)
	XA<-sapply(1:n,Agreg,Xcov=Xcov,Ragreg=Ragreg)
	if (!is.null(dim(XA))) XA=data.frame(XA) else XA=data.frame(t(XA))
	rownames(XA)<-paste('Region',Ragreg$start,sep="")
	colnames(XA)<-paste('Patient',1:n,sep="")

	XA=cbind("OS"=SurvivalPatients$time,"Event"=SurvivalPatients$status,t(XA))
	XA=data.frame(XA)
	return(XA)

}

XagregTest<-function(Ragreg,Xcov,SurvivalPatientsTest)
{
	XA<-sapply(1:n,AgregTest,Xcov=Xcov,Ragreg=Ragreg)
	if (!is.null(dim(XA))) XA=data.frame(XA) else XA=data.frame(t(XA))
	rownames(XA)<-paste('Region',Ragreg$start,sep="")
	colnames(XA)<-paste('PatientTest',1:n,sep="")

	XA=cbind("OS"=SurvivalPatientsTest$time,"Event"=SurvivalPatientsTest$status,t(XA))
	XA=data.frame(XA)
	return(XA)

}


EventAUC<-function(time,OS,Event)
{
	EvTime=(OS<time)
	track=(Event==0)
	rm=which(EvTime*track==1)
	return(list(EvTime=EvTime,rm=rm))
}

AUC<-function(scaledR,A)
{
	nb=length(A)/2
	Auc=c()
	for (i in 1:nb)
	{
		rm=A[[2*i]]
		if (length(rm)==0) 
		{
			Eve=A[[2*i-1]] 
			Pred=scaledR
		}	else 
		{
			Eve=A[[2*i-1]][-rm]
			Pred=scaledR[-A[[2*i]]]
		}	
		Auc=c(Auc,auc(roc(Eve,Pred)))
	}
	return(Auc)
}
