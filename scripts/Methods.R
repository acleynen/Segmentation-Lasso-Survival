################################ Estimation of univariate cox-model p-values

CoxUni<-function(X)
{
	OS=X$OS
	Event=X$Event
	SOS=Surv(OS, Event)
	
	Xp=X
	Xp$SurvObj<-SOS

	ag=sapply(3:ncol(X),function(x) {summary(coxph(SurvObj~Xp[,x],data=Xp))})
	coefg=sapply(1:(ncol(X)-2),function(x) {ag[,x]$coefficients[1]})
	pvalg=sapply(1:(ncol(X)-2),function(x) {ag[,x]$coefficients[5]})

	return(data.frame(coxcoeff=coefg,coxpval=pvalg))
}


################################ Functions for segmentation of survival odds coefficients
SlopeHeur<-function(Seg,greatjump)
{
    saut <- function(Lv, pen, Kseq, seuil = sqrt(n)/log(n), biggest = TRUE) {
        J = -Lv
        Kmax = length(J)
        k = 1
        kv = c()
        dv = c()
        pv = c()
        dmax = 1
        while (k < Kmax) {
            pk = (J[(k + 1):Kmax] - J[k])/(pen[k] - pen[(k + 
                1):Kmax])
            pm = max(pk)
            dm = which.max(pk)
            dv = c(dv, dm)
            kv = c(kv, k)
            pv = c(pv, pm)
            if (dm > dmax) {
                dmax = dm
                kmax = k
                pmax = pm
            }
            k = k + dm
        }
        if (biggest) {
            pv = c(pv, 0)
            kv = c(kv, Kmax)
            dv = diff(kv)
            dmax = max(dv)
            rt = max(dv)
            rt = which(dv == rt)
            pmax = pv[rt[length(rt)]]
            alpha = 2 * pmax
            km = kv[alpha >= pv]
            Kh = Kseq[km[1]]
            return(c(Kh, alpha))
        }
        else {
            paux <- pv[which(kv <= seuil)]
            alpha <- 2 * min(paux)
            km = kv[alpha >= pv]
            Kh = Kseq[km[1]]
            return(c(Kh, alpha))
        }
    }
  GB<-GetBreaks(Seg$BestB)
	Km=nrow(Seg$BestB)
	N=ncol(Seg$BestB)
  Kseq = 1:Km
  pen = Kseq * (2 * log(N/Kseq) + 5)

  K = saut(-Seg$BestCost[,N], pen, Kseq,biggest=greatjump)
  K <- K[1]
	return(K)
 
}

MBicCrit<-function(Seg)
{
	GB<-GetBreaks(Seg$BestB)
	Km=nrow(Seg$BestB)
	N=ncol(Seg$BestB)
	sizenr <- function(k) {
        vec <- sum(log(diff(GB[k,1:k+1])))}
	p1<-sapply(1:Km,sizenr)
	p2<-((1:Km)-0.5)*log(N)
	p3<-Seg$BestCost[,N]
	return(which.min(p1+p2+p3))
}

GetBreaks<-function(MatBreaks)
{
	Km=nrow(MatBreaks)
	Segmentation= matrix(nrow=Km,ncol=ncol(MatBreaks)+1)
	Segmentation[,1]<-1
	Segmentation[1,2]<-ncol(MatBreaks) 
	for (K in 2:Km)
	{
		Segmentation[K,K+1]<-ncol(MatBreaks)
		for (k in K:2)
			Segmentation[K,k]=MatBreaks[k,Segmentation[K,k+1]-1]
	}
	return(Segmentation)
}

OptSeg<-function(K,MatBreaks,MatParam,MatCost)
{
	Ns=ncol(MatBreaks)
	Segmentation<-MatBreaks[K,Ns]
	Likelihood<-MatCost[K,Ns]
	Param<-MatParam[Segmentation,Ns]
	PlotParam<-rep(MatParam[Segmentation,Ns],(Ns-Segmentation[1]+1))
	if (K>2)
	{
		for (k in (K-1):2)
		{
			Segmentation<-c(MatBreaks[k,Segmentation[1]-1],Segmentation)
			Param<-c(MatParam[Segmentation[1],Segmentation[2]],Param)
			PlotParam<-c(rep(MatParam[Segmentation[1],Segmentation[2]],(Segmentation[2]-Segmentation[1])),PlotParam)
		}
	}
	Param<-c(MatParam[1,Segmentation[1]-1],Param)
	PlotParam<-c(rep(MatParam[1,Segmentation[1]-1],Segmentation[1]-1),PlotParam)
	return(list(Segmentation=Segmentation,Likelihood=Likelihood,Param=Param,PlotParam=PlotParam))
}


SegMatCost<-function(matcost,Kmax)
{
	Ns=ncol(matcost)
	BestCost<-matrix(nrow=Kmax,ncol=Ns)
	BestB<-matrix(nrow=Kmax,ncol=Ns)
	BestCost[1,]<-matcost[1,]
	BestB[1,]<-1:Ns
	for (k in 2:Kmax)
		for (E in k:Ns)
		{
			C=sapply((k-1):(E-1),function(x) {BestCost[k-1,x]+matcost[x+1,E]})
			if (length(which.min(C))>0)
			{
				BestB[k,E]=k-1+which.min(C)
				BestCost[k,E]=min(C,na.rm=TRUE)
			} else
			{
				BestB[k,E]=k
				BestCost[k,E]=NA		
			}	
		}	
	return(list(BestB=BestB,BestCost=BestCost))
}


MatCost<-function(Matalpha,Alpha)
{
	Ns=ncol(Matalpha)
	CostMat=matrix(ncol=Ns,nrow=Ns)
	Param=matrix(ncol=Ns,nrow=Ns)
	for (i in 1:Ns)
	{
		if((i/10)==round(i/10))
		{
			print(Sys.time())
			print(i)
		}	
		if (i!=Ns)
		{
			subMat=Matalpha[,i:Ns]
			Valpha=apply(subMat,1,cumsum)
			CostMat[i,i:Ns]<-apply(Valpha,1,min)
			Param[i,i:Ns]<-Alpha[apply(Valpha,1,which.min)]
		} else
		{
			CostMat[Ns,Ns]<-min(Matalpha[,Ns])
			Param[Ns,Ns]<-Alpha[which.min(Matalpha[,Ns])]		
		}	
	}
	return(list(CostMat=CostMat,Param=Param))
}

MatCostAll<-function(Xcov,Alpha)
{
	Ns=nrow(Xcov)
	CostMat=matrix(ncol=Ns,nrow=Ns)
	Param=matrix(ncol=Ns,nrow=Ns)
	for (i in 1:Ns)
	{
		if((i/10)==round(i/10))
		{
			print(Sys.time())
			print(i)
		}	
		subMat=sapply(i:Ns,MinTable,segstart=i,Alpha=Alpha,Xcov=Xcov)
		CostMat[i,i:Ns]<-subMat[1,]
		Param[i,i:Ns]<-subMat[2,]
	}
	return(list(CostMat=CostMat,Param=Param))
}


Mlist<-function(i,Xcov,Alpha)
{
	Ns=nrow(Xcov)
	if((i/10)==round(i/10))
	{
		print(Sys.time())
		print(i)
	}	
	subMat=sapply(i:Ns,MinTable,segstart=i,Alpha=Alpha,Xcov=Xcov)
	return(subMat)
}



MinTable<-function(segend,segstart,Alpha,Xcov)
{
	Tab<-CreateTableSeg(segstart,segend,Xcov=Xcov)
	Ta=sapply(Alpha,TabAlpha,Tab=Tab)
	Cost=min(Ta)
	if (Cost==0) Param=1 else 	Param=Alpha[which.min(Ta)]
	return(c(Cost,Param))
}

TabAlpha<-function(x,Tab)
{
	Res=0
	if (nrow(Tab)>0)
	{
		ni=Tab[,"n.tot"]; n1i=Tab[,"n.risk.0"]; di=Tab[,"e.tot"]; d1i=Tab[,"n.event.0"]
		for (i in 1:(length(ni)-1))
		{
			n=ni[i]; n1=n1i[i]; d=di[i+1]; d1=d1i[i+1]
			p1=max(0,d-n+n1)
			p2=min(n1,d)
			if ((p2-p1)>0)
			{
				xa=x
				lP0=log(sum(sapply(p1:p2,function(p) {xa^p*choose(n1,p)*choose(n-n1,d-p)})))
				Res=Res-d1*log(x)+lP0	
			}	
		}	
	}
	return(Res) 	
}

CreateTableSeg<-function(i,j,Xcov)
{
	if (i==j) Xseg=Xcov[i:j,] else 	Xseg<-apply(Xcov[i:j,],2,prod)
	tp=as.factor(Xseg!=0)
	my.fitos=survfit(Sos~tp,data=Data)


	data.cumevents <- ggsurvtable(my.fitos,data=Data,break.time.by=1)$cumevents
	Tab=data.cumevents$data[,c("strata","time","n.risk","n.event")]
	Table0=Tab[Tab[,"strata"]=="tp=FALSE",]
	Table1=Tab[Tab[,"strata"]=="tp=TRUE",]

	TableR=cbind(Table0[,c("time","strata","n.risk","n.event")],Table1[,c("strata","n.risk","n.event")])
	TableR=cbind(TableR,n.tot=TableR[,3]+TableR[,6],e.tot=TableR[,4]+TableR[,7])
	TableR=unique(TableR[,c(3,4,8,9)])
	colnames(TableR)=c("n.risk.0","n.event.0","n.tot","e.tot")
	
	return(TableR)

}



Vecx<-function(x,TableSeg)
{
	ResVec=c()
	for (segind in 1:N)
	{
		Tab=TableSeg[[segind]]
		Res=0
		if (nrow(Tab)>0)
		{
			ni=Tab[,"n.tot"]; n1i=Tab[,"n.risk.0"]; di=Tab[,"e.tot"]; d1i=Tab[,"n.event.0"]
			for (i in 1:(length(ni)-1))
			{
				n=ni[i]; n1=n1i[i]; d=di[i+1]; d1=d1i[i+1]
				p1=max(0,d-n+n1)
				p2=min(n1,d)
				if ((p2-p1)>0)
				{
					xa=x
					lP0=log(sum(sapply(p1:p2,function(p) {xa^p*choose(n1,p)*choose(n-n1,d-p)})))
					Res=Res-d1*log(x)+lP0	
				}	
			}	
		} 	
		ResVec=c(ResVec,Res)
	}
	return(ResVec)
}

CreateTable<-function(i,Data,Xcov)
{
	tp=as.factor(Xcov[i,]==1)
	my.fitos=survfit(Sos~tp,data=Data)


	data.cumevents <- ggsurvtable(my.fitos,data=Data,break.time.by=1)$cumevents
	Tab=data.cumevents$data[,c("strata","time","n.risk","n.event")]
	Table0=Tab[Tab[,"strata"]=="tp=FALSE",]
	Table1=Tab[Tab[,"strata"]=="tp=TRUE",]

	TableR=cbind(Table0[,c("time","strata","n.risk","n.event")],Table1[,c("strata","n.risk","n.event")])
	TableR=cbind(TableR,n.tot=TableR[,3]+TableR[,6],e.tot=TableR[,4]+TableR[,7])
	colnames(TableR)=c("time","Group0","n.risk.0","n.event.0","Group1","n.risk.1","n.event.1","n.tot","e.tot")
	
	return(TableR)

}

	ChangeCostChr<-function(x,MatCost,Index)
	{
		V=MatCost[x,]
		a=min(Index[x<=Index])
		l=sum(x<=Index)
		if (a==Inf) return(V)
		V[a:length(V)]<-V[a:length(V)]+1500000*l*l
		return(V)
		
	}

