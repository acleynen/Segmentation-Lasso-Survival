library(extraDistr)
library(pdfCluster)
set.seed(1)



### True regions for Cox model

N=30000 #number of covariates
D=11 # number of relevant regions
sD=c(5,15,750,30,120,5,15,30,60,120,15) # length of each relevant region
betaD=c(sort(rsign(9)*rnorm(9,1,0.4)),1,1) # values of beta for Cox model for whole region
startD=sort(c(500+2500*c(0:8,10),4250)) # location start of each region

# summary of all true regions of interest
Regions=data.frame(start=startD,length=sD,beta=betaD,ebeta=exp(betaD)) 


################ Test Datasets
set.seed(1234)

### Patient covariates

n=500 # number of patients
FreqR=sort(c(startD[-3]-15,500,1000,5400,7200,8000,13000,18200,22200))
LR=c(150,150,300,300,100,100,200,200,30,200,50,100,150,200,500,15,50,50)

# first half of patients, share long events

Rshare=c(3500,7000,12500,21000) 
Lrshare=1500
restoftheworld1=(1:N)[-c(unlist(sapply(Rshare,function(x){x:(x+Lrshare)})),unlist(sapply(FreqR,function(x){(x-15):(x+15)})))]

Nsh=sample(1:4,n/2,replace=T) # number of shared long events for each patient
Nad=rpois(n/2,1) # number of random event for each patient

Covhalf1<-function(p)
{
	Xp=rep(0,N)
	share=sort(sample(1:4,Nsh[p],replace=F)) # which long event
	for (i in 1:Nsh[p])
	{
		uni=runif(1)
		if (uni<0.6) # 60% patient have whole-length event
		{
			Xp[Rshare[share[i]]:(Rshare[share[i]]+Lrshare)]<-1
		} else # 40% patient have only a subset of event
		{
			if ((share[i]==1)|| share[i]==3)
				prop=round(Lrshare*0.5*runif(1)) else prop=round(Lrshare*(0.5+abs(runif(1)-0.5)))
			Xp[Rshare[share[i]]:(Rshare[share[i]]+prop)]<-1
		}
	}
	if(sample(c(TRUE,FALSE),1)) # if true patient also has one short frequent event
	{
		rs=sample(1:length(FreqR),1)
		start=FreqR[rs]+sample((-15):15,1)
		Xp[start:(start+round(LR[rs]*(0.5+abs(runif(1)-0.5))))]<-1
	}
	if (Nad[p]!=0)
	{
		rs=sample(restoftheworld1,Nad[p])
		for (i in 1:length(rs))
			Xp[rs[i]:min(N,(rs[i]+sample(c(15,50,75),1)))]<-1
	}			
	return(Xp)	
}


# second half of patients, have many short events

restoftheworld2=(1:N)[-unlist(sapply(FreqR,function(x){(x-15):(x+15)}))]

Nfreq=rpois(n/2,2) # number of frequent short events, for each patient
Nad2=rpois(n/2,2) # number of random events, for each patient

Covhalf2<-function(p)
{
	Xp=rep(0,N)

	if(Nfreq[p]!=0)
	{
		for (i in 1:length(Nfreq[p]))
		{
			rs=sample(1:length(FreqR),1)
			ls=rpois(1,15) # 
			start=FreqR[rs]+sample((-ls):ls,1) # random starting point of event
			Xp[start:(start+round(LR[rs]*(0.5+abs(runif(1)-0.5))))]<-1 # length of event
		}	
	}
	if (Nad2[p]!=0)
	{
		for (i in 1:length(Nad2[p]))
		{	
			ls=rpois(1,30)
			rs=sample(restoftheworld2,1)
			Xp[rs:min(N,(rs+ls))]<-1
		}	
	}			
	return(Xp)	
}


##### survival functions

# extract coefficient values for each patient based on CNV event
coxb<-function(p, XCov)
{
	pb=0
	for (i in 1:D)
	{
		ind=Regions$start[i]:(Regions$start[i]+Regions$length[i]-1)
		if (sum(XCov[ind,p])>0)
			pb=pb+Regions$beta[i]
	}	
	return(pb)
}

# Simulation of lifetime duration from Cox model with exponential censoring
simulCox <- function(lambda, rho, beta, rateC) 
{
  Nsimul=length(beta)
  
  # Weibull latent event times
  v <- runif(Nsimul)
  Tlat <- (- log(v) / (lambda * exp(beta)))^(1 / rho)

  # censoring times
  C <- rexp(Nsimul, rate=rateC)

  # follow-up times and event indicators
  time <- pmin(Tlat, C)
  status <- as.numeric(Tlat <= C)

  # data set
  data.frame(id=1:Nsimul,
             time=time,
             status=status,
             beta=beta)
}


lambda=0.01 # scale parameter of baseline function
rho=0.75 # shape parameter of baseline function
rateC=0.0005 # rate parameter of censoring variable



##############################################################

NTest=100
XcovTest<-c()
SurvivalPatientsTest<-c()

for (ntest in 1:NTest)
{
	if ((ntest %% 10)==0) print(ntest)
	XcovTest[[ntest]]=cbind(sapply(1:(n/2),Covhalf1),sapply(1:(n/2),Covhalf2)) # la matrice de design (présence du CNV) des patients
	rownames(XcovTest[[ntest]])<-paste('Region',1:N,sep="")
	colnames(XcovTest[[ntest]])<-paste('PatientTest',1:n,sep="")


	CoxebetaTest=sapply(1:n,coxb,XCov=XcovTest[[ntest]])

	SurvivalPatientsTest[[ntest]]=simulCox(lambda,rho,CoxebetaTest,rateC)
}

#########
save(XcovTest,SurvivalPatientsTest,Regions, file="Test30000.RData")
#########

##################################################################
# Downscale dataset

Nd=1000
fact=30000/Nd

RegionsNd=data.frame(start=floor(Regions$start/fact), length=floor((Regions$start+Regions$length-1)/fact)+1-floor(Regions$start/fact), beta=Regions$beta, ebeta=Regions$ebeta)


Downscale<-function(start,Xcov)
{
	Xtemp<-Xcov[start:(start+29),]
	V=apply(Xtemp,2, sum)
	return(as.numeric(V>14))
}

start=seq(0,29970,by=30)+1

XcovTest2<-c()
for (i in 1:length(XcovTest))
{
	if ((i %% 10)==0) print(i)
	XcovNd<-sapply(start,Downscale,Xcov=XcovTest[[i]])
	XcovTest2[[i]]<-t(XcovNd)
}

########

XcovTest=XcovTest2
N=Nd
Regions=RegionsNd
save(XcovTest,SurvivalPatientsTest,Regions, file="Test1000_downscaled.RData")

###################################################################################################################################################################
## Train datasets were simulated with exact same codes replacing the seed with
set.seed(12345)




