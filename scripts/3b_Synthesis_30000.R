
#
GoldR2=c()
GoldNcoeff=c()
GoldARI=c()
GoldId=c()
GoldAUC=c()

LsminR2=c()
LsminNcoeff=c()
LsminARI=c()
LsminId=c()
TLsminAUC=c()

LseR2=c()
LseNcoeff=c()
LseARI=c()
LseId=c()
TLseAUC=c()

PvalR2=c()
PvalNcoeff=c()
PvalARI=c()
PvalId=c()
TPvalAUC=c()

AdjPvR2=c()
AdjPvNcoeff=c()
AdjPvARI=c()
AdjPvId=c()
AdjPvAUC=c()

SOddsR2=c()
SOddsNcoeff=c()
SOddsARI=c()
SOddsId=c()
TSOddsAUC=c()

SLasminR2=c()
SLasminNcoeff=c()
SLasminARI=c()
SLasminId=c()
SLasminAUC=c()

SLasseR2=c()
SLasseNcoeff=c()
SLasseARI=c()
SLasseId=c()
SLasseAUC=c()


for (i in 1:100)
{
	load(paste("SimulationLong",i,".RData",sep=""))

	GoldR2=c(GoldR2,Gold$R2)
	GoldNcoeff=c(GoldNcoeff,Gold$Ncoeff)
	GoldARI=c(GoldARI,Gold$ARI)
	GoldId=cbind(GoldId,Gold$Id)
	GoldAUC=cbind(GoldAUC,AUCGold)

	

	LsminR2=c(LsminR2,Lasso$R2L1)
	LsminNcoeff=c(LsminNcoeff,Lasso$NcoeffL1)
	LsminARI=c(LsminARI,Lasso$ARIL1)
	LsminId=cbind(LsminId,Lasso$IdL1)
	TLsminAUC=cbind(TLsminAUC,LsminAUC)
	LseR2=c(LseR2,Lasso$R2L2)
	LseNcoeff=c(LseNcoeff,Lasso$NcoeffL2)
	LseARI=c(LseARI,Lasso$ARIL2)
	LseId=cbind(LseId,Lasso$IdL2)
	TLseAUC=cbind(TLseAUC,LseAUC)



	PvalR2=c(PvalR2,SegPval$R2)
	PvalNcoeff=c(PvalNcoeff,SegPval$Ncoeff)
	PvalARI=c(PvalARI,SegPval$ARI)
	PvalId=cbind(PvalId,SegPval$Id)
	TPvalAUC=cbind(TPvalAUC,PvalAUC)
		


	AdjPvR2=c(AdjPvR2,AdjPv$R2)
	AdjPvNcoeff=c(AdjPvNcoeff,AdjPv$Ncoeff)
	AdjPvARI=c(AdjPvARI,AdjPv$ARI)
	AdjPvId=cbind(AdjPvId,AdjPv$Id)
	AdjPvAUC=cbind(AdjPvAUC,AUCAdjPv)				
	
	if (i<21)
	{
		SR2=SegOdds$R2
		SN=SegOdds$Ncoeff
		SARI=SegOdds$ARI
		SId=SegOdds$Id
		SLR2=SegOddsLasso$R2L1
		SLN=SegOddsLasso$NcoeffL1
		SLARI=SegOddsLasso$ARIL1
		SLId=SegOddsLasso$IdL1
		

		SOddsR2=c(SOddsR2,SR2)
		SOddsNcoeff=c(SOddsNcoeff,SN)
		SOddsARI=c(SOddsARI,SARI)
		SOddsId=cbind(SOddsId,SId)
		TSOddsAUC=cbind(TSOddsAUC,SFKAUC)				
		
		SLasminR2=c(SLasminR2,SLR2)
		SLasminNcoeff=c(SLasminNcoeff,SLN)
		SLasminARI=c(SLasminARI,SLARI)
		SLasminId=cbind(SLasminId,SLId)
		SLasminAUC=cbind(SLasminAUC,AUCSegLassL1)	
	
	}
		
	
}

save(Sim,TPvalAUC,PvalId,PvalARI,PvalNcoeff,PvalR2,TLseAUC,LseId, LseARI,LseNcoeff,LseR2,TLsminAUC,LsminId,LsminARI,LsminNcoeff,LsminR2,GoldAUC,GoldId,GoldARI,GoldNcoeff,GoldR2, AdjPvR2,AdjPvNcoeff,AdjPvARI,AdjPvId,AdjPvAUC,SOddsR2,SOddsNcoeff,SOddsARI,SOddsId,TSOddsAUC,SLasminR2,SLasminNcoeff,SLasminARI,SLasminId,SLasminAUC,SLasseR2,SLasseNcoeff,SLasseARI,SLasseId,SLasseAUC, file="NewEvalLongAll.RData")


