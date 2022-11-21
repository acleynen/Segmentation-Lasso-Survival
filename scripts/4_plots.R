library(glmnet)
library(survival)
library(survminer)
library(RColorBrewer)
library(edgeR)
library(hrbrthemes)
library(plotly)
library(ggpubr)
library(ggrepel)
library(Segmentor3IsBack)
library(pdfCluster)
library(zoo)
library(extrafont)
library(xtable)
library(pROC)
library(Rmisc)
library(scales)


source("scripts/Methods.R")
source("scripts/Criteria.R")
source("scripts/RunMethods.R")
set.seed(12341234)


hex <- hue_pal()(6)

scale_fill_chris <- function(...){
    ggplot2:::manual_scale(
        'color', 
        values = setNames(hex, c("Gold","Lasso","Adj-Pval","Seg Pval","Seg Odds","Seg+Lasso")),
        ...
    )
}
scale_fill_chris2 <- function(...){
    ggplot2:::manual_scale(
        'fill', 
        values = setNames(hex, c("Gold","Lasso","Adj-Pval","Seg Pval","Seg Odds","Seg+Lasso")),
        ...
    )
}

#####################
 ## Downscaled
####################

load("NewEvalDownscaledAll.RData")
Nexp=length(Sim)

SumSFAUC<-apply(TSOddsAUC,1,CI)
SumPvalAUC<-apply(TPvalAUC,1,CI)
SumAdjPvAUC<-apply(AdjPvAUC,1,CI)
SumLseAUC<-apply(TLseAUC,1,CI)
SumLsminAUC<-apply(TLsminAUC,1,CI)
SumGoldAUC<-apply(GoldAUC,1,CI)
SumSFLassAUC<-apply(SLasminAUC,1,CI)
SumSFLassL2AUC<-apply(SLasseAUC,1,CI)

neach=length(seq(50,2000,by=25))

Auc2=data.frame(Time=rep(seq(50,2000,by=25),6), Auc=c(SumGoldAUC[2,],SumLsminAUC[2,],SumAdjPvAUC[2,],SumPvalAUC[2,],SumSFAUC[2,],SumSFLassAUC[2,]),  upper=c(SumGoldAUC[1,],SumLsminAUC[1,],SumAdjPvAUC[1,],SumPvalAUC[1,],SumSFAUC[1,],SumSFLassAUC[1,]), lower=c(SumGoldAUC[3,],SumLsminAUC[3,],SumAdjPvAUC[3,],SumPvalAUC[3,],SumSFAUC[3,],SumSFLassAUC[3,]),  Method=factor(c(rep("Gold",neach),rep("Lasso",neach),rep("Adj-Pval",neach),rep("Seg Pval",neach),rep("Seg Odds",neach), rep("Seg+Lasso",neach)), levels=c("Gold","Lasso","Adj-Pval","Seg Pval","Seg Odds","Seg+Lasso")))


GAucDown<-ggplot(Auc2, aes(x=Time,y=Auc,group=Method,color=Method))+
		geom_line(size=2, alpha=0.9) +
		geom_line(aes(y=upper,color=Method,group=Method), size=1, alpha=0.9,linetype="dotted") +
		geom_line(aes(y=lower,color=Method,group=Method), size=1, alpha=0.9,linetype="dotted") +
		theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold"),
        legend.key.size = unit(1.5, 'cm'),
        legend.text = element_text(size=18),
        legend.position="none") +
		ylim(0.5,0.625) +
		ylab("AUC") 

GAucDown





AllR2D=data.frame(R2=c(GoldR2,  LsminR2, AdjPvR2,PvalR2, SOddsR2,SLasminR2), Method=factor(c(rep("Gold",Nexp),rep("Lasso",Nexp),rep("Adj-Pval",Nexp),rep("Seg Pval",Nexp), rep("Seg Odds",Nexp),rep("Seg+Lasso",Nexp)), levels=c("Gold","Lasso","Adj-Pval", "Seg Pval","Seg Odds","Seg+Lasso")))
AllARID=data.frame(ARI=c(GoldARI,  LsminARI, AdjPvARI, PvalARI, SOddsARI,SLasminARI), Method=factor(c(rep("Gold",Nexp),rep("Lasso",Nexp),rep("Adj-Pval",Nexp),rep("Seg Pval",Nexp), rep("Seg Odds",Nexp),rep("Seg+Lasso",Nexp)), levels=c("Gold","Lasso","Adj-Pval", "Seg Pval","Seg Odds","Seg+Lasso")))
AllNcoeffD=data.frame(Ncoeff=c(GoldNcoeff,LsminNcoeff, AdjPvNcoeff, PvalNcoeff,SOddsNcoeff,SLasminNcoeff), Method=factor(c(rep("Gold",Nexp),rep("Lasso",Nexp),rep("Adj-Pval",Nexp),rep("Seg Pval",Nexp), rep("Seg Odds",Nexp),rep("Seg+Lasso",Nexp)), levels=c("Gold","Lasso","Adj-Pval", "Seg Pval","Seg Odds","Seg+Lasso")))
AllIDD=data.frame(Id=c(apply(GoldId,1,sum),  apply(LsminId,1,sum), apply(AdjPvId,1,sum), apply(PvalId,1,sum), apply(SOddsId,1,sum), apply(SLasminId,1,sum)), Method=factor(c(rep("Gold",11),rep("Lasso",11), rep("Adj-Pval",11), rep("Seg Pval",11),rep("Seg Odds",11),rep("Seg+Lasso",11)),levels=c("Gold","Lasso","Adj-Pval","Seg Pval","Seg Odds","Seg+Lasso")), Regions=factor(rep(paste("Region",Regions$start,sep=""),6),levels=paste("Region",Regions$start,sep="")))

GR2D<-ggplot(AllR2D, aes(x=Method, y=R2, color=Method)) +
	geom_boxplot() +
	theme(axis.text=element_text(size=18),
				legend.key.size = unit(2, 'cm'), legend.text = element_text(size=16), legend.box = "horizontal",
				legend.position="none",
				axis.title.y = element_text(size = 18),
				axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
				ylim(0,40) 

GARID<-ggplot(AllARID, aes(x=Method, y=ARI, color=Method)) +
  geom_boxplot()+
	theme(axis.text=element_text(size=18),
			axis.title.y = element_text(size = 18),
			legend.position="none",
				#legend.position="bottom", legend.key.size = unit(2, 'cm'), legend.text = element_text(size=16),
				axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

GNcoeffD<-ggplot(AllNcoeffD, aes(x=Method, y=Ncoeff, color=Method)) +
  geom_boxplot() +
	theme(axis.text=element_text(size=18),
				legend.position="none",
				axis.title.y = element_text(size = 18),
				axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
				ylab("Nb of coeff")

GIDD<-ggplot(AllIDD, aes(x=Regions, y=Id, fill=Method)) +
  geom_bar(stat="identity",position=position_dodge(0.8), width=0.6)  +
  theme(axis.text=element_text(size=18), axis.title.y = element_text(size = 18),legend.position="none") +
  ylab("Percent Identification")


GG1D<- ggarrange(GR2D,GARID,GNcoeffD,
	nrow=1, # align="hv",
	common.legend = TRUE,
	font.label = list(size = 18, color = "black", face = "bold.italic", family = NULL),
	labels=c("a","b","c"))
	
GG2D<-ggarrange(GG1D,GIDD,
	nrow=2,# align="hv",
	#common.legend = TRUE,
	heights=c(1.6,1),
	font.label = list(size = 18, color = "black", face = "bold.italic", family = NULL),
	labels=c("","d"))


GG2D

pdf('Figure2.pdf',width=20,height=12)
(GG2D)
dev.off()


#####################
 ## Long
####################

load("NewEvalLongAll.RData")
Nexp=length(Sim)
Nlong=20


SumPvalAUC<-apply(TPvalAUC,1,CI)
SumAdjPvAUC<-apply(AdjPvAUC,1,CI)
SumLseAUC<-apply(TLseAUC,1,CI)
SumLsminAUC<-apply(TLsminAUC,1,CI)
SumGoldAUC<-apply(GoldAUC,1,CI)
SumSFAUC<-apply(TSOddsAUC,1,CI)
SumSFLassAUC<-apply(SLasminAUC,1,CI)

neach=length(seq(50,2000,by=25))

AucL=data.frame(Time=rep(seq(50,2000,by=25),6), Auc=c(SumGoldAUC[2,],SumLsminAUC[2,],SumAdjPvAUC[2,],SumPvalAUC[2,],SumSFAUC[2,],SumSFLassAUC[2,]),  upper=c(SumGoldAUC[1,],SumLsminAUC[1,],SumAdjPvAUC[1,],SumPvalAUC[1,],SumSFAUC[2,],SumSFLassAUC[2,]), lower=c(SumGoldAUC[3,],SumLsminAUC[3,],SumAdjPvAUC[3,],SumPvalAUC[3,],SumSFAUC[2,],SumSFLassAUC[2,]),  Method=factor(c(rep("Gold",neach),rep("Lasso",neach),rep("Adj-Pval",neach),rep("Seg Pval",neach),rep("Seg Odds",neach), rep("Seg+Lasso",neach)), levels=c("Gold","Lasso","Adj-Pval","Seg Pval","Seg Odds","Seg+Lasso")))


GAucL<-ggplot(AucL, aes(x=Time,y=Auc,group=Method,color=Method))+
		geom_line(size=2, alpha=0.9, show.legend = FALSE) +
		geom_line(aes(y=upper,color=Method,group=Method), size=1, alpha=0.9,linetype="dotted") +
		geom_line(aes(y=lower,color=Method,group=Method), size=1, alpha=0.9,linetype="dotted")  +
		scale_fill_chris() +
		theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold"),
        legend.position="none") +
   ylab("AUC")       +
		ylim(0.5,0.625)

GAucL

GGAUC<-ggarrange(GAucDown,GAucL,
	nrow=1, align="hv",
	labels=c("Prediction downscaled","Prediction long"),
	common.legend = TRUE,
	font.label = list(size = 22))

GGAUC

pdf('Figure4.pdf',width=20,height=12)
(GGAUC)
dev.off()


###########


AllR2L=data.frame(R2=c(GoldR2,  LsminR2, AdjPvR2,PvalR2,SOddsR2,SLasminR2), Method=factor(c(rep("Gold",Nexp),rep("Lasso",Nexp),rep("Adj-Pval",Nexp),rep("Seg Pval",Nexp), rep("Seg Odds",Nlong),rep("Seg+Lasso",Nlong)), levels=c("Gold","Lasso","Adj-Pval", "Seg Pval","Seg Odds","Seg+Lasso")))
AllARIL=data.frame(ARI=c(GoldARI,LsminARI,AdjPvARI,  PvalARI, SOddsARI,SLasminARI), Method=factor(c(rep("Gold",Nexp),rep("Lasso",Nexp),rep("Adj-Pval",Nexp),rep("Seg Pval",Nexp), rep("Seg Odds",Nlong),rep("Seg+Lasso",Nlong)), levels=c("Gold","Lasso","Adj-Pval", "Seg Pval","Seg Odds","Seg+Lasso")))
AllNcoeffL=data.frame(Ncoeff=c(GoldNcoeff, LsminNcoeff, AdjPvNcoeff, PvalNcoeff,SOddsNcoeff,SLasminNcoeff), Method=factor(c(rep("Gold",Nexp),rep("Lasso",Nexp),rep("Adj-Pval",Nexp),rep("Seg Pval",Nexp), rep("Seg Odds",Nlong),rep("Seg+Lasso",Nlong)), levels=c("Gold","Lasso","Adj-Pval", "Seg Pval","Seg Odds","Seg+Lasso")))
AllIDL=data.frame(Id=c(apply(GoldId,1,sum), apply(LsminId,1,sum),  apply(AdjPvId,1,sum),  apply(PvalId,1,sum), 100/Nlong*apply(SOddsId,1,sum), 100/Nlong*apply(SLasminId,1,sum)), Method=factor(c(rep("Gold",11),rep("Lasso",11),rep("Adj-Pval",11), rep("Seg Pval",11),rep("Seg Odds",11),rep("Seg+Lasso",11)), levels=c("Gold","Lasso","Adj-Pval", "Seg Pval","Seg Odds","Seg+Lasso")), Regions=factor(rep(paste("Region",Regions$start,sep=""),6),levels=paste("Region",Regions$start,sep="")))

GR2L<-ggplot(AllR2L, aes(x=Method, y=R2, color=Method)) +
  geom_boxplot()+
		scale_fill_chris() +
theme(axis.text=element_text(size=18),
				legend.key.size = unit(2, 'cm'), legend.text = element_text(size=16), legend.box = "horizontal",
				legend.position="none",
				axis.title.y = element_text(size = 18),
				axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
				ylim(0,40) 


GARIL<-ggplot(AllARIL, aes(x=Method, y=ARI, color=Method)) +
  geom_boxplot()+
		scale_fill_chris() +
theme(axis.text=element_text(size=18),
			axis.title.y = element_text(size = 18),
			legend.position="none",
				#legend.position="bottom", legend.key.size = unit(2, 'cm'), legend.text = element_text(size=16),
				axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


GNcoeffL<-ggplot(AllNcoeffL, aes(x=Method, y=Ncoeff, color=Method)) +
  geom_boxplot()+
		scale_fill_chris() +
	theme(axis.text=element_text(size=18),
				legend.position="none",
				axis.title.y = element_text(size = 18),
				axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
				ylab("Nb of coeff")


GIDL<-ggplot(AllIDL, aes(x=Regions, y=Id, fill=Method)) +
  geom_bar(stat="identity",position=position_dodge(0.8), width=0.6) +
		scale_fill_chris2() +
		 theme(axis.text=element_text(size=18), axis.title.y = element_text(size = 18),legend.position="none") +
		 ylab("Percent Identification")


GG1L<- ggarrange(GR2L,GARIL,GNcoeffL,
	nrow=1, #align="hv",
	common.legend = TRUE,
	font.label = list(size = 18, color = "black", face = "bold.italic", family = NULL),
	labels=c("a","b","c"))
	
GG2L<-ggarrange(GG1L,GIDL,
	nrow=2, #align="hv",
	#common.legend = TRUE,
	heights=c(1.6,1),
	font.label = list(size = 18, color = "black", face = "bold.italic", family = NULL),
	labels=c("","d"))

GG2L


pdf("Images/Figure3.pdf",width=20,height=12)
(GG2L)
dev.off()



