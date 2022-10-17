library(cmdstanr)
library(powdist)
library(HDInterval)
library(tidyverse)
library(ggpubr)
library(patchwork)


###############################################################################
source("auxiliar_functions.R")

###### reading the fitted models ########
dflogit2 = readRDS("logitapplication.rds")
dfpcpres2 = readRDS("logisticpcpcoverage.rds")
dfunifres2 = readRDS("logisticuniformcoverage.rds")
dfnormres2 = readRDS("logisticnormalcoverage.rds")

dflogit = dflogit2%>%filter(.chain == 1)
dfpcpres = dfpcpres2%>%filter(.chain == 1)
dfunifres= dfunifres2%>%filter(.chain == 1)
dfnormres = dfnormres2%>%filter(.chain == 1)


################################################################
####### Computing the Bayesian estimates #######################
################################################################

thetalogit = cbind(dflogit$beta0, dflogit$`beta[1]`, dflogit$`beta[2]`,dflogit$`beta[3]`,
                    dflogit$`beta[4]`,dflogit$`beta[5]`)
thetapcpres = cbind(dfpcpres$beta0, dfpcpres$`beta[1]`, dfpcpres$`beta[2]`,dfpcpres$`beta[3]`,
                    dfpcpres$`beta[4]`,dfpcpres$`beta[5]`,dfpcpres$alpha1)

thetaunifres = cbind(dfunifres$beta0, dfunifres$`beta[1]`, dfunifres$`beta[2]`,dfunifres$`beta[3]`,
                     dfunifres$`beta[4]`,dfunifres$`beta[5]`,dfunifres$alpha1)

thetanormres = cbind(dfnormres$beta0, dfnormres$`beta[1]`, dfnormres$`beta[2]`,dfnormres$`beta[3]`,
                     dfnormres$`beta[4]`,dfnormres$`beta[5]`,dfnormres$alpha1)

thetalogitest = c(median(dflogit$beta0), median(dflogit$`beta[1]`), median(dflogit$`beta[2]`),
                   median(dflogit$`beta[3]`),median(dflogit$`beta[4]`),median(dflogit$`beta[5]`))

thetapcpresest = c(median(dfpcpres$beta0), median(dfpcpres$`beta[1]`), median(dfpcpres$`beta[2]`),
                   median(dfpcpres$`beta[3]`),median(dfpcpres$`beta[4]`),
                   median(dfpcpres$`beta[5]`),median(dfpcpres$alpha1))

thetaunifresest = c(median(dfunifres$beta0), median(dfunifres$`beta[1]`), median(dfunifres$`beta[2]`),
                    median(dfunifres$`beta[3]`),median(dfunifres$`beta[4]`),
                    median(dfunifres$`beta[5]`),median(dfunifres$alpha1))

thetanormresest = c(median(dfnormres$beta0), median(dfnormres$`beta[1]`), median(dfnormres$`beta[2]`),
                    median(dfnormres$`beta[3]`),median(dfnormres$`beta[4]`),
                    median(dfnormres$`beta[5]`),median(dfnormres$alpha1))



estpcp= t(rbind(round(thetapcpresest,3)))
estuniform =  t(rbind(round(thetaunifresest,3)))
estnormal = t(rbind(round(thetanormresest,3)))
estlogit =  t(rbind(round(thetalogitest,3)))

#################################
##### Estimated parameters ######
################################

estimates = cbind(estpcp,estnormal,estuniform)
colnames(estimates) = c("PCP", "Normal", "Uniform")

estimates

intpcp = t(round(apply(thetapcpres,2,FUN=hdi,credMass=0.95),3))
intuniform = t(round(apply(thetaunifres,2,FUN=hdi,credMass=0.95),3))
intnormal = t(round(apply(thetanormres,2,FUN=hdi,credMass=0.95),3))
intlogit = t(round(apply(thetalogit,2,FUN=hdi,credMass=0.95),3))

#########################################
######Credibility intervals #############
#########################################

intervals = cbind(intpcp,intnormal,intuniform)
colnames(intervals) = c("lowerPCP", "upperPCP","lowerNormal", "upperNormal",
                        "lowerUniform","upperUniform")
intervals
####################################################################################
######## Real data: Computing model selection criteria and prediction measures #####
####################################################################################

dados = read.table("coverageX.txt",header=T)

n = length(dados$y)
p1 = mean(dados$y)
p0 = 1-p1
set.seed(1234)
pos1 = sample(which(dados$y ==1),round(0.1*p1*n))
pos0 = sample(which(dados$y ==0),round(0.1*p0*n))

#### I am gonna need this latter for the diagnostic plots ##############
dadostot = dados
ytot = dadostot$y
x1tot= dadostot$MEN
x2tot= dadostot$URBAN
x3tot = dadostot$PRIVATE
x4tot = dadostot$AGE
x5tot = dadostot$SENIORITY

xtot = cbind(1,x1tot,x2tot,x3tot,x4tot,x5tot)

##############################################

dados = dadostot[-c(pos0,pos1),]
dadostest = dadostot[c(pos0,pos1),]

y = dados$y
x1= dados$MEN
x2= dados$URBAN
x3 = dados$PRIVATE
x4 = dados$AGE
x5 = dados$SENIORITY




x = cbind(1,x1,x2,x3,x4,x5)
y = y

###################################################################333

###############################################################################
########## Computing model selection criteria #################################
###############################################################################

##### CRITERIA PC PRIOR

terms=-2*apply(as.matrix(thetapcpres),1,FUN=liklogitpcp,matx=x,y=y)
logdens=apply(as.matrix(thetapcpres),1,FUN=termilogit,matx=x,y=y)
lppd=sum(apply(logdens,1,mean))
varlog=sum(apply(logdens,1,var))
Dhat = -2*liklogitpcp(as.numeric(thetaunifresest),matx=x,y=y)


p=length(thetapcpres[1,])
DIClogis=2*mean(terms)- Dhat
BIClogis=mean(terms)+p*log(n)
AIClogis=mean(terms)+2*p
WAIClogis= Dhat + 2*varlog

### CRITERIA UNIFORM

terms1=-2*apply(as.matrix(thetaunifres),1,FUN=liklogitpcp,matx=x,y=y)
logdens1=apply(as.matrix(thetaunifres),1,FUN=termilogit,matx=x,y=y)
lppd1=sum(apply(logdens,1,mean))
varlog1=sum(apply(logdens1,1,var))
Dhat  = - 2*liklogitpcp(as.numeric(thetaunifresest),matx=x,y=y)


p=length(thetaunifres[1,])
DICpriorunif=2*mean(terms1) - Dhat
BICpriorunif=mean(terms1)+p*log(n)
AICpriorunif=mean(terms1)+2*p
WAICpriorunif=Dhat + 2*varlog1

### CRITERIA NORMAL

terms1=-2*apply(as.matrix(thetanormres),1,FUN=liklogitpcp,matx=x,y=y)
logdens1=apply(as.matrix(thetanormres),1,FUN=termilogit,matx=x,y=y)
lppd1=sum(apply(logdens,1,mean))
varlog1=sum(apply(logdens1,1,var))
Dhat = - 2*liklogitpcp(as.numeric(thetanormresest),matx=x,y=y)

p=length(thetanormres[1,])
DICpriornorm=2*mean(terms1) - Dhat
BICpriornorm=mean(terms1)+p*log(n)
AICpriornorm=mean(terms1)+2*p
WAICpriornorm=Dhat + 2*varlog


### CRITERIA TRADITIONAL LOGISTIC #############

terms1=-2*apply(as.matrix(thetalogit),1,FUN=liklogit,matx=x,y=y)
logdens1=apply(as.matrix(thetalogit),1,FUN=termilogit1,matx=x,y=y)
lppd1=sum(log(apply(exp(logdens1),1,mean)))
varlog1=sum(apply(logdens1,1,var))


p=length(thetalogit[1,])
DIClogit=2*mean(terms1) +2*liklogit(as.numeric(thetalogitest),matx=x,y=y)
BIClogit=mean(terms1)+p*log(n)
AIClogit=mean(terms1)+2*p
WAIClogit=-2*(lppd1-varlog1)

####summary criteria###3
DIC=c(DIClogis,DICpriornorm,DICpriorunif,DIClogit)
AIC=c(AIClogis,AICpriornorm,AICpriorunif,AIClogit)
BIC=c(BIClogis,BICpriornorm,BICpriorunif,BIClogit)
WAIC=c(WAIClogis,WAICpriornorm,WAICpriorunif,WAIClogit)

### COMPUTED CRITERIA ####

criteria = round(rbind(DIC,AIC,BIC,WAIC),3)

colnames(criteria) = c("PCP", "Normal", "Uniform", "Trad. Logit")

criteria


##### getting the test data
ytest = dadostest$y
x1test = dadostest$MEN
x2test = dadostest$URBAN
x3test = dadostest$PRIVATE
x4test = dadostest$AGE
x5test = dadostest$SENIORITY

xtest = cbind(1,x1test,x2test,x3test,x4test,x5test)
y = ytest

############################

plogisnormad=apply(thetanormres,1,FUN=adjplogis,matx=xtest)
plogisunifad=apply(thetaunifres,1,FUN=adjplogis,matx=xtest)
plogispcpad=apply(thetapcpres,1,FUN=adjplogis,matx=xtest)
logisticad=apply(thetalogit,1,FUN=adjlogis,matx=xtest)

pnormad = apply(plogisnormad,1,mean)
punifad = apply(plogisunifad,1,mean)
ppcpad = apply(plogispcpad,1,mean)
logisad = apply(logisticad,1,mean)


library(pROC)

##################################
#### PREDICTION CRITERIA##########
##################################

#### AUC values
auxnorm=roc(predictor=pnormad,response=ytest)
auxunif=roc(predictor=punifad,response=ytest)
auxpcp=roc(predictor=ppcpad,response=ytest)
auxlogis=roc(predictor=logisad,response=ytest)

###### different criteria values (ACC, TPR, TNR, CSI, SSI, FAITH, PDIF)
claslogis(ppcpad,0.329,ytest)
claslogis(logisad,0.329,ytest)


#######################################
#### independent chains Plots ########
######################################

labelsa = c("Chain 1","Chain 2")
dfpcpres2$iter = rep(1:1800,2)

 ##### plots
 library(tidyverse)
 ### pc prior
 g0 = dfpcpres2%>%ggplot(aes(x=iter,y = beta0,color=factor(.chain))) + geom_line() +
   theme_bw()+
   scale_color_manual(name = "Chains", values  = rep(c("black","gray")), labels = labelsa)+
   ylab(expression(beta[0])) + xlab("Index")

 g1 = dfpcpres2%>%ggplot(aes(x=iter,y = `beta[1]`,color=factor(.chain))) + geom_line() +
   theme_bw()+
   scale_color_manual(name = "Chains", values  = rep(c("black","gray")), labels = labelsa)+
   ylab(expression(beta[1])) + xlab("Index")

 g2 = dfpcpres2%>%ggplot(aes(x=iter,y = `beta[2]`,color=factor(.chain))) + geom_line() +
   theme_bw()+
   scale_color_manual(name = "Chains", values  = rep(c("black","gray")), labels = labelsa)+
   ylab(expression(beta[2])) + xlab("Index")

 g3 = dfpcpres2%>%ggplot(aes(x=iter,y = `beta[3]`,color=factor(.chain))) + geom_line() +
   theme_bw()+
   scale_color_manual(name = "Chains", values  = rep(c("black","gray")), labels = labelsa)+
   ylab(expression(beta[3])) + xlab("Index")

 g4 = dfpcpres2%>%ggplot(aes(x=iter,y = `beta[4]`,color=factor(.chain))) + geom_line() +
   theme_bw()+
   scale_color_manual(name = "Chains", values  = rep(c("black","gray")), labels = labelsa)+
   ylab(expression(beta[4])) + xlab("Index")

 g5 = dfpcpres2%>%ggplot(aes(x=iter,y = `beta[5]`,color=factor(.chain))) + geom_line() +
   theme_bw()+
   scale_color_manual(name = "Chains", values  = rep(c("black","gray")), labels = labelsa)+
   ylab(expression(beta[5])) + xlab("Index")


 g6 = dfpcpres2%>%ggplot(aes(x=iter,y = `alpha1`,color=factor(.chain))) + geom_line() +
   theme_bw()+
   scale_color_manual(name = "Chains", values  = rep(c("black","gray")), labels = labelsa)+
   ylab(expression(alpha)) + xlab("Index")



 g0 + g1+g2 +g3 +g4 + g5 +g6 + plot_layout(guides = 'collect')



 ggsave("/home/ordonezjosealejandro/Escritorio/resultadoscan/chainsbetasapp2.eps",width = 10, height =7, device=cairo_ps)

#################################
#### Density plot ##############
################################
alphapc =cbind(dfpcpres2$alpha1, "PCP")
alphanormal = cbind(dfnormres2$alpha1,"Normal")
alphaunif = cbind(dfunifres2$alpha1,"Uniform")
alphas = rbind(alphapc,alphanormal,alphaunif)
colnames(alphas) = c("alpha", "Prior")
alphas = data.frame(alphas)
alphas$alpha = as.numeric(alphas$alpha)

## density plot
alphas%>%ggplot(aes(x=alpha,fill=Prior)) + geom_density(alpha=0.4) + scale_fill_manual(values=c("gray60", "gray1", "white"))+
  theme_bw() + coord_cartesian(x=c(2,20),y=c(0,0.3))+
  labs(title="",y="Density",x = expression(alpha))

ggsave("/home/ordonezjosealejandro/Escritorio/resultadoscan/posteriorappalpha.eps",width = 6, height =4, device=cairo_ps)

###############
#### ACF plots
##############
 library(ggfortify)

 p1 <- autoplot(acf(dfpcpres2$beta0, plot = FALSE,lag.min=0,lag.max=30), conf.int.value = 0.95
                ,conf.int.fill=1,,conf.int.col=1)
 p1 = p1+geom_hline(yintercept = 0)+ coord_cartesian(x=c(0,28.5))+
   ggtitle(expression(beta[0]))+
   theme(plot.title = element_text(hjust = 0.5)) + theme_bw()

 p2 <- autoplot(acf(dfpcpres2$`beta[1]`, plot = FALSE,lag.min=0,lag.max=30), conf.int.value = 0.95
                ,conf.int.fill=1,,conf.int.col=1)
 p2 = p2+geom_hline(yintercept = 0)+ coord_cartesian(x=c(0,28.5))+
   ggtitle(expression(beta[1]))+
   theme(plot.title = element_text(hjust = 0.5)) + theme_bw()


 p3 <- autoplot(acf(dfpcpres2$`beta[2]`, plot = FALSE,lag.min=0,lag.max=30), conf.int.value = 0.95
                ,conf.int.fill=1,,conf.int.col=1)
 p3 = p3+geom_hline(yintercept = 0)+ coord_cartesian(x=c(0,28.5))+
   ggtitle(expression(beta[2]))+
   theme(plot.title = element_text(hjust = 0.5)) + theme_bw()

 p4 <- autoplot(acf(dfpcpres2$`beta[3]`, plot = FALSE,lag.min=0,lag.max=30), conf.int.value = 0.95
                ,conf.int.fill=1,,conf.int.col=1)
 p4 = p4+geom_hline(yintercept = 0)+ coord_cartesian(x=c(0,28.5))+
   ggtitle(expression(beta[2]))+
   theme(plot.title = element_text(hjust = 0.5)) + theme_bw()

 p5 <- autoplot(acf(dfpcpres2$`beta[4]`, plot = FALSE,lag.min=0,lag.max=30), conf.int.value = 0.95
                ,conf.int.fill=1,,conf.int.col=1)
 p5 = p5+geom_hline(yintercept = 0)+ coord_cartesian(x=c(0,28.5))+
   ggtitle(expression(beta[2]))+
   theme(plot.title = element_text(hjust = 0.5)) + theme_bw()

 p6 <- autoplot(acf(dfpcpres2$`beta[5]`, plot = FALSE,lag.min=0,lag.max=30), conf.int.value = 0.95
                ,conf.int.fill=1,,conf.int.col=1)
 p6 = p6+geom_hline(yintercept = 0)+ coord_cartesian(x=c(0,28.5))+
   ggtitle(expression(beta[2]))+
   theme(plot.title = element_text(hjust = 0.5)) + theme_bw()



 p7 <- autoplot(acf(dfpcpres2$`alpha1`, plot = FALSE,lag.min=0,lag.max=30), conf.int.value = 0.95
                ,conf.int.fill=1,,conf.int.col=1)
 p7 = p7+geom_hline(yintercept = 0)+ coord_cartesian(x=c(0,28.5))+
   ggtitle(expression(alpha))+
   theme(plot.title = element_text(hjust = 0.5)) + theme_bw()

 ### ACF PLOT ###
 p1+p2 +p3 +p4 + p5 + p6 + p7 + plot_layout(guides = 'collect')

 ggsave("/home/ordonezjosealejandro/Escritorio/resultadoscan/acfbetasapp2.eps",width = 10, height =7, device=cairo_ps)

####################################
#### DIAGNOTSTIC PLOTS############
##################################
plogispcpadtot=apply(thetapcpres,1,FUN=adjplogis,matx=xtot)

ppcpadtot = apply(plogispcpadtot,1,mean)

matpcpnorm=cbind(ytot,ppcpadtot)

respcpnorm=apply(matpcpnorm,1,FUN=quantres)

library(ggpubr)
a1 = ggqqplot(respcpnorm,ggtheme=theme_gray(),xlab="Theoretical Quantiles",ylab="Quantile residuals")+
   theme(axis.text=element_text(size=13),axis.title=element_text(size=13))+
   theme_bw()

respcp3 = data.frame(1:length(respcpnorm),respcpnorm)
names(respcp3) = c("index","residuals")
a2 = ggplot(data=respcp3, aes(residuals)) +
  geom_histogram(col=1,fill="gray",aes(y = stat(count) / sum(count))) +
   coord_cartesian(y=c(0,0.11))+
   labs(title="",x="Quantile residuals",y="Frequency")+
   theme(axis.text=element_text(size=13),axis.title=element_text(size=13))+
   theme_bw()

a3 = ggplot(data=respcp3, aes(y=residuals,x=index)) + geom_point()+
   coord_cartesian(y=c(-4,4))+
   labs(title="",x="index",y="Quantile residuals")+
   theme(axis.text=element_text(size=13),axis.title=element_text(size=13))+
   theme_bw()

### DIAGNOSTIC PLOT ####
 ((a2 / a3 + plot_layout(guides = 'keep')) | a1)


 ggsave("/home/ordonezjosealejandro/Escritorio/resultadoscan/residualsapp2.eps",width = 6, height =3, device=cairo_ps)

