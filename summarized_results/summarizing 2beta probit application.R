library(cmdstanr)
library(powdist)
library(tidyverse)
library(HDInterval)
library(patchwork)

###############################################################################
source("auxiliar_functions.R")

###### reading the fitted models ########
dfprobit2 = readRDS("probitapplication.rds")
dfpcpres2 = readRDS("pnormalpcpapplication.rds")
dfunifres2 = readRDS("pnormaluniformapplication.rds")
dfnormres2 = readRDS("pnormalnormalapplication.rds")


dfprobit = dfprobit2%>%filter(.chain == 1)
dfpcpres = dfpcpres2%>%filter(.chain == 1)
dfunifres= dfunifres2%>%filter(.chain == 1)
dfnormres = dfnormres2%>%filter(.chain == 1)

################################################################
####### Computing the Bayesian estimates #######################
################################################################

thetaprobit = cbind(dfprobit$beta0, dfprobit$`beta[1]`, dfprobit$`beta[2]`,dfprobit$`beta[3]`,
                    dfprobit$`beta[4]`,dfprobit$`beta[5]`)
thetapcpres = cbind(dfpcpres$beta0, dfpcpres$`beta[1]`, dfpcpres$`beta[2]`,dfpcpres$`beta[3]`,
                    dfpcpres$`beta[4]`,dfpcpres$`beta[5]`,dfpcpres$alpha1)

thetaunifres = cbind(dfunifres$beta0, dfunifres$`beta[1]`, dfunifres$`beta[2]`,dfunifres$`beta[3]`,
                     dfunifres$`beta[4]`,dfunifres$`beta[5]`,dfunifres$alpha1)

thetanormres = cbind(dfnormres$beta0, dfnormres$`beta[1]`, dfnormres$`beta[2]`,dfnormres$`beta[3]`,
                     dfnormres$`beta[4]`,dfnormres$`beta[5]`,dfnormres$alpha1)

thetaprobitest = c(median(dfprobit$beta0), median(dfprobit$`beta[1]`), median(dfprobit$`beta[2]`),
                   median(dfprobit$`beta[3]`),median(dfprobit$`beta[4]`),median(dfprobit$`beta[5]`))

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

#################################
##### Estimated parameters ######
################################

estimates = cbind(estpcp,estnormal,estuniform)
colnames(estimates) = c("PCP", "Normal", "Uniform")

estimates

intpcp = t(round(apply(thetapcpres,2,FUN=hdi,credMass=0.9),3))
intuniform = t(round(apply(thetaunifres,2,FUN=hdi,credMass=0.9),3))
intnormal = t(round(apply(thetanormres,2,FUN=hdi,credMass=0.9),3))

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

terms=-2*apply(as.matrix(thetapcpres),1,FUN=likprobitpcp,matx=x,y=y)
logdens=apply(as.matrix(thetapcpres),1,FUN=termiprobit,matx=x,y=y)
lppd=sum(apply(logdens,1,mean))
varlog=sum(apply(logdens,1,var))
Dhat = -2*likprobitpcp(as.numeric(thetaunifresest),matx=x,y=y)

p=length(thetapcpres[1,])
DICpcp=2*mean(terms)- Dhat
BICpcp=mean(terms)+p*log(n)
AICpcp=mean(terms)+2*p
WAICpcp= Dhat + 2*varlog

### CRITERIA UNIFORM

terms1=-2*apply(as.matrix(thetaunifres),1,FUN=likprobitpcp,matx=x,y=y)
logdens1=apply(as.matrix(thetaunifres),1,FUN=termiprobit,matx=x,y=y)
lppd1=sum(apply(logdens,1,mean))
varlog1=sum(apply(logdens1,1,var))
Dhat  = - 2*likprobitpcp(as.numeric(thetaunifresest),matx=x,y=y)

p=length(thetaunifres[1,])
DICpriorunif=2*mean(terms1) - Dhat
BICpriorunif=mean(terms1)+p*log(n)
AICpriorunif=mean(terms1)+2*p
WAICpriorunif=Dhat + 2*varlog1


### CRITERIA NORMAL

terms1=-2*apply(as.matrix(thetanormres),1,FUN=likprobitpcp,matx=x,y=y)
logdens1=apply(as.matrix(thetanormres),1,FUN=termiprobit,matx=x,y=y)
lppd1=sum(apply(logdens,1,mean))
varlog1=sum(apply(logdens1,1,var))
Dhat = - 2*likprobitpcp(as.numeric(thetanormresest),matx=x,y=y)

p=length(thetanormres[1,])
DICpriornorm=2*mean(terms1) - Dhat
BICpriornorm=mean(terms1)+p*log(n)
AICpriornorm=mean(terms1)+2*p
WAICpriornorm=Dhat + 2*varlog



### CRITERIA TRADITIONAL PROBIT #############


terms1=-2*apply(as.matrix(thetaprobit),1,FUN=likprobit,matx=x,y=y)
logdens1=apply(as.matrix(thetaprobit),1,FUN=termiprobit1,matx=x,y=y)
lppd1=sum(log(apply(exp(logdens1),1,mean)))
varlog1=sum(apply(logdens1,1,var))

p=length(thetaprobit[1,])
DICprobit=2*mean(terms1) +2*likprobit(as.numeric(thetaprobitest),matx=x,y=y)
BICprobit=mean(terms1)+p*log(n)
AICprobit=mean(terms1)+2*p
WAICprobit=-2*(lppd1-varlog1)

####summary criteria###3
DIC=c(DICpcp,DICpriornorm,DICpriorunif,DICprobit)
AIC=c(AICpcp,AICpriornorm,AICpriorunif,AICprobit)
BIC=c(BICpcp,BICpriornorm,BICpriorunif,BICprobit)
WAIC=c(WAICpcp,WAICpriornorm,WAICpriorunif,WAICprobit)

criteria = round(rbind(DIC,AIC,BIC,WAIC),3)

colnames(criteria) = c("PCP", "Normal", "Uniform", "Trad. Probito")

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

pnormnormad=apply(thetanormres,1,FUN=adjpnorm,matx=xtest)
pnormunifad=apply(thetaunifres,1,FUN=adjpnorm,matx=xtest)
pnormpcpad=apply(thetapcpres,1,FUN=adjpnorm,matx=xtest)
probitoad=apply(thetaprobit,1,FUN=adjprobit,matx=xtest)

pnormad = apply(pnormnormad,1,median)
punifad = apply(pnormunifad,1,median)
ppcpad = apply(pnormpcpad,1,median)
probitad = apply(probitoad,1,median)


library(pROC)


##################################
#### PREDICTION CRITERIA##########
##################################

#### AUC values
auxnorm=roc(predictor=pnormad,response=ytest)
auxunif=roc(predictor=punifad,response=ytest)
auxpcp=roc(predictor=ppcpad,response=ytest)
auxprobit=roc(predictor=probitad,response=ytest)

###### different criteria values (ACC, TPR, TNR, CSI, SSI, FAITH, PDIF)
claslogis(ppcpad,0.329,ytest)
claslogis(probitad,0.329,ytest)











