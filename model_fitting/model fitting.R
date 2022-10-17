library(foreach)
library(doParallel)
library(ranger);library(palmerpenguins);library(kableExtra)
library(powdist)
library(cmdstanr)
library(tidyverse)
library(powdist)

dados = read.table("dataset/coverageX.txt",header=T)

### declaring models #######

logisticpcp <- cmdstan_model("stan_routines/logistic2application.stan", cpp_options = list(stan_threads = TRUE))
logisticuniform <- cmdstan_model("stan_routines/uniform priorapplication.stan", cpp_options = list(stan_threads = TRUE))
logisticnormal <- cmdstan_model("stan_routines/normal priorapplication.stan", cpp_options = list(stan_threads = TRUE))

pnormalpcp <- cmdstan_model("stan_routines/pnormalpcpapplication.stan", cpp_options = list(stan_threads = TRUE))
pnormaluniform <- cmdstan_model("stan_routines/normal ppnormapplication.stan", cpp_options = list(stan_threads = TRUE))
pnormalnormal <- cmdstan_model("stan_routines/uniform ppnormapplication.stan", cpp_options = list(stan_threads = TRUE))

### usual links #### 
probit <- cmdstan_model("stan_routines/probitalpha1.stan", cpp_options = list(stan_threads = TRUE))
logit <- cmdstan_model("stan_routines/logitalpha1application.stan", cpp_options = list(stan_threads = TRUE))


n = length(dados$y)
p1 = mean(dados$y)
p0 = 1-p1
set.seed(1234)

####Getting the test dataset
pos1 = sample(which(dados$y ==1),round(0.1*p1*n))
pos0 = sample(which(dados$y ==0),round(0.1*p0*n))

## Train dataset 
dados = dados[-c(pos0,pos1),]

y = dados$y
x= dados$MEN
z= dados$URBAN
w = dados$PRIVATE
x1 = dados$AGE
x2 = dados$SENIORITY

my_data=list(N=length(x),x=x,y=y,z=z,w=w,x1 = x1, x2=x2, x_r=numeric(0),x_i=numeric(0),grainsize = 1)

#### Number of iterations#####
iter2 = 10000



#### traditional links
logitstan = logit$sample(my_data,
                         chains = 2,
                         parallel_chains = 2,
                         threads_per_chain = 10,
                         refresh = 1000,iter_warmup =0.1*iter2,
                         iter_sampling = 0.9*iter2,
                         init = function() list(beta0 = 0.1,beta = c(0.1,0.1)),
                         adapt_delta = 0.8,thin = 5) 




probitstan = probit$sample(my_data,
                           chains = 2,
                           parallel_chains = 2,
                           threads_per_chain = 10,
                           refresh = 1000,iter_warmup =0.1*iter2,
                           iter_sampling = 0.9*iter2,
                           init = function() list(beta0 = 0.1,beta = c(0.1,0.1)),
                           adapt_delta = 0.8,thin = 5) 




##### uniform prior ############

logisticuniformstan = logisticuniform$sample(my_data,
                                     chains = 2,
                                     parallel_chains = 2,
                                     threads_per_chain = 10,
                                     refresh = 1000,iter_warmup =0.1*iter2,
                                     iter_sampling = 0.9*iter2,
                                     init = function() list(logalpha1= -0.6931472,
                                                            q = 0.5,beta = c(0.1,0.1)),
                                     adapt_delta = 0.8,thin = 5) 


pnormaluniformstan = pnormaluniform$sample(my_data,
                                           chains = 2,
                                           parallel_chains = 2,
                                           threads_per_chain = 10,
                                           refresh = 1000,iter_warmup =0.1*iter2,
                                           iter_sampling = 0.9*iter2,
                                           init = function() list(logalpha1= -0.6931472,
                                                                  q = 0.5,beta = c(0.1,0.1)),
                                           adapt_delta = 0.8,thin = 5) 



##### normal prior ######

logisticnormalstan = logisticnormal$sample(my_data,
                                             chains = 2,
                                             parallel_chains = 2,
                                             threads_per_chain = 10,
                                             refresh = 1000,iter_warmup =0.1*iter2,
                                             iter_sampling = 0.9*iter2,
                                             init = function() list(alpha1 = 0.4,
                                                                    q = 0.5,beta = c(0.1,0.1)),
                                             adapt_delta = 0.8,thin = 5) 





pnormalnormalstan = pnormalnormal$sample(my_data,
                                           chains = 2,
                                           parallel_chains = 2,
                                           threads_per_chain = 10,
                                           refresh = 1000,iter_warmup =0.1*iter2,
                                           iter_sampling = 0.9*iter2,
                                           init = function() list(alpha1 = 0.4,
                                                                  q = 0.5,beta = c(0.1,0.1)),
                                           adapt_delta = 0.8,thin = 5) 



#### pcp prior###############

pnormalpcpstan = pnormalpcp$sample(my_data,
                                   chains = 2,
                                   parallel_chains = 2,
                                   threads_per_chain = 10,
                                   refresh = 1000,iter_warmup =0.1*iter2,
                                   iter_sampling = 0.9*iter2,
                                   init = function() list(alpha1 = 0.4,
                                                          q = 0.5,beta = c(0.1,0.1)),
                                   adapt_delta = 0.8,thin = 5) 





logisticpcpstan = logisticpcp$sample(my_data,
                                     chains = 2,
                                     parallel_chains = 2,
                                     threads_per_chain = 10,
                                     refresh = 1000,iter_warmup =0.1*iter2,
                                     iter_sampling = 0.9*iter2,
                                     init = function() list(alpha1 = 0.4,
                                                            q = 0.5,beta = c(0.1,0.1)),
                                     adapt_delta = 0.8,thin = 5) 





