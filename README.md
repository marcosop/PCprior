PENALIZED COMPLEXITY PRIORS FOR THE SKEWNESS PARAMETER OF POWER LINKS
MAIN AUTHOR: JOSE ALEJANDRO ORDONEZ
email: ordonezjosealejandro@gmail.com
University: State University of Campinas.


These are the instructions for running the codes attached to the folder. 

The scripts are organized in four internal folders:

*dataset
*model_fitting 
*stan_routines
*summarized_results


dataset: Includes the real dataset used for the application section of the manuscript. The data consist of a 
sample of 4000 motor insurance policyholders randomly selected by gender. The response variable is whether 
the client buys a full coverage plan (Y = 1) or not (Y = 0). For more details, see the main manuscript.


model_fitting: Contain the script used for running the models. The sample obtained for each model were saved as .rds archives, and their
names are listed below, for each link and each prior 

*** Generalized logistic link***

** 'logisticpcpcoverage.rds' contains the sample obtained under the penalized complexity prior for the skewness alpha.
** 'logisticuniformcoverage.rds' contains the sample obtained under the uniform prior for the skewness alpha.
** 'logisticnormalcoverage.rds' contains the sample obtained under the uniform prior for the skewness alpha.

*** Power normal link***

** 'pnormalpcpapplication.rds' contains the sample obtained under the penalized complexity prior for the skewness alpha.
** 'pnormaluniformapplication.rds' contains the sample obtained under the uniform prior for the skewness alpha.
** 'pnormalnormalapplication.rds' contains the sample obtained under the uniform prior for the skewness alpha.


*** Traditional links ***

** 'logitapplication.rds' contains the sample obtained for the logit link.
** 'probitation.rds' contains the sample obtained for the probit link.



stan routines: contains the stan routines used for the Bayesian model fitting for the generalized logistic,
power normal and traditional links under different priors. Its contents are listed below.

*** Generalized logistic link *****

** 'logistic2application.stan' was the stan routine used for running this model under the pc prior.
** 'normal priorapplication.stan' was the stan routine used for running this model under the normal prior.
** 'uniform priorapplication.stan' was the stan routine used for running this model under the uniform prior.


*** Power normal link *****

** 'pnormalpcpapplication.stan' was the stan routine used for running this model under the pc prior.
** 'normal ppnormapplication.stan' was the stan routine used for running this model under the normal prior.
** 'uniform ppnormapplication.stan' was the stan routine used for running this model under the uniform

*** Traditional links ***

** 'logitalpha1application.stan' was the stan routine used for running the traditional logistic model under the pc prior.
** 'probitalpha1.stan' was the stan routine used for running the traditional probit model under the pc prior.




summarized results: contains the .rds archives as well as the scripts used for processing 
the saved samples (.rds archives), its contents are listed below.

** 'summarizing 2beta logit application.R' was the script used for processing the results for the logistic
	type links under the normal, uniform and pc prior.

** 'summarizing 2beta probit application.R' was the script used for processing the results for the normal
	type links under the normal, uniform and pc prior.

** .rds archives generated during the model fitting stage (described above).

** auxiliar functions.R: Auxiliar functions used for the computation of different model selection criteria.








