functions{
  real pplogis(real x, real alpha){
    real dens = pow(logistic_cdf(x,0,1),alpha);
    return dens;
  } 
  
  real [] pplogisvec(vector x, real alpha){    
      int N = num_elements(x);
      real  dens[N];
      for( n in 1:N){
     dens[n] = pow(logistic_cdf(x[n],0,1),alpha);
      }
    return dens;
  } 

// parallelizing function
  
 real partial_sum_lpmf(int [] y, int start, int end,
                        real alpha,
                      vector x,
                     vector z,
                      vector w,
                      vector x1,
                      vector x2,
                      vector beta,
                      real beta0) {
                       return bernoulli_lupmf(y |pplogisvec(beta0 + beta[1] * x[start:end] + beta[2] * z[start:end] + 
                      beta[3]*w[start:end] + beta[4]*x1[start:end] + beta[5]*x2[start:end],alpha));
                  
  }
}

data {
  int<lower=0> N;
  vector[N] x;
  vector[N] z;
  vector[N] w;
   vector[N] x1;
  vector[N] x2;
  int<lower=0,upper=1> y[N];
  int<lower=1> grainsize;
  real x_r[0];
  int x_i[0];
}
parameters {
real<lower=0,upper= 1> q;
  vector[5] beta;
  real  logalpha1;
}
transformed parameters{
real<lower=0> alpha1;
alpha1= exp(logalpha1);
real beta0 = log(q^(1/alpha1)) - log(1-(q^(1/alpha1)));
}

model {
  real lambda = 0.05;
  logalpha1 ~ uniform(-2.58938,2.58938);
   q~ beta (1/alpha1,1);
  beta[1] ~ normal (0,100); // prior distribution for beta0
  beta[2] ~ normal (0,100); // prior distribution for beta1
   beta[3] ~ normal (0,100); // prior distribution for beta0
  beta[4] ~ normal (0,100); // prior distribution for beta1
  beta[5] ~ normal (0,100); // prior distribution for beta1
   target += reduce_sum(partial_sum_lupmf, y, grainsize, alpha1,
                       x, z,w,x1,x2, beta,beta0);
}
