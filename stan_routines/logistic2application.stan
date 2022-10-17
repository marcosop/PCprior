functions{
  //CDF power logistic

  real pplogis(real x, real alpha){             
    real dens = pow(logistic_cdf(x,0,1),alpha);
    return dens;
  } 
  

 // tetragamma function (not available in Rstan)
real tetragamma (real x){ 
  real d = 0.0001;
  int r = 4;
  int v = 2;
  real eps = 0.0001;
  vector[r] a;
  vector[r-1] b3;
  vector[r-2] b2;
  real b1;
  real zerotol = 0.00001781029;
  real h = fabs(d*x) + eps * (fabs(x) < zerotol);
  for(i in 1:r)  { 
      a[i] = (trigamma(x+h) - trigamma(x-h))/(2*h);
         h = h/v;
       }	
  b3 = (a[2:4]*(4)-a[1:3])/(3);
  b2 = (b3[2:3]*(pow(4,2))-b3[1:2])/(pow(4,2)-1);
  b1 = (b2[2]*(pow(4,3))-b2[1])/(pow(4,3)-1);
  return b1;
}
  

    
  
// pdf power logistic 
  real dplogisstan (real x, real mu, real sigma, real alpha){
     real dens = alpha*pow(logistic_cdf(x,mu,sigma),alpha-1)*exp(logistic_lpdf(x|mu,sigma));
     return dens;
}

//  integrand of the theoretical KLD expectation term
real eapp_integrand (real x, real xc, real[] theta, real[] x_r,int[] x_i){
  real alpha = theta[1];
  real ez = digamma(alpha)-digamma(1);
  real vz = trigamma(alpha)+trigamma(1);
  real sigmaalpha = sqrt(1/vz);
  real sigma1 = sqrt(0.5/trigamma(1));
  real mualpha = -sigmaalpha*(digamma(alpha)-digamma(1));
  real integrand = log (1+ exp(-x/sigma1))*dplogisstan(x,mualpha,sigmaalpha,alpha);
  return integrand;
 }
 
//  theoretical KLD expectation term
real eapp_integral(real alpha, data real[] x_r, data int[] x_i){
   real integral=integrate_1d(eapp_integrand, -20,20, {alpha}, x_r, x_i,1e-8);
   return integral;
   }
 
//  theoretical KLD for the power logistic distribution (actually is kind of a hybrid)  
real kldplogis (real alpha, data real[] x_r, data int[] x_i){
  real ez = digamma(alpha)-digamma(1);
  real vz = trigamma(alpha)+trigamma(1);
  real sigmaalpha = sqrt(1/vz);
  real sigma1 = sqrt(0.5/trigamma(1));
  real mualpha = -sigmaalpha*(digamma(alpha)-digamma(1));
  real integral=integrate_1d(eapp_integrand, -20,20, {alpha}, x_r, x_i,1e-8);
  real kldteo = log(alpha*sigma1/sigmaalpha)-((alpha+1)/alpha) + (mualpha/sigmaalpha)+2*integral;
  return kldteo;
}


// derivate of the expectation term used to compute the theoretical KLD.
real deapp(real x, data real[] x_r, data int[] x_i){
  real d = 0.0001;
  int r = 4;
  int v = 2;
  real eps = 0.0001;
  vector[r] a;
  vector[r-1] b3;
  vector[r-2] b2;
  real b1;
  real zerotol = 0.00001781029;
  real h = fabs(d*x) + eps * (fabs(x) < zerotol);
  for(i in 1:r)  { 
      a[i] = (eapp_integral(x+h,x_r,x_i) - eapp_integral(x,x_r,x_i))/(2*h);
         h = h/v;
       }	
  b3 = (a[2:4]*(4)-a[1:3])/(3);
  b2 = (b3[2:3]*(pow(4,2))-b3[1:2])/(pow(4,2)-1);
  b1 = (b2[2]*(pow(4,3))-b2[1])/(pow(4,3)-1);
  return 2*b1;
  }

// derivate of the theoretical KLD.
real dkldplogis(real alpha, data real[] x_r, data int[] x_i){
  real ez = digamma(alpha)-digamma(1);
  real vz = trigamma(alpha)+trigamma(1);
  real sigmaalpha = sqrt(1/(trigamma(alpha)+trigamma(1)));
  real mualpha = -sigmaalpha*(digamma(alpha)-digamma(1));
  real sigma1 = sqrt(0.5/trigamma(1));
  real tetra = tetragamma(alpha);
  real dmualpha = -(sigmaalpha)*(trigamma(alpha)+(sigmaalpha*tetra*mualpha)/2);
  real dsigmaalpha = -0.5*sigmaalpha*tetra/vz;

  real t1 = (1/(alpha*sigmaalpha))*(sigmaalpha-(dsigmaalpha*alpha));
  real t2 = 1/(alpha*alpha);
  real t3 = (1/(sigmaalpha^2))*(sigmaalpha*dmualpha- mualpha*dsigmaalpha);
  real t4 = deapp(alpha,x_r,x_i);
  real dkldteo2 = t1+t2+t3+2*t4;
  return dkldteo2;

}

real pcalpha_lpdf (real alpha,real lambda, data real[] x_r, data int[] x_i){
  real logprior;
 if (alpha>=0.95 && alpha<=1.05){
      logprior = -0.5498479;
   }else if(alpha>=18){
    logprior = -10000;
   }
   else{
     real t11=(0.5*sqrt(2))/sqrt(kldplogis(alpha,x_r,x_i));
  real t1=t11*fabs(dkldplogis(alpha,x_r,x_i));
  real dist=sqrt(2*kldplogis(alpha,x_r,x_i));
    logprior=log(lambda)-lambda*dist+ log(t1);
  }
  return logprior;
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
                      beta[3] * w[start:end] +beta[4] * x1[start:end] +beta[5] * x2[start:end] ,alpha));
                     
  }

  
  ////


}

data {
  int<lower=0> N; // number of observations
  vector[N] x; // independent variable x.
  vector[N] z; // independent variable z;
  vector[N] w; // independent variable w;
  vector[N] x1; // independent variable w;
  vector[N] x2; // independent variable w;
  int<lower=1> grainsize;
  int<lower=0,upper=1> y[N]; // response variable y;
  real x_r[0]; // object of length 0 necessary to use the function integrate1d
  int x_i[0]; // object of length 0 necessary to use the function integrate1d
}
// model parameters specification : in this case q (commented) is the quantile
// used for the idea of Harvard
parameters { 
  real<lower=0,upper= 1> q;
  vector[5] beta;
  real<lower=0> alpha1;
}

transformed parameters{
real beta0 = log(q^(1/alpha1)) - log(1-(q^(1/alpha1)));
}

// model specification
model { 
 real lambda =  2.003899;
  // prior distribution for q (commented)
  q ~ beta(1/alpha1,1);
  alpha1~pcalpha(lambda,x_r,x_i); // pcprior distribution for alpha
  beta[1] ~ normal (0,100); // prior distribution for beta0
  beta[2] ~ normal (0,100); // prior distribution for beta1
  beta[3] ~ normal (0,100); // prior distribution for beta1
  beta[4] ~ normal (0,100); // prior distribution for beta1
  beta[5] ~ normal (0,100); // prior distribution for beta1
   target += reduce_sum(partial_sum_lupmf, y, grainsize, alpha1,
                       x, z,w,x1,x2,beta,beta0);
  
}



