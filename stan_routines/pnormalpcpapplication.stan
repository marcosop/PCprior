functions{
   //CDF power normal
  real Phialpha(real x, real alpha){
    real dens = pow (Phi(x),alpha);
    return dens;
  } 
// pdf power normal 
  real dpnormstan (real x, real mu, real sigma, real alpha){
    real dens = alpha*pow(normal_cdf(x,mu,sigma),alpha-1)*exp(normal_lpdf(x|mu,sigma));
    return dens;
  } 

// finding the two first moments of the power normal (0,1)//
  
  real intg1alpha(real x, real xc, real[] theta, real[] x_r,int[] x_i){ 
    real alpha = theta[1];
    real intg = x*dpnormstan(x,0,1,alpha);
    return intg; 
  }

  real intg2alpha(real x, real xc, real[] theta, real[] x_r,int[] x_i){ 
    real alpha = theta[1];
    real intg = pow(x,2)*dpnormstan(x,0,1,alpha);
    return intg; 
  }

 //  theoretical momments
  real ezpnorm(real alpha, data real[] x_r, data int[] x_i){
    real integral=integrate_1d(intg1alpha, -20,20, {alpha}, x_r, x_i,1e-8);
    return integral;
   }

  real ezzpnorm(real alpha, data real[] x_r, data int[] x_i){
    real integral=integrate_1d(intg2alpha, -20,20, {alpha}, x_r, x_i,1e-8);
    return integral;
   }



//  theoretical KLD for the power logistic distribution (actually is kind of a hybrid)  
real kldpnorm (real alpha, data real[] x_r, data int[] x_i){
  real ez = ezpnorm(alpha,x_r,x_i);
  real ezz = ezzpnorm(alpha,x_r,x_i); 
  real vz = ezz- pow(ez,2);
  real sigmaalpha = sqrt(1/vz);
  real mualpha = -sigmaalpha*ez;
  real kldteo = log(alpha)-((alpha-1)/alpha) + 0.5 - log(sigmaalpha)- ((1+pow(mualpha,2))/(2*pow(sigmaalpha,2)));
  return kldteo;
}


// derivates of ezpnorm and ezzpnorm

real dezpnorm(real x, data real[] x_r, data int[] x_i){
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
      a[i] = (ezpnorm(x+h,x_r,x_i) - ezpnorm(x,x_r,x_i))/(2*h);
         h = h/v;
       }	
  b3 = (a[2:4]*(4)-a[1:3])/(3);
  b2 = (b3[2:3]*(pow(4,2))-b3[1:2])/(pow(4,2)-1);
  b1 = (b2[2]*(pow(4,3))-b2[1])/(pow(4,3)-1);
  return 2*b1;
  }


real dezzpnorm(real x, data real[] x_r, data int[] x_i){
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
      a[i] = (ezzpnorm(x+h,x_r,x_i) - ezzpnorm(x,x_r,x_i))/(2*h);
         h = h/v;
       }	
  b3 = (a[2:4]*(4)-a[1:3])/(3);
  b2 = (b3[2:3]*(pow(4,2))-b3[1:2])/(pow(4,2)-1);
  b1 = (b2[2]*(pow(4,3))-b2[1])/(pow(4,3)-1);
  return 2*b1;
  }


// derivate of the theoretical KLD.
real dkldpnorm(real alpha, data real[] x_r, data int[] x_i){
  real ez = ezpnorm(alpha,x_r,x_i);
  real ezz = ezzpnorm(alpha,x_r,x_i); 
  real vz = ezz- pow(ez,2);
  real dvz = dezzpnorm(alpha,x_r,x_i)- (2*ezpnorm(alpha,x_r,x_i)*dezpnorm(alpha,x_r,x_i));
  real sigmaalpha = sqrt(1/vz);
  real mualpha = -sigmaalpha*ez;
  real dmualpha = ((0.5*dvz*ez/sqrt(vz))- (dezpnorm(alpha,x_r,x_i)*sqrt(vz)))/vz;
  real dsigmaalpha =-0.5*pow(vz,-1.5)*dvz;
  
  real t1= (alpha-1)/(alpha^2);
  real t2=(1/sigmaalpha)*dsigmaalpha;
  real t3=((mualpha^2)+1)/(sigmaalpha^2)-1;
  real t4=mualpha*dmualpha/(sigmaalpha^2);
  real derivhybteo=t1+(t2*t3)-t4;
  return derivhybteo;
}

// logarithm of the penalised complexity prior for alpha using the power logistic distribution.
real pcalphapnorm_lpdf (real alpha,real lambda, data real[] x_r, data int[] x_i){
  real logprior;
   if (alpha>=0.95 && alpha<=1.05){
      logprior = -1.750838;
   }else if (alpha>=16){
    logprior = -10000;
   }else{
  real t11=(0.5*sqrt(2))/sqrt(kldpnorm(alpha,x_r,x_i));
  real t1=t11*fabs(dkldpnorm(alpha,x_r,x_i));
  real dist=sqrt(2*kldpnorm(alpha,x_r,x_i));
 
    logprior=log(lambda)-lambda*dist+ log(t1);
    }
  return logprior;
}


real [] ppnormvec(vector x, real alpha){    
      int N = num_elements(x);
      real  dens[N];
      for( n in 1:N){
     dens[n] = pow(Phi(x[n]),alpha);
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
                       return bernoulli_lupmf(y |ppnormvec(beta0 + beta[1] * x[start:end] + beta[2] * z[start:end] + 
                      beta[3]*w[start:end] + beta[4]*x1[start:end] + beta[5]*x2[start:end],alpha));
                     // + beta[3] * w[start:end]
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
  real<lower=0,upper = 30> alpha1;
}

transformed parameters{
real beta0 = inv_Phi(q^(1/alpha1));
}

// model specification
model { 
 real lambda =  2.903899;
  q ~ beta(1/alpha1,1);   // prior distribution for q (commented)
alpha1~pcalphapnorm(lambda,x_r,x_i); // pcprior distribution for alpha
  beta[1] ~ normal (0,100); // prior distribution for beta0
  beta[2] ~ normal (0,100); // prior distribution for beta1
    beta[3] ~ normal (0,100); // prior distribution for beta1
      beta[4] ~ normal (0,100); // prior distribution for beta1
        beta[5] ~ normal (0,100); // prior distribution for beta1
   target += reduce_sum(partial_sum_lupmf, y, grainsize, alpha1,
                       x, z,w,x1,x2, beta,beta0);
  
}



