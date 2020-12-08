functions{
  real Jpv(real v){
    real y;
    y =trigamma(v/2) -trigamma((v+1)/2) - 2*(v+3)/(v*(v+1)*(v+1));
    y = (v/(v+3))*y;
    return sqrt(y);
  }
    real asymf(real u, real gamma,int asym){
    real y;
    if(asym == 0) y = u;
    if(asym == 1) y = inv_logit(-gamma*u);
    if(asym == 2) y = 1-exp(-gamma*pow(u,2));
    return(pow(u,2)*y);
  }
}
data {
  // Model data
  int<lower=0> n;            // number of data items
  int<lower=0> s;            // number of predictors  arch
  int<lower=0> k;            // number of predictions garch
  int<lower=0> h;            // number of predictions mgarch
  int<lower=0> p;            // number of predictors  ar
  int<lower=0> q;            // number of predictions ma
  vector[n] y;               // outcome time series
  int<lower=0,upper=1> genT; // Generalized t-student
  int<lower=0>d1;            // number of independent variables
  matrix[n,d1] xreg;         // matrix with independent variables
  // prior data
  vector[4] prior_mu0;       // prior location parameter
  vector[4] prior_sigma0;    // prior scale parameter
  vector[4]   prior_dfv;     // prior defree freedom genT
  matrix[p,4] prior_ar;      // ar location hyper parameters
  matrix[q,4] prior_ma;      // ma location hyper parameters
  matrix[s,4] prior_arch;    // prior arch hyper parameters
  matrix[k,4] prior_garch;   // prior ma hyper parameters
  matrix[h,4] prior_mgarch;  // prior ma hyper parameters
  matrix[d1,4] prior_breg;   // prior ma hyper parameters

  // asymmetric garch
  int<lower=0,upper=2> asym;    // Asymmetric garch
  int<lower=0,upper=1> asym1;   // Boolean Asymmetric garch
  matrix[asym1*2,4] prior_gamma;// prior gamma asymmetric parameters
}

parameters{
  real mu0;
  real<lower=0> sigma0;               // Variance parameter
  vector[d1] breg;                    // regression parameters
  vector<lower=-1,upper=1>[p] ar0;   // ar parameters
  vector<lower=-1,upper=1>[q] ma0; // ma parameters
  vector<lower=0,upper=1>[s] arch;   // arch parameters
  vector<lower=0,upper=1>[k] garch;    // garch parameters
  vector[h] mgarch;                   // mean garch parameters
  vector<lower=2.01>[genT*1] v;       // Degree fredom
  vector<lower=1>[genT*n] lambda;     // lambda parameter
  vector<lower=0>[asym1*2] gamma;     // parameter for asymmetric garch;

}
transformed parameters{
  vector[p] ar;          // ar parameters
  vector[q] ma;        // ma parameters

  // Temporal mean and residuals
  vector[n] mu;                     // Mean Parameter
  vector[n] epsilon;                // error parameter
  vector<lower=0>[n] sigma;         // Variance parameter


  //***********************************************
  //         Transformation coeficients
  //***********************************************


  for( i in 1:p){
    if(prior_ar[i,4]== 1) ar[i] = ar0[i];
    else ar[i] = 2*ar0[i] - 1;
  }
  for(i in 1:q){
    if(prior_ma[i,4] == 1) ma[i] = ma0[i];
    else ma[i] = 2*ma0[i]-1;
  }
  // regression estimation
  if(d1 > 0) mu = xreg*breg;
  else mu = rep_vector(0,n);

  //***********************************************
  //         ARMA estimations
  //***********************************************

  for(i in 1:n){
     mu[i] += mu0;
     sigma[i] = sigma0;
    //  ar Estimation
    if(p > 0) for(j in 1:p) if(i > j) mu[i] += y[i-j]*ar[j];
    // ma estimation
    if(q > 0) for(j in 1:q) if(i > j) mu[i] += epsilon[i-j]*ma[j];
    epsilon[i] = y[i] - mu[i];
    // Garch Iteration
    if(s >= k){
       // arch estimation
      if(s > 0) for(j in 1:s) if(i > j) sigma[i] += arch[j]*pow(epsilon[i-j],2);
       // garch estimation
      if(k > 0) for(j in 1:k) if(i > j)  sigma[i] += garch[j]*pow(sigma[i-j], 2);
    }
    // Degrees freedom t-student innovations
    if(genT == 1 ) sigma[i] = sqrt((v[1]-2)*lambda[i]*sigma[i]/v[1]);
    else sigma[i] = sqrt(sigma[i]);
    // mgarch estimation
    if(h > 0) for(j in 1:h)if(i > j) mu[i] += mgarch[j]*sigma[i-j+1];

    // asymetric garch
    if(asym1 == 1) if(i > 1)sigma[i] += gamma[1]*asymf(epsilon[i-1],gamma[2],asym);
  }
}
model {
  //  prior for \mu0
  if(prior_mu0[4] == 1)    target += normal_lpdf(mu0|prior_mu0[1],prior_mu0[2]);
  else if(prior_mu0[4]==2) target += beta_lpdf(mu0|prior_mu0[1],prior_mu0[2]);
  else if(prior_mu0[4]==3) target += uniform_lpdf(mu0|prior_mu0[1],prior_mu0[2]);
  else if(prior_mu0[4]==4) target += student_t_lpdf(mu0|prior_mu0[3],prior_mu0[1],prior_mu0[2]);
  else if(prior_mu0[4]==5) target += cauchy_lpdf(mu0|prior_mu0[1],prior_mu0[2]);
  else if(prior_mu0[4]==6) target += inv_gamma_lpdf(mu0|prior_mu0[1],prior_mu0[2]);
  else if(prior_mu0[4]==7) target += inv_chi_square_lpdf(mu0|prior_mu0[3]);
  else if(prior_mu0[4]==8) target += -log(sigma0);
  else if(prior_mu0[4]==9) target += gamma_lpdf(mu0|prior_mu0[1],prior_mu0[2]);
  else if(prior_mu0[4]==10)target += exponential_lpdf(mu0|prior_mu0[2]);
  else if(prior_mu0[4]==11)target += chi_square_lpdf(mu0|prior_mu0[3]);
  else if(prior_mu0[4]==12)target += double_exponential_lpdf(mu0|prior_mu0[1],prior_mu0[2]);

  // Prior sigma
  if(prior_sigma0[4] == 1)    target += normal_lpdf(sigma0|prior_sigma0[1],prior_sigma0[2]);
  else if(prior_sigma0[4]==2) target += beta_lpdf(sigma0|prior_sigma0[1],prior_sigma0[2]);
  else if(prior_sigma0[4]==3) target += uniform_lpdf(sigma0|prior_sigma0[1],prior_sigma0[2]);
  else if(prior_sigma0[4]==4) target += student_t_lpdf(sigma0|prior_sigma0[3],prior_sigma0[1],prior_sigma0[2]);
  else if(prior_sigma0[4]==5) target += cauchy_lpdf(sigma0|prior_sigma0[1],prior_sigma0[2]);
  else if(prior_sigma0[4]==6) target += inv_gamma_lpdf(sigma0|prior_sigma0[1],prior_sigma0[2]);
  else if(prior_sigma0[4]==7) target += inv_chi_square_lpdf(sigma0|prior_sigma0[3]);
  else if(prior_sigma0[4]==8) target += -log(sigma0);
  else if(prior_sigma0[4]==9) target += gamma_lpdf(sigma0|prior_sigma0[1],prior_sigma0[2]);
  else if(prior_sigma0[4]==10)target += exponential_lpdf(sigma0|prior_sigma0[2]);
  else if(prior_sigma0[4]==11)target += chi_square_lpdf(sigma0|prior_sigma0[3]);
  else if(prior_sigma0[4]==12)target += double_exponential_lpdf(sigma0|prior_sigma0[1],prior_sigma0[2]);

  // prior breg
  if(d1 > 0){
    for(i in 1:d1){
      if(prior_breg[i,4] == 1)    target += normal_lpdf(breg[i]|prior_breg[i,1],prior_breg[i,2]);
      else if(prior_breg[i,4]==2) target += beta_lpdf(breg[i]|prior_breg[i,1],prior_breg[i,2]);
      else if(prior_breg[i,4]==3) target += uniform_lpdf(breg[i]|prior_breg[i,1],prior_breg[i,2]);
      else if(prior_breg[i,4]==4) target += student_t_lpdf(breg[i]|prior_breg[i,3],prior_breg[i,1],prior_breg[i,2]);
      else if(prior_breg[i,4]==5) target += cauchy_lpdf(breg[i]|prior_breg[i,1],prior_breg[i,2]);
      else if(prior_breg[i,4]==6) target += inv_gamma_lpdf(breg[i]|prior_breg[i,1],prior_breg[i,2]);
      else if(prior_breg[i,4]==7) target += inv_chi_square_lpdf(breg[i]|prior_breg[i,3]);
      else if(prior_breg[i,4]==8) target += -log(sigma0);
      else if(prior_breg[i,4]==9) target += gamma_lpdf(breg[i]|prior_breg[i,1],prior_breg[i,2]);
      else if(prior_breg[i,4]==10)target += exponential_lpdf(breg[i]|prior_breg[i,2]);
      else if(prior_breg[i,4]==11)target += chi_square_lpdf(breg[i]|prior_breg[i,3]);
      else if(prior_breg[i,4]==12)target += double_exponential_lpdf(breg[i]|prior_breg[i,1],prior_breg[i,2]);
    }
  }

  // prior ar
  if(p > 0){
    for(i in 1:p){
     if(prior_ar[i,4]==1) target += normal_lpdf(ar0[i]|prior_ar[i,1],prior_ar[i,2]);
     else if(prior_ar[i,4]==2) target += beta_lpdf(fabs(ar0[i])|prior_ar[i,1],prior_ar[i,2]);
     else if(prior_ar[i,4]==3) target += uniform_lpdf(ar0[i]|prior_ar[i,1],prior_ar[i,2]);
    }
  }
  // prior ma
  if(q > 0){
    for(i in 1:q){
      if(prior_ma[i,4]==1) target += normal_lpdf(ma0[i]|prior_ma[i,1],prior_ma[i,2]);
      else if(prior_ma[i,4]==2) target += beta_lpdf(fabs(ma0[i])|prior_ma[i,1],prior_ma[i,2]);
      else if(prior_ma[i,4]==3) target += uniform_lpdf(ma0[i]|prior_ma[i,1],prior_ma[i,2]);
    }
  }

  // prior arch
  if(s > 0){
    for(i in 1:s){
     if(prior_arch[i,4]==1) target += normal_lpdf(arch[i]|prior_arch[i,1],prior_arch[i,2]);
     else if(prior_arch[i,4]==2) target += beta_lpdf(arch[i]|prior_arch[i,1],prior_arch[i,2]);
     else if(prior_arch[i,4]==3) target += uniform_lpdf(arch[i]|prior_arch[i,1],prior_arch[i,2]);

    }
  }
  // prior garch
  if(k > 0){
    for(i in 1:k){
     if(prior_garch[i,4]==1) target += normal_lpdf(garch[i]|prior_garch[i,1],prior_garch[i,2]);
     else if(prior_garch[i,4]==2) target += beta_lpdf(garch[i]|prior_garch[i,1],prior_garch[i,2]);
     else if(prior_garch[i,4]==3) target += uniform_lpdf(garch[i]|prior_garch[i,1],prior_garch[i,2]);
    }
  }
  // prior mean_garch
  if(h > 0){
    for(i in 1:h){
      if(prior_mgarch[i,4]== 1)     target += normal_lpdf(mgarch[i]|prior_mgarch[i,1],prior_mgarch[i,2]);
      else if(prior_mgarch[i,4]==2) target += beta_lpdf(mgarch[i]|prior_mgarch[i,1],prior_mgarch[i,2]);
      else if(prior_mgarch[i,4]==3) target += uniform_lpdf(mgarch[i]|prior_mgarch[i,1],prior_mgarch[i,2]);
      else if(prior_mgarch[i,4]==4) target += student_t_lpdf(mgarch[i]|prior_mgarch[i,3],prior_mgarch[i,1],prior_mgarch[i,2]);
      else if(prior_mgarch[i,4]==5) target += cauchy_lpdf(mgarch[i]|prior_mgarch[i,1],prior_mgarch[i,2]);
      else if(prior_mgarch[i,4]==6) target += inv_gamma_lpdf(mgarch[i]|prior_mgarch[i,1],prior_mgarch[i,2]);
      else if(prior_mgarch[i,4]==7) target += inv_chi_square_lpdf(mgarch[i]|prior_mgarch[i,3]);
      else if(prior_mgarch[i,4]==8) target += -log(sigma0);
      else if(prior_mgarch[i,4]==9) target += gamma_lpdf(mgarch[i]|prior_mgarch[i,1],prior_mgarch[i,2]);
      else if(prior_mgarch[i,4]==10)target += exponential_lpdf(mgarch[i]|prior_mgarch[i,2]);
      else if(prior_mgarch[i,4]==11)target += chi_square_lpdf(mgarch[i]|prior_mgarch[i,3]);
      else if(prior_mgarch[i,4]==12)target += double_exponential_lpdf(mgarch[i]|prior_mgarch[i,1],prior_mgarch[i,2]);
    }
  }
  if(genT == 1){
    // Prior dfv
    if(prior_dfv[4] == 1) target += normal_lpdf(v[1]|prior_dfv[1],prior_dfv[2]);
    else if(prior_dfv[4]==2) target += beta_lpdf(v[1]|prior_dfv[1],prior_dfv[2]);
    else if(prior_dfv[4]==3) target += uniform_lpdf(v[1]|prior_dfv[1],prior_dfv[2]);
    else if(prior_dfv[4]==4) target += student_t_lpdf(v[1]|prior_dfv[3],prior_dfv[1],prior_dfv[2]);
    else if(prior_dfv[4]==5) target += cauchy_lpdf(v[1]|prior_dfv[1],prior_dfv[2]);
    else if(prior_dfv[4]==6) target += inv_gamma_lpdf(v[1]|prior_dfv[1],prior_dfv[2]);
    else if(prior_dfv[4]==7) target += inv_chi_square_lpdf(v[1]|prior_dfv[3]);
    else if(prior_dfv[4] == 8) target += log(Jpv(v[1]));
    else if(prior_dfv[4]==9) target += gamma_lpdf(v[1]|prior_dfv[1],prior_dfv[2]);
    else if(prior_dfv[4]==10)target += exponential_lpdf(v[1]|prior_dfv[2]);
    else if(prior_dfv[4]==11)target += chi_square_lpdf(v[1]|prior_dfv[3]);
    else if(prior_dfv[4]==12)target += double_exponential_lpdf(v[1]|prior_dfv[1],prior_dfv[2]);
  }
  // asymmetric garch
  if(asym1 == 1){
    for(i in 1:2){
      if(prior_gamma[i,4]== 1)     target += normal_lpdf(gamma[i]|prior_gamma[i,1],prior_gamma[i,2]);
      else if(prior_gamma[i,4]==2) target += beta_lpdf(gamma[i]|prior_gamma[i,1],prior_gamma[i,2]);
      else if(prior_gamma[i,4]==3) target += uniform_lpdf(gamma[i]|prior_gamma[i,1],prior_gamma[i,2]);
      else if(prior_gamma[i,4]==4) target += student_t_lpdf(gamma[i]|prior_gamma[i,3],prior_gamma[i,1],prior_gamma[i,2]);
      else if(prior_gamma[i,4]==5) target += cauchy_lpdf(gamma[i]|prior_gamma[i,1],prior_gamma[i,2]);
      else if(prior_gamma[i,4]==6) target += inv_gamma_lpdf(gamma[i]|prior_gamma[i,1],prior_gamma[i,2]);
      else if(prior_gamma[i,4]==7) target += inv_chi_square_lpdf(gamma[i]|prior_gamma[i,3]);
      else if(prior_gamma[i,4]==8) target += -log(sigma0);
      else if(prior_gamma[i,4]==9) target += gamma_lpdf(gamma[i]|prior_gamma[i,1],prior_gamma[i,2]);
      else if(prior_gamma[i,4]==10)target += exponential_lpdf(gamma[i]|prior_gamma[i,2]);
      else if(prior_gamma[i,4]==11)target += chi_square_lpdf(gamma[i]|prior_gamma[i,3]);
      else if(prior_gamma[i,4]==12)target += double_exponential_lpdf(gamma[i]|prior_gamma[i,1],prior_gamma[i,2]);
    }
  }

  // Likelihood
  if(genT == 1)  target+= inv_gamma_lpdf(lambda|v[1]/2,v[1]/2);
  target += normal_lpdf(epsilon|0,sigma);
}
generated quantities{
  real loglik = 0;
  vector[n] log_lik;
  vector[n] fit;
  vector[n] residuals;

  for(i in 1:n){
    if(genT == 1){
      residuals[i] = student_t_rng(v[1],epsilon[i],sigma[i]);
     log_lik[i] = student_t_lpdf(y[i]|v[1],mu[i],sigma[i]);
     loglik += log_lik[i];
    }
    else{
      residuals[i] = normal_rng(epsilon[i],sigma[i]);
      log_lik[i] = normal_lpdf(y[i]|mu[i],sigma[i]);
      loglik += log_lik[i];
    }
  }
  fit = y - residuals;
}
