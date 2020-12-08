functions{
  real Jpv(real v){
    real y;
    y =trigamma(v/2) -trigamma((v+1)/2) - 2*(v+3)/(v*(v+1)*(v+1));
    y = (v/(v+3))*y;
    return sqrt(y);
  }
}
data {
  int<lower=0> n;             // Time series length
  int<lower=0>d1;             // number of independent variables
  matrix[n,d1] xreg;          // matrix with independent variables
  int<lower=0> period;        // Time series period
  vector[n] y;                // TIme series data
  // conditionals
  int<lower=0,upper=1> is_td; // Has a trend
  int<lower=0,upper=1> is_dp; // Has damped trend
  int<lower=0,upper=1> is_ss; // Has seasonality
  int<lower=0,upper=1> genT;  // Generalized t-student
  // Priors
  vector[4] prior_sigma0;    // prior scale parameter
  vector[4] prior_dfv;       // prior defree freedom genT
  vector[4] prior_level;     // prior defree freedom genT
  vector[4] prior_trend;     // prior defree freedom genT
  vector[4] prior_damped;    // prior defree freedom genT
  vector[4] prior_seasonal;  // prior defree freedom genT
  matrix[d1,4] prior_breg;   // prior ma hyper parameters
  // priors initial values
  vector[4] prior_level1;    // prior inital level
  vector[4] prior_trend1;    // prior inital trend
  vector[4] prior_seasonal1; // prior inital seasonal
}
transformed data{
  int m = period;
  if(is_ss == 0) m = 1;
}
parameters {
  vector[d1] breg;                        // regression parameters
  real<lower=0> sigma0;                   // Scale parameter
  real<lower=0,upper=1> level;            // The level parameter
  vector<lower=0,upper=1>[is_td]trend;    // The Trend parameter
  vector<lower=0,upper=1>[is_dp] damped;  // The damped Trend parameter
  vector<lower=0,upper=1>[is_ss] seasonal;// The seasonal parameter
  vector<lower=2.01>[genT*1] v;           // Degree fredom
  vector<lower=1>[genT*1] lambda;         // lambda parameter

  // Initial values
  real level1;
  vector[is_td] trend1;
  vector[m*is_ss] seasonal1;
}
transformed parameters{
  vector[n] mu;              // the location paramters
  vector[n] epsilon;         // The models error
  vector[n] l;               // The levels   ts
  vector[is_td*n] b;         // The trends   ts
  vector[is_ss*n] s;         // The seasonal ts
  vector<lower=1>[genT*1] v1;// Degree fredom

  // regression estimation
  if(d1 > 0) mu = xreg*breg;
  else mu = rep_vector(0,n);

  // initial values
  l[1] = level1;mu[1] += level1;
  epsilon[1] = y[1] - mu[1];
  if(is_td == 1) b[1] = trend1[1];
  if(is_ss == 1) for(i in 1:m) s[i] = seasonal1[i];

  // Degrees freedom t-student innovations
  if(genT == 1) v1[1] = sqrt((v[1]-2)*lambda[1]/v[1]);


  // ets filter
  for (i in 2:n){

    // local level
    if(is_ss == 1 && i > m) l[i] = level*(y[i] - s[i-m]) + (1 - level)*l[i-1];
    else  l[i] = level*y[i] + (1 - level)*l[i-1];

    mu[i] += l[i-1];

    // Trend component
    if(is_td == 1){
      // damped correction
      if(is_dp == 1) b[i-1] = damped[1]*b[i-1];
      b[i] = trend[1]*(l[i] - l[i-1]) + (1 - trend[1])*b[i-1];
      mu[i] += b[i-1];
    }

    // Seasonal component
    if( is_ss == 1){
       if(i > m) s[i] = seasonal[1]*(y[i] - l[i]) + (1 - seasonal[1])*s[i-m];
      mu[i] += s[i];
    }
    // Errors estimate
    epsilon[i] = y[i] - mu[i];
  }
}
model {
   // priors

  // Prior sigma0
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
  // prior level
  if(prior_level[4] == 1) target += normal_lpdf(level|prior_level[1],prior_level[2]);
  else if(prior_level[4]==2) target += beta_lpdf(level|prior_level[1],prior_level[2]);
  else if(prior_level[4]==3) target += uniform_lpdf(level|prior_level[1],prior_level[2]);

  //  prior for initial level
  if(prior_level1[4] == 1)    target += normal_lpdf(level1|prior_level1[1],prior_level1[2]);
  else if(prior_level1[4]==2) target += beta_lpdf(level1|prior_level1[1],prior_level1[2]);
  else if(prior_level1[4]==3) target += uniform_lpdf(level1|prior_level1[1],prior_level1[2]);
  else if(prior_level1[4]==4) target += student_t_lpdf(level1|prior_level1[3],prior_level1[1],prior_level1[2]);
  else if(prior_level1[4]==5) target += cauchy_lpdf(level1|prior_level1[1],prior_level1[2]);
  else if(prior_level1[4]==6) target += inv_gamma_lpdf(level1|prior_level1[1],prior_level1[2]);
  else if(prior_level1[4]==7) target += inv_chi_square_lpdf(level1|prior_level1[3]);
  else if(prior_level1[4]==8) target += -log(sigma0);
  else if(prior_level1[4]==9) target += gamma_lpdf(level1|prior_level1[1],prior_level1[2]);
  else if(prior_level1[4]==10)target += exponential_lpdf(level1|prior_level1[2]);
  else if(prior_level1[4]==11)target += chi_square_lpdf(level1|prior_level1[3]);
  else if(prior_level1[4]==12)target += double_exponential_lpdf(level1|prior_level1[1],prior_level1[2]);

  // Prior trend
  if(is_td == 1){
    if(prior_trend[4] == 1) target += normal_lpdf(trend[1]|prior_trend[1],prior_trend[2]);
    else if(prior_trend[4]==2) target += beta_lpdf(trend[1]|prior_trend[1],prior_trend[2]);
    else if(prior_trend[4]==3) target += uniform_lpdf(trend[1]|prior_trend[1],prior_trend[2]);

    //  prior for initial trend
    if(prior_trend1[4] == 1)    target += normal_lpdf(trend1|prior_trend1[1],prior_trend1[2]);
    else if(prior_trend1[4]==2) target += beta_lpdf(trend1|prior_trend1[1],prior_trend1[2]);
    else if(prior_trend1[4]==3) target += uniform_lpdf(trend1|prior_trend1[1],prior_trend1[2]);
    else if(prior_trend1[4]==4) target += student_t_lpdf(trend1|prior_trend1[3],prior_trend1[1],prior_trend1[2]);
    else if(prior_trend1[4]==5) target += cauchy_lpdf(trend1|prior_trend1[1],prior_trend1[2]);
    else if(prior_trend1[4]==6) target += inv_gamma_lpdf(trend1|prior_trend1[1],prior_trend1[2]);
    else if(prior_trend1[4]==7) target += inv_chi_square_lpdf(trend1|prior_trend1[3]);
    else if(prior_trend1[4]==8) target += -log(sigma0);
    else if(prior_trend1[4]==9) target += gamma_lpdf(trend1|prior_trend1[1],prior_trend1[2]);
    else if(prior_trend1[4]==10)target += exponential_lpdf(trend1|prior_trend1[2]);
    else if(prior_trend1[4]==11)target += chi_square_lpdf(trend1|prior_trend1[3]);
    else if(prior_trend1[4]==12)target += double_exponential_lpdf(trend1|prior_trend1[1],prior_trend1[2]);
  }

  // Prior Damped Trend
   if(is_dp == 1){
    if(prior_damped[4] == 1) target += normal_lpdf(damped[1]|prior_damped[1],prior_damped[2]);
    else if(prior_damped[4]==2) target += beta_lpdf(damped[1]|prior_damped[1],prior_damped[2]);
    else if(prior_damped[4]==3) target += uniform_lpdf(damped[1]|prior_damped[1],prior_damped[2]);
   }

   // Prior Damped Trend
   if(is_ss == 1){
    if(prior_seasonal[4] == 1) target += normal_lpdf(seasonal[1]|prior_seasonal[1],prior_seasonal[2]);
    else if(prior_seasonal[4]==2) target += beta_lpdf(seasonal[1]|prior_seasonal[1],prior_seasonal[2]);
    else if(prior_seasonal[4]==3) target += uniform_lpdf(seasonal[1]|prior_seasonal[1],prior_seasonal[2]);

    //  prior for initial Seasonal
    if(prior_seasonal1[4] == 1)    target += normal_lpdf(seasonal1|prior_seasonal[1],prior_seasonal[2]);
    else if(prior_seasonal1[4]==2) target += beta_lpdf(seasonal1|prior_seasonal1[1],prior_seasonal1[2]);
    else if(prior_seasonal1[4]==3) target += uniform_lpdf(seasonal1|prior_seasonal1[1],prior_seasonal1[2]);
    else if(prior_seasonal1[4]==4) target += student_t_lpdf(seasonal1|prior_seasonal1[3],prior_seasonal1[1],prior_seasonal1[2]);
    else if(prior_seasonal1[4]==5) target += cauchy_lpdf(seasonal1|prior_seasonal1[1],prior_seasonal1[2]);
    else if(prior_seasonal1[4]==6) target += inv_gamma_lpdf(seasonal1|prior_seasonal1[1],prior_seasonal1[2]);
    else if(prior_seasonal1[4]==7) target += inv_chi_square_lpdf(seasonal1|prior_seasonal1[3]);
    else if(prior_seasonal1[4]==8) target += -log(sigma0);
    else if(prior_seasonal1[4]==9) target += gamma_lpdf(seasonal1|prior_seasonal1[1],prior_seasonal1[2]);
    else if(prior_seasonal1[4]==10)target += exponential_lpdf(seasonal1|prior_seasonal1[2]);
    else if(prior_seasonal1[4]==11)target += chi_square_lpdf(seasonal1|prior_seasonal1[3]);
    else if(prior_seasonal1[4]==12)target += double_exponential_lpdf(seasonal1|prior_seasonal1[1],prior_seasonal1[2]);
   }

  // Likelihood
  if(genT == 1){
        // Prior dfv
    if(prior_dfv[4] == 1) target += normal_lpdf(v[1]|prior_dfv[1],prior_dfv[2]);
    else if(prior_dfv[4]==2) target += beta_lpdf(v[1]|prior_dfv[1],prior_dfv[2]);
    else if(prior_dfv[4]==3) target += uniform_lpdf(v[1]|prior_dfv[1],prior_dfv[2]);
    else if(prior_dfv[4]==4) target += student_t_lpdf(v[1]|prior_dfv[3],prior_dfv[1],prior_dfv[2]);
    else if(prior_dfv[4]==5) target += cauchy_lpdf(v[1]|prior_dfv[1],prior_dfv[2]);
    else if(prior_dfv[4]==6) target += inv_gamma_lpdf(v[1]|prior_dfv[1],prior_dfv[2]);
    else if(prior_dfv[4]==7) target += inv_chi_square_lpdf(v[1]|prior_dfv[3]);
    else if(prior_dfv[4]==8) target += log(Jpv(v[1]));
    else if(prior_dfv[4]==9) target += gamma_lpdf(v[1]|prior_dfv[1],prior_dfv[2]);
    else if(prior_dfv[4]==10)target += exponential_lpdf(v[1]|prior_dfv[2]);
    else if(prior_dfv[4]==11)target += chi_square_lpdf(v[1]|prior_dfv[3]);
    else if(prior_dfv[4]==12)target += double_exponential_lpdf(v[1]|prior_dfv[1],prior_dfv[2]);

    // Likelihood
    target += gamma_lpdf(v[1]|2,0.1);
    target += inv_gamma_lpdf(lambda[1]|v[1]/2,v[1]/2);
    target += normal_lpdf(epsilon|0,v1[1]*sigma0);
  }
  else
    target += normal_lpdf(epsilon|0,sigma0);
}
generated quantities{
  real loglik = 0;
  vector[n] log_lik;
  vector[n] fit;
  vector[n] residuals;

  for(i in 1:n){
    if(genT == 1){
      residuals[i] = student_t_rng(v[1],epsilon[i],sigma0);
     log_lik[i] = student_t_lpdf(y[i]|v[1],mu[i],sigma0);
     loglik += log_lik[i];
    }
    else{
      residuals[i] = normal_rng(epsilon[i],sigma0);
      log_lik[i] = normal_lpdf(y[i]|mu[i],sigma0);
      loglik += log_lik[i];
    }
  }
  fit = y - residuals;
}
