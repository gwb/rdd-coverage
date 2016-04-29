data {
  int<lower=1> N1;
  int<lower=1> N2;

  vector[N1] x1;
  vector[N1] y1;

  vector[N2] x2;
  vector[N2] y2;
}

transformed data {
  real<lower=0> nug_w;
  vector[N1] mu_w1;
  vector[N2] mu_w2;

  for (i in 1:N1) {
    mu_w1[i] <- 0;
  }
  for (i in 1:N2) {
    mu_w2[i] <- 0;
  }

  nug_w <- 0.001;
}

parameters {
  real<lower=0> phi;
  real<lower=0> sigmasq;
  real<lower=0> etasq;
  real mu1;
  real mu2;
  vector[N1] w_s1;
  vector[N2] w_s2;
}

transformed parameters {
  vector[N1] mu_vec1;
  vector[N2] mu_vec2;

  for(i in 1:N1){
    mu_vec1[i] <- mu1 + w_s1[i]; 
  }
  for(i in 1:N2){
    mu_vec2[i] <- mu2 + w_s2[i]; 
  }
}

model {
  matrix[N1,N1] Kappa1;
  matrix[N2,N2] Kappa2;
  matrix[N1,N1] Sigma_I1;
  matrix[N2,N2] Sigma_I2;      
  
  // off-diagonal elements
  for(i in 1:(N1-1)){
    for(j in (i+1):N1){
      Kappa1[i,j] <- sigmasq*exp( -phi*pow(x1[i] - x1[j],2));
      Kappa1[j,i] <- Kappa1[i,j];
      Sigma_I1[i,j] <- 0;
      Sigma_I1[j,i] <- 0;
    }
  }
  for(i in 1:(N2-1)){
    for(j in (i+1):N2){
      Kappa2[i,j] <- sigmasq*exp( -phi*pow(x2[i] - x2[j],2));
      Kappa2[j,i] <- Kappa2[i,j];
      Sigma_I2[i,j] <- 0;
      Sigma_I2[j,i] <- 0;
    }
  }

  // diagonal elements
  for(k in 1:N1){
    Kappa1[k,k] <- sigmasq + nug_w;
    Sigma_I1[k,k] <- etasq;
  }
  for(k in 1:N2){
    Kappa2[k,k] <- sigmasq + nug_w;
    Sigma_I2[k,k] <- etasq;
  }
  
  mu1 ~ normal(0, 100);
  mu2 ~ normal(0, 100);

  phi ~ cauchy(0,5);
  sigmasq ~ cauchy(0,5);
  etasq ~ cauchy(0,5);
  
  
  w_s1 ~ multi_normal(mu_w1, Kappa1);
  w_s2 ~ multi_normal(mu_w2, Kappa2);
  
  y1 ~ multi_normal(mu_vec1, Sigma_I1);
  y2 ~ multi_normal(mu_vec2, Sigma_I2);
}
