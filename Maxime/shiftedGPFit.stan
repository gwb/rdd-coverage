data {
  int<lower=1> N;
  vector[N] x;
  vector[N] y;
}

transformed data {
  vector[N] mu_w;
  real<lower=0> nug_w;

  for (i in 1:N) {
    mu_w[i] <- y[i]/4;
  }

  nug_w <- 0.001;
}

parameters {
  real shift;
  real beta;
  real<lower=0> sigmasq;
  real<lower=0> phi;
  real<lower=0> tausq;
  
  vector[N] w_s;
}

transformed parameters {
  vector[N] mu_vec;
  
  for(i in 1:N){
    mu_vec[i] <- if_else(x[i] < 5, 0, shift) + x[i] * beta + w_s[i]; 
  }
}

model {
  matrix[N,N] Kappa;
  matrix[N,N] Sigma_I;      
  
  // off-diagonal elements
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      Kappa[i,j] <- sigmasq * exp( -phi * pow(x[i] - x[j],2));
      Kappa[j,i] <- Kappa[i,j];
      
      Sigma_I[i,j] <- 0;
      Sigma_I[j,i] <- 0;
    }
  }
  
  // diagonal elements
  for(k in 1:N){
    Sigma_I[k,k] <- tausq;
    Kappa[k,k] <- sigmasq + nug_w;
  }
  
  sigmasq ~ inv_gamma(2, 0.1);
  phi ~ gamma(2, 0.1);
  tausq ~ inv_gamma(2, 0.1);
  shift ~ normal(0, 10);
  beta ~ normal(0, 10);
  
  
  w_s ~ multi_normal(mu_w, Kappa);
  
  y ~ multi_normal(mu_vec, Sigma_I);
}
