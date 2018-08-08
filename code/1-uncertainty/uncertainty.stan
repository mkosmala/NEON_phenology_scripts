data {
  int<lower=1> N; // number of data points
  int<lower=1> K; // number of site-spp
  int sitesp[N]; // site-species
  vector[N] y; // uncertainty (days)
}
parameters {
  real<lower=0> mu;
  vector[K] gamma;
  real<lower=0> tau;
  real<lower=0> sigma;
}
transformed parameters {
  vector[N] yhat;

  for (i in 1:N) {
    yhat[i] = mu + gamma[sitesp[i]];
  }
}
model {
  y ~ normal(yhat, sigma);
  gamma ~ normal(0, tau);
}

