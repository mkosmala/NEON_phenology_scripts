// partially pool variances across years instead of treating
// each year at a site as independent

data {
  int<lower=1> N; // number of data points
  int<lower=1> K; // number of site-spp-yr-pheno combinations
  int<lower=1> L; // number of site-spp-pheno combinations
  int<lower=1> M; // number of phenophases
  int sitespyrpheno[N]; // site-species-yr-pheno lookup
  int sitesppheno[N]; // site-species-pheno lookup
  int pheno[N]; // phenophase lookup
  vector[N] y; // date estimate (doy)
}
transformed data {
  int pheno_lookup[K];

  // create lookup table for site-spp-yr-pheno to pheno
  for (i in 1:N) {
    pheno_lookup[sitespyrpheno[i]] = pheno[i];
  }
}
parameters {

  // one entry for each site-spp-yr-pheno 
  vector<lower=1,upper=366>[K] mu_j;
  vector<lower=0>[K] sigma_j;

  // one entry for each site-spp-pheno
  // this is the mean and stdev of the stdev of the _j vars
  vector<lower=0>[M] mu_k;
  vector<lower=0>[M] sigma_k;
}
transformed parameters 
{
  vector[N] y_hat;
  vector[N] sigma_hat;
  vector[K] pheno_mu;
  vector[K] pheno_sigma;

  for (i in 1:N) {
    y_hat[i] = mu_j[sitespyrpheno[i]];  // leave unmodeled
    sigma_hat[i] = exp(sigma_j[sitespyrpheno[i]]); 
  }

  for (i in 1:K) {
    pheno_mu[i] = mu_k[pheno_lookup[i]]; // leave unmodeled
    pheno_sigma[i] = exp(sigma_k[pheno_lookup[i]]); // leave unmodeled
  }
}
model {

  // priors on all the mus and sigmas
  for(i in 1:K) {
    mu_j[i] ~ uniform(1,366);
    // sigma_j[i] ~ cauchy(0,10);
  }
  for(i in 1:M) {
    mu_k[i] ~ cauchy(0,10);
    sigma_k[i] ~ cauchy(0,10);
  }

  // model
  y ~ normal(y_hat, sigma_hat);

  // partial pooling of variance across phenophases
  sigma_j ~ normal(pheno_mu,pheno_sigma);
}


