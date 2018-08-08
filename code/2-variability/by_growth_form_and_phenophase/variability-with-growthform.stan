// Partially pool all site-spp-year units
// letting means float freely and variances be partially pooled
// by growth form and phenophase
// No assumptions about relationships among growth forms
// or among phenophases

data {
  int<lower=1> N; // number of data points
  int<lower=1> K; // number of site-spp-yr-gf-pheno combinations
  int<lower=1> M; // number of growth-form-phenophase combinations
  int sitespyrgfpheno[N]; // site-species-yr-gf-pheno lookup
  int gfpheno[N]; // growth-form-phenophase lookup
  vector[N] y; // date estimate (doy)
}
transformed data {
  int gfpheno_lookup[K];

  // create lookup table for site-spp-yr-gf-pheno to gf-pheno
  for (i in 1:N) {
    gfpheno_lookup[sitespyrgfpheno[i]] = gfpheno[i];
  }
}
parameters {

  // one entry for each site-spp-yr-gf-pheno 
  vector<lower=1,upper=366>[K] mu_j;
  vector<lower=0>[K] sigma_j;

  // one entry for each gf-pheno
  // this is the mean and stdev of the stdev of the _j vars
  vector<lower=0>[M] mu_k;
  vector<lower=0>[M] sigma_k;
}
transformed parameters 
{
  vector[N] y_hat;
  vector[N] sigma_hat;
  vector[K] gfpheno_mu;
  vector[K] gfpheno_sigma;

  for (i in 1:N) {
    y_hat[i] = mu_j[sitespyrgfpheno[i]];  // leave unmodeled
    sigma_hat[i] = exp(sigma_j[sitespyrgfpheno[i]]); 
  }

  for (i in 1:K) {
    gfpheno_mu[i] = mu_k[gfpheno_lookup[i]]; // leave unmodeled
    gfpheno_sigma[i] = exp(sigma_k[gfpheno_lookup[i]]); // leave unmodeled
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
  sigma_j ~ normal(gfpheno_mu,gfpheno_sigma);
}


