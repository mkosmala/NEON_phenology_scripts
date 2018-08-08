// we will run this separately for each phenophase, since we have no
// a priori reason to believe that different phenophases are affected
// by the same factors in the same direction

data {
  int<lower=1> N; // number of data points
  int<lower=1> K; // number of site-spp-yr combinations
//  int<lower=1> L; // number of site-spp combinations
  int sitespyr[N]; // site-species-yr lookup
//  int sitesp[N]; // site-species lookup
  vector[N] y; // date estimate (doy), centered
  vector[N] elev; // elevation (meters), centered
  vector[N] slope; // slope (degrees from horizontal, non-neg)
  vector[N] aspectNS; // north-slope aspect (-1 to 1)
  vector[N] aspectEW; // east-west aspect (-1 to 1)
  vector[N] pheight; // plant height, normalized (0 to 1)
  vector[N] disease; // whether diseaseed (0 or 1)
  vector[N] canopy_full_sun;         // canopy position, dummy-coded (0 or 1)
  vector[N] canopy_mostly_shaded;    // canopy position, dummy-coded (0 or 1)
  vector[N] canopy_open_grown;       // canopy position, dummy-coded (0 or 1)
  vector[N] canopy_partially_shaded; // canopy position, dummy-coded (0 or 1)
}
parameters {
  vector[K] alpha;  // so I don't have to center everything

  vector[K] beta_elev_raw; // slopes for effect of elevation
  vector[K] beta_slope_raw; 
  vector[K] beta_aspectNS_raw;
  vector[K] beta_aspectEW_raw;
  vector[K] beta_slope_aspectNS_raw;
  vector[K] beta_slope_aspectEW_raw;
  //vector[K] beta_apsectNSEW_raw;
  vector[K] beta_pheight_raw;
  vector[K] beta_disease_raw;
  vector[K] beta_canopy_full_sun_raw;
  vector[K] beta_canopy_mostly_shaded_raw;
  vector[K] beta_canopy_open_grown_raw;
  vector[K] beta_canopy_partially_shaded_raw;

  real<lower=0> sigma;
  real gamma_elev;
  real gamma_slope;
  real gamma_aspectNS;
  real gamma_aspectEW;
  real gamma_slope_aspectNS;
  real gamma_slope_aspectEW;
//  real gamma_aspectNSEW;
  real gamma_pheight;
  real gamma_disease;
  real gamma_canopy_full_sun;
  real gamma_canopy_mostly_shaded;
  real gamma_canopy_open_grown;
  real gamma_canopy_partially_shaded;

  real<lower=0> tau_elev;
  real<lower=0> tau_slope;  
  real<lower=0> tau_aspectNS;  
  real<lower=0> tau_aspectEW;  
  real<lower=0> tau_slope_aspectNS;  
  real<lower=0> tau_slope_aspectEW;  
//  real<lower=0> tau_aspectNSEW;  
  real<lower=0> tau_pheight; 
  real<lower=0> tau_disease;
  real<lower=0> tau_canopy_full_sun;
  real<lower=0> tau_canopy_mostly_shaded;
  real<lower=0> tau_canopy_open_grown;
  real<lower=0> tau_canopy_partially_shaded;
 
}
transformed parameters {
  vector[N] yhat;
  vector[K] beta_elev;
  vector[K] beta_slope; 
  vector[K] beta_aspectNS;
  vector[K] beta_aspectEW;
  vector[K] beta_slope_aspectNS;
  vector[K] beta_slope_aspectEW;
  //vector[K] beta_apsectNSEW;
  vector[K] beta_pheight;
  vector[K] beta_disease;
  vector[K] beta_canopy_full_sun;
  vector[K] beta_canopy_mostly_shaded;
  vector[K] beta_canopy_open_grown;
  vector[K] beta_canopy_partially_shaded;

  // leave alpha unmodeled
  // betas are normally distributed
  beta_elev = gamma_elev + tau_elev * beta_elev_raw;
  beta_slope = gamma_slope + tau_slope * beta_slope_raw;
  beta_aspectNS = gamma_aspectNS + tau_aspectNS * beta_aspectNS_raw;
  beta_aspectEW = gamma_aspectEW + tau_aspectEW * beta_aspectEW_raw;
  beta_slope_aspectNS = gamma_slope_aspectNS + tau_slope_aspectNS * beta_aspectNS_raw;
  beta_slope_aspectEW = gamma_slope_aspectEW + tau_slope_aspectEW * beta_aspectEW_raw;
  //beta_aspectNSEW = gamma_aspectNSEW + tau_aspectNSEW * beta_aspectNSEW_raw;
  beta_pheight = gamma_pheight + tau_pheight * beta_pheight_raw;
  beta_disease = gamma_disease + tau_disease * beta_disease_raw;
  beta_canopy_full_sun = gamma_canopy_full_sun + tau_canopy_full_sun * beta_canopy_full_sun_raw;
  beta_canopy_mostly_shaded = gamma_canopy_mostly_shaded + tau_canopy_mostly_shaded * beta_canopy_mostly_shaded_raw;
  beta_canopy_open_grown = gamma_canopy_open_grown + tau_canopy_open_grown * beta_canopy_open_grown_raw;
  beta_canopy_partially_shaded = gamma_canopy_partially_shaded + tau_canopy_partially_shaded * beta_canopy_partially_shaded_raw;



  for (i in 1:N) {
    yhat[i] = alpha[sitespyr[i]] + 
              beta_elev[sitespyr[i]] * elev[i] +
              beta_slope[sitespyr[i]] * slope[i] +
              beta_aspectNS[sitespyr[i]] * aspectNS[i] +
              beta_aspectEW[sitespyr[i]] * aspectEW[i] +
              beta_slope_aspectNS[sitespyr[i]] * aspectNS[i] * slope[i] +
              beta_slope_aspectEW[sitespyr[i]] * aspectEW[i] * slope[i] +
              beta_pheight[sitespyr[i]] * pheight[i] +
              beta_disease[sitespyr[i]] * disease[i] +
              beta_canopy_full_sun[sitespyr[i]] * canopy_full_sun[i] +
              beta_canopy_mostly_shaded[sitespyr[i]] * canopy_mostly_shaded[i] +
              beta_canopy_open_grown[sitespyr[i]] * canopy_open_grown[i] +
              beta_canopy_partially_shaded[sitespyr[i]] * canopy_partially_shaded[i];  
  }
}
model {
  // priors
  //for(i in 1:K) {
  //   alpha[i] ~ uniform(1,366);
  //}
  sigma ~ cauchy(0,10);
  gamma_elev ~ normal(0,10);
  gamma_slope ~ normal(0,10);
  gamma_aspectNS ~ normal(0,10);
  gamma_aspectEW ~ normal(0,10);
  gamma_slope_aspectNS ~ normal(0,10);
  gamma_slope_aspectEW ~ normal(0,10);
//  gamma_aspectNSEW ~ normal(0,10);
  gamma_pheight ~ normal(0,10);
  gamma_disease ~ normal(0,10);
  gamma_canopy_full_sun ~ normal(0,10);
  gamma_canopy_mostly_shaded ~ normal(0,10);
  gamma_canopy_open_grown ~ normal(0,10);
  gamma_canopy_partially_shaded ~ normal(0,10);

  tau_elev ~ cauchy(0,10);  
  tau_slope ~ cauchy(0,10);  
  tau_aspectNS ~ cauchy(0,10);  
  tau_aspectEW ~ cauchy(0,10);  
  tau_slope_aspectNS ~ cauchy(0,10);  
  tau_slope_aspectEW ~ cauchy(0,10); 
//  tau_aspectNSEW ~ cauchy(0,10);  
  tau_pheight ~ cauchy(0,10);  
  tau_disease ~ cauchy(0,10);  
  tau_canopy_full_sun ~ cauchy(0,10);  
  tau_canopy_mostly_shaded ~ cauchy(0,10);  
  tau_canopy_open_grown ~ cauchy(0,10);  
  tau_canopy_partially_shaded ~ cauchy(0,10);  

  // model
  y ~ normal(yhat, sigma);

  // partial pooling of the effect of elevation
  beta_elev_raw ~ normal(0,1);
  beta_slope_raw ~ normal(0,1);
  beta_aspectNS_raw ~ normal(0,1);
  beta_aspectEW_raw ~ normal(0,1);
  beta_slope_aspectNS_raw ~ normal(0,1);
  beta_slope_aspectEW_raw ~ normal(0,1);
//  beta_aspectNSEW_raw ~ normal(0,1);
  beta_pheight_raw ~ normal(0,1);
  beta_disease_raw ~ normal(0,1);
  beta_canopy_full_sun_raw ~ normal(0,1);
  beta_canopy_mostly_shaded_raw ~ normal(0,1);
  beta_canopy_open_grown_raw ~ normal(0,1);
  beta_canopy_partially_shaded_raw ~ normal(0,1);
}



