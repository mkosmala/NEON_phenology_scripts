# This file runs an "uncertainty" analysis that examines to what degree
# of precision each phenophase transition date is known. The result
# is an estimate of uncertainty for each site-species for each phenophase
# transition date (gammas), as well as an overall uncertainty estimate 
# for each phenophase transition date (mus).

setwd("~/Harvard/NEON_phenology/code/1-uncertainty")
options(max.print = 100000)
library("rstan") 

# read in data
data = read.csv("../../data/all_phenology_dates_18_days_or_less_diff_20180319.csv")

# add site-species info by phenophase types

# compile the stan model
# It will return an error, but it successfully compiles
fit = stan(file = "uncertainty.stan", iter=0)

# break data frame up into one data frame for each phenophase
list_by_phenophase = split(data, f = data$phenophase)
num_phenophases = length(list_by_phenophase)

# and run the model for each type
all_fits = list()
match_site_spp = list()
all_lookups = list()
for (i in 1:num_phenophases) {

  ndata = list_by_phenophase[[i]]
  
  # speces lookup
  nsitespp = unique(ndata["site.species"])
  nsppnum = order(nsitespp)
  all_lookups[[i]] = data.frame(seq(1,length(nsppnum)),nsitespp[nsppnum,])
  colnames(all_lookups[[i]]) = c("index","site.species")  
  
  ndatamerge = merge(ndata,all_lookups[[i]])
  nsite_spp_nums = ndatamerge$index
  
  # run the model
  model_dat = list(N = nrow(ndata), 
                   K = max(nsite_spp_nums),
                   sitesp = nsite_spp_nums,
                   y = unlist(ndata["doy.diff"]))

  all_fits[[i]] = stan(fit=fit, data = model_dat, 
                     iter = 1000, chains = 4)

}

# analyze
# output to file
sink('uncertainty-analysis-output.txt')
for (i in 1:num_phenophases) {
  print(names(list_by_phenophase[i]))
  print(all_fits[[i]],pars = c("mu","tau","sigma","gamma"))
  print(all_lookups[[i]],row.names = FALSE)
  print("")
}
sink()



