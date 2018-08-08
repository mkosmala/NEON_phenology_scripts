setwd("/media/mkosmala/Data/Dropbox/Ecology/Harvard/PhenologyScaling/Analysis/1-uncertainty")
#setwd("D:/Dropbox/Ecology/Harvard/PhenologyScaling/Analysis/1-uncertainty")
options(max.print = 100000)

library("rstan") # observe startup messages

# parallel
#rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())

# read in data
data = read.csv("../all_phenology_dates_18_days_or_less_diff_20180319.csv")

# set up and run the model
model_dat = list(J = nrow(data), y = unlist(data["doy.diff"]))

fit = stan(file = "uncertainty_mean_only.stan", data = model_dat, 
           iter = 1000, chains = 4)

# analyze
print(fit)

# do by phenophase types
# break data frame up into one data frame for each phenophase
list_by_phenophase = split(data, f = data$phenophase)
num_phenophases = length(list_by_phenophase)

# and run the model for each type
all_fits = list()
for (i in 1:num_phenophases) {

  ndata = list_by_phenophase[[i]]
  model_dat = list(J = nrow(ndata), y = unlist(ndata["doy.diff"]))

  all_fits[[i]] = stan(fit=fit, data = model_dat, 
                     iter = 1000, chains = 4)
}

# analyze
for (i in 1:num_phenophases) {
  print(names(list_by_phenophase[i]))
  print(all_fits[[i]])
}


# add site-species info

# specifically map site-spp to integers for reverse lookup
sitespp = unique(data["site.species"])
sppnum = order(sitespp)
sitespp_lookup = data.frame(seq(1,length(sppnum)),sitespp[sppnum,])
colnames(sitespp_lookup) = c("index","site-spp")
write.csv(sitespp_lookup,"site-spp_lookup.csv",quote=FALSE,row.names=FALSE)

# run the model
site_spp_nums = as.numeric(unlist(data["site.species"]))
model_dat = list(N = nrow(data), 
                 K = max(site_spp_nums),
                 sitesp = site_spp_nums,
                 y = unlist(data["doy.diff"]))

fit = stan(file = "uncertainty.stan", data = model_dat, 
           iter = 1000, chains = 4)


# analyze
print(fit,pars = c("mu","tau","sigma","gamma"))
plot(fit, pars = c("mu","tau","sigma","gamma"))


###### ONLY THIS ANALYSIS MATTERS

# add site-species info by phenophase types

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


