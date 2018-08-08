#setwd("D:/Dropbox/Ecology/Harvard/PhenologyScaling/Analysis/2-variability")
setwd("/media/mkosmala/Data/Dropbox/Ecology/Harvard/PhenologyScaling/Analysis/2-variability/by_phenophase")

options(max.print = 100000)

library("rstan") # observe startup messages

# parallel
#rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())

# read in data
data = read.csv("../all_phenology_dates_18_days_or_less_diff_20180319.csv")

# need to break data up into site-spp-phenophase categories, by year
data$site.spp.pheno = paste(data$site.species,data$phenophase,sep=" ")
data$site.spp.pheno.year = paste(paste(data$site.species,data$phenophase,sep=" "),data$year,sep=" ")

# specifically map site-spp-pheno-yr to integers for reverse lookup
counts = as.data.frame(table(data$site.spp.pheno.year))
counts$index = seq(1,dim(counts)[1])
counts$mean = tapply(data$estimate,data$site.spp.pheno.year,mean)
counts$stdev = tapply(data$estimate,data$site.spp.pheno.year,sd)
colnames(counts) = c("site-spp-pheno-yr","n","index","mean","stdev")
write.csv(counts,"site-spp-pheno-yr_lookup.csv",quote=FALSE,row.names=FALSE)

# map site-spp-pheno to integers
sitesppheno = as.data.frame(table(data$site.spp.pheno))
sitesppheno$index = seq(1,dim(sitesppheno)[1])
sitesppheno$mean = tapply(data$estimate,data$site.spp.pheno,mean)
sitesppheno$stdev = tapply(data$estimate,data$site.spp.pheno,sd)
colnames(sitesppheno) = c("site-spp-pheno","n","index","mean","stdev")
write.csv(sitesppheno,"site-spp-pheno_lookup.csv",quote=FALSE,row.names=FALSE)

# map phenophases to integers
pheno = as.data.frame(table(data$phenophase))
pheno$index = seq(1,dim(pheno)[1])
pheno$mean = tapply(data$estimate,data$phenophase,mean)
pheno$stdev = tapply(data$estimate,data$phenophase,sd)
colnames(pheno) = c("pheno","n","index","mean","stdev")
write.csv(pheno,"pheno_lookup.csv",quote=FALSE,row.names=FALSE)

# join everything to get data frame for model
d1 = merge(data,pheno,by.x="phenophase",by.y="pheno")
names(d1)[names(d1) == 'index'] <- 'pheno_index'
d2 = merge(d1,sitesppheno,by.x="site.spp.pheno",by.y="site-spp-pheno")
names(d2)[names(d2) == 'index'] <- 'sitesppheno_index'
d3 = merge(d2,counts,by.x="site.spp.pheno.year",by.y="site-spp-pheno-yr")
names(d3)[names(d3) == 'index'] <- 'sitespphenoyr_index'

# set up and run the model
model_dat = list(N = nrow(data), 
                 K = nrow(counts),
                 L = nrow(sitesppheno),
                 M = nrow(pheno),
                 sitespyrpheno = d3$sitespphenoyr_index,
                 sitesppheno = d3$sitesppheno_index,
                 pheno = d3$pheno_index,
                 y = d3$estimate)

fit = stan(file = "variability-pooled-by-site-spp-then-by-pheno.stan", 
           data = model_dat, 
           iter = 1000, chains = 4)

# analyze
#print(fit)
print(fit,pars=c("mu_j","sigma_j","mu_k","sigma_k"))

#iters = as.data.frame(fit,pars=c("mu_j","sigma_j","mu_k","sigma_k"))

fit_summary = summary(fit,pars=c("mu_j","sigma_j","mu_k","sigma_k"))$summary
#print(fit_summary)

# write results to file
write.csv(as.data.frame(fit_summary),"parameter_values_model_partially_pooled_by_pheno.csv",quote=FALSE)


