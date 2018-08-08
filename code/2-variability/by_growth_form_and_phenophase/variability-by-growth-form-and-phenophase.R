#setwd("D:/Dropbox/Ecology/Harvard/PhenologyScaling/Analysis/2-variability")
setwd("/media/mkosmala/Data/Dropbox/Ecology/Harvard/PhenologyScaling/Analysis/2-variability/by_growth_form_and_phenophase")

options(max.print = 100000)

library("rstan") # observe startup messages
library("reshape")

# parallel
#rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())

# read in data
data = read.csv("../../all_phenology_dates_18_days_or_less_diff_20180319.csv")

# need to break data up into site-spp-phenophase categories, by year
#data$site.spp.pheno = paste(data$site.species,data$phenophase,sep=" ")
data$site.spp.year.gf.pheno = paste(paste(paste(data$site.species,data$year,sep="-"),data$growthform,sep="-"),data$phenophase,sep="-")

# specifically map site-spp-yr-gf-pheno to integers for reverse lookup
counts = as.data.frame(table(data$site.spp.year.gf.pheno))
counts$index = seq(1,dim(counts)[1])
counts$mean = tapply(data$estimate,data$site.spp.year.gf.pheno,mean)
counts$stdev = tapply(data$estimate,data$site.spp.year.gf.pheno,sd)
colnames(counts) = c("site-spp-yr-gf-pheno","n","index","mean","stdev")
write.csv(counts,"site-spp-yr-gf-pheno_lookup.csv",quote=FALSE,row.names=FALSE)

# map site-spp-pheno to integers
#sitesppheno = as.data.frame(table(data$site.spp.pheno))
#sitesppheno$index = seq(1,dim(sitesppheno)[1])
#sitesppheno$mean = tapply(data$estimate,data$site.spp.pheno,mean)
#sitesppheno$stdev = tapply(data$estimate,data$site.spp.pheno,sd)
#colnames(sitesppheno) = c("site-spp-pheno","n","index","mean","stdev")
#write.csv(sitesppheno,"site-spp-pheno_lookup.csv",quote=FALSE,row.names=FALSE)

# map growthform-phenophases to integers
gfpheno = as.data.frame(table(data$growthform,data$phenophase))
colnames(gfpheno) = c('growthform','phenophase','n')
  
tmp = as.data.frame(with(data, tapply(estimate, list(growthform, phenophase), mean)))
tmp$growthform = rownames(tmp)
means = melt(tmp,id="growthform")
colnames(means) = c('growthform','phenophase','mean')

tmp = as.data.frame(with(data, tapply(estimate, list(growthform, phenophase), sd)))
tmp$growthform = rownames(tmp)
stdev = melt(tmp,id="growthform")
colnames(stdev) = c('growthform','phenophase','stdev')

gfpheno = merge(merge(gfpheno,means),stdev)

gfpheno = gfpheno[gfpheno$n > 0,]  # get rid of rows with no data
gfpheno$index = seq(1,dim(gfpheno)[1])

write.csv(gfpheno,"growthform-phenophase_lookup.csv",quote=FALSE,row.names=FALSE)


# join everything to get data frame for model
d1 = merge(data,gfpheno)
names(d1)[names(d1) == 'index'] <- 'gfpheno_index'

d3 = merge(d1,counts,by.x="site.spp.year.gf.pheno",by.y="site-spp-yr-gf-pheno")
names(d3)[names(d3) == 'index'] <- 'sitespyrgfpheno_index'

# set up and run the model
model_dat = list(N = nrow(data), 
                 K = nrow(counts),
                 #L = nrow(sitesppheno),
                 M = nrow(gfpheno),
                 sitespyrgfpheno = d3$sitespyrgfpheno_index,
                 #sitesppheno = d3$sitesppheno_index,
                 gfpheno = d3$gfpheno_index,
                 y = d3$estimate)

fit = stan(file = "variability-with-growthform.stan", 
           data = model_dat, 
           iter = 1000, chains = 4)

# analyze
#print(fit)
#print(fit,pars=c("mu_j","sigma_j","mu_k","sigma_k"))

#iters = as.data.frame(fit,pars=c("mu_j","sigma_j","mu_k","sigma_k"))

fit_summary = summary(fit,pars=c("mu_j","sigma_j","mu_k","sigma_k"))$summary
#print(fit_summary)

# write results to file
write.csv(as.data.frame(fit_summary),"parameter_values_model_partially_pooled_by_growthform_and_pheno.csv",quote=FALSE)


