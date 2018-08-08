#setwd("D:/Dropbox/Ecology/Harvard/PhenologyScaling/Analysis/4-factors")
setwd("/media/mkosmala/Data/Dropbox/Ecology/Harvard/PhenologyScaling/Analysis/4-factors")
options(max.print = 100000)

library("rstan") # observe startup messages

# parallel
#rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())

# read in data
pheno = "first_color_NOT_deciduous_broadleaf"
filename = paste(paste("../plants_and_factors_",pheno,sep=""),".csv",sep="")
data = read.csv(filename)
nrows = dim(data)[1]

# remove rows with NA for plant height or canopy position
data = data[! is.na(data$plant_height),]
data = data[data$canopy_position != "",]

# fix capitalization issues and change expression for dummying
levels(data$canopy_position) = tolower(levels(data$canopy_position))
data$canopy_position <- as.character(data$canopy_position)
data[data[,"canopy_position"]=="full shade",]$canopy_position = "_full_shade"
data[data[,"canopy_position"]=="mostly shaded",]$canopy_position = "_mostly_shaded"
data[data[,"canopy_position"]=="partially shaded",]$canopy_position = "_partially_shaded"
data[data[,"canopy_position"]=="full sun",]$canopy_position = "_full_sun"
data[data[,"canopy_position"]=="open grown",]$canopy_position = "_open_grown"
data$canopy_position <- as.factor(data$canopy_position)



# center doy by site-species-year
# then we can get rid of the alpha paramater, which should make the model run better
mean_doy = as.data.frame(table(data$site_species_year))
mean_doy$mean_doy = tapply(data$estimate,data$site_species_year,mean)
mean_doy$stdev_doy = tapply(data$estimate,data$site_species_year,sd)
mean_doy$min_out_doy = mean_doy$mean_doy - 3*mean_doy$stdev_doy
mean_doy$max_out_doy = mean_doy$mean_doy + 3*mean_doy$stdev_doy
colnames(mean_doy) = c("site-spp-yr","n","mean_doy","stdev_doy","min_out_doy","max_out_doy")
data = merge(data,mean_doy,by.x="site_species_year",by.y="site-spp-yr")
data$centered_doy = data$estimate - data$mean_doy

# center elevation data by site-species-year
mean_elev = as.data.frame(table(data$site_species_year))
mean_elev$mean_elev = tapply(data$elevation,data$site_species_year,mean)
colnames(mean_elev) = c("site-spp-yr","n1","mean_elev")
data = merge(data,mean_elev,by.x="site_species_year",by.y="site-spp-yr")
data$centered_elev = data$elevation - data$mean_elev

# normalize plant height data by site-species (range 0-1)
mean_pheight = as.data.frame(table(data$site_species_year))
mean_pheight$min_pheight = tapply(data$plant_height,data$site_species_year,min)
mean_pheight$max_pheight = tapply(data$plant_height,data$site_species_year,max)
mean_pheight$ave_pheight = tapply(data$plant_height,data$site_species_year,mean)
mean_pheight$sd_pheight = tapply(data$plant_height,data$site_species_year,sd)
mean_pheight$min_out_pheight = mean_pheight$ave_pheight - 3*mean_pheight$sd_pheight
mean_pheight$max_out_pheight = mean_pheight$ave_pheight + 3*mean_pheight$sd_pheight
mean_pheight$height_sum = tapply(data$plant_height,data$site_species_year,sum)
mean_pheight$normalized_mean = (mean_pheight$height_sum - (mean_pheight$Freq * mean_pheight$min_pheight)) /
                               (mean_pheight$Freq * 
                               (mean_pheight$max_pheight - mean_pheight$min_pheight + 0.000001)) 
colnames(mean_pheight) = c("site-spp-yr","n2","min_pheight","max_pheight","mean_pheight","stdev_pheight",
                           "min_out_pheight","max_out_pheight","pheight_sum","pheight_normalized_mean")
data = merge(data,mean_pheight,by.x="site_species_year",by.y="site-spp-yr")
# tiny addition to avoid division by zero
data$normalized_pheight = ((data$plant_height - data$min_pheight) / 
                          (data$max_pheight - data$min_pheight + 0.000001))
data$normalized_centered_pheight = data$normalized_pheight - data$pheight_normalized_mean
data$centered_absolute_pheight = data$plant_height - data$mean_pheight

# and yeah, we can get rid of site-species-year combos with only one individual represented
data = data[!data$n==1,]


# get rid of lost factors
data$site = factor(data$site)
data$site_species = factor(data$site_species)
data$site_species_year = factor(data$site_species_year)
data$canopy_position = factor(data$canopy_position)
data$site = factor(data$site)
data$species = factor(data$species)
data$growth_form = factor(data$growth_form)
data$site_species = factor(data$site_species)
data$site_species_year = factor(data$site_species_year)



# how many rows were lost
nrows - dim(data)[1]

# fraction of lost rows
(nrows - dim(data)[1])/nrows



# error check for weird outliers in doy and plant height - put in outlier file
outliers = rbind(data[data$estimate < data$min_out_doy,],
                 rbind(data[data$estimate > data$max_out_doy,],
                       rbind(data[data$plant_height < data$min_out_pheight,],
                             data[data$plant_height > data$max_out_pheight,])))
write.csv(outliers,paste(pheno,"_outliers.csv",sep=""),quote=FALSE,row.names=FALSE)


# divide slope by 45 so we keep small numbers for all factors
data$normalized_slope = data$slope/45.0

# create column for disease_yn
data$disease = ifelse(is.na(data$disease),"",data$disease)
data$disease_yn = ifelse(data$disease == "",0,1)

# dummy variablize canopy position
# we'll drop the "full shade" factor, so we haven't overspecified the model
dummy = model.matrix(~canopy_position-1,data)
data = cbind(data,dummy)


# specifically map site-spp-yr to integers for reverse lookup
counts = as.data.frame(table(data$site_species_year))
counts$index = seq(1,dim(counts)[1])
counts$mean = tapply(data$estimate,data$site_species_year,mean)
counts$stdev = tapply(data$estimate,data$site_species_year,sd)
colnames(counts) = c("site-spp-yr","n","index","mean","stdev")
write.csv(counts,paste(pheno,"_site-spp-yr_lookup.csv",sep=""),quote=FALSE,row.names=FALSE)

# map site-spp to integers for reverse lookup -- maybe not necessary
#sitesp = as.data.frame(table(data$site_species))
#sitesp$index = seq(1,dim(sitesp)[1])
#sitesp$mean = tapply(data$estimate,data$site_species,mean)
#sitesp$stdev = tapply(data$estimate,data$site_species,sd)
#colnames(sitesp) = c("site-spp","n","index","mean","stdev")
#write.csv(sitesp,paste(pheno,"_site-spp_lookup.csv",sep=""),quote=FALSE,row.names=FALSE)

# join everything to get data frame for model
d1 = merge(data,counts,by.x="site_species_year",by.y="site-spp-yr")
names(d1)[names(d1) == 'index'] <- 'sitespyr_index'
#d2 = merge(d1,sitesp,by.x="site_species",by.y="site-spp")
#names(d2)[names(d2) == 'index'] <- 'sitesp_index'
d2=d1

# output the data
write.csv(d1,paste(pheno,"_all_data.csv",sep=""),quote=TRUE,row.names=FALSE)


# set up and run the model
model_dat = list(N = nrow(data),
		             K = nrow(counts),
#                 L = nrow(sitesp),
                 sitespyr = d2$sitespyr_index,
#                 sitesp = d2$sitesp_index,
                 y = d2$estimate,
                 elev = d2$centered_elev,
                 slope = d2$slope,
                 aspectNS = d2$aspect_NS,
                 aspectEW = d2$aspect_EW,
#                 pheight = d2$normalized_centered_pheight,
                 pheight = d2$centered_absolute_pheight,
		             disease = d2$disease_yn,
		             canopy_full_sun = d2$canopy_position_full_sun,
		             canopy_mostly_shaded = d2$canopy_position_mostly_shaded,
		             canopy_open_grown = d2$canopy_position_open_grown,
		             canopy_partially_shaded = d2$canopy_position_partially_shaded)


fit = stan(file = "all-factors.stan", 
           data = model_dat, 
           iter = 1000, chains = 4,
           control = list(adapt_delta = 0.99))


#list_of_draws <- extract(fit)
#print(names(list_of_draws))

#pairs(fit, pars = c("beta_elev[1]", "beta_elev[9]", "beta_elev[13]", "tau", "yhat[1]", "lp__"), las = 1)
#pairs(fit, pars = c("beta_elev[1]", "beta_elev[10]", "sigma", "gamma", "tau", "yhat[1]", "lp__"), las = 1)


# analyze
#print(fit)
#print(fit,pars=c("beta_pheight"))


# plot results to file
plotfilename = paste(paste("gamma_values_plants_and_factors_",pheno,sep=""),".png",sep="")
png(plotfilename,width=900, height=600, units="px")
plot(fit,plotfun = "stan_dens",
             pars=c("gamma_elev",
                    "gamma_slope",
                    "gamma_aspectNS",
                    "gamma_aspectEW",
                    "gamma_slope_aspectNS",
                    "gamma_slope_aspectEW",
                    "gamma_pheight",
                    "gamma_disease",
                    "gamma_canopy_full_sun",
                    "gamma_canopy_mostly_shaded",
                    "gamma_canopy_open_grown",
                    "gamma_canopy_partially_shaded"))       
dev.off()


#plot(fit,pars=c("beta_pheight"))


fit_summary = summary(fit,pars=c("alpha",
                                 "beta_elev",
                                 "beta_slope",
                                 "beta_aspectNS",
                                 "beta_aspectEW",
                                 "beta_slope_aspectNS",
                                 "beta_slope_aspectEW",
                                 "beta_pheight",
                                 "beta_disease",
                                 "beta_canopy_full_sun",
                                 "beta_canopy_mostly_shaded",
                                 "beta_canopy_open_grown",
                                 "beta_canopy_partially_shaded",
                                 "sigma",
                                 "gamma_elev",
                                 "gamma_slope",
                                 "gamma_aspectNS",
                                 "gamma_aspectEW",
                                 "gamma_slope_aspectNS",
                                 "gamma_slope_aspectEW",
                                 "gamma_pheight",
                                 "gamma_disease",
                                 "gamma_canopy_full_sun",
                                 "gamma_canopy_mostly_shaded",
                                 "gamma_canopy_open_grown",
                                 "gamma_canopy_partially_shaded",       
                                 "tau_elev",
                                 "tau_slope",
                                 "tau_aspectNS",
                                 "tau_aspectEW",
                                 "tau_slope_aspectNS",
                                 "tau_slope_aspectEW",
                                 "tau_pheight",
                                 "tau_disease",
                                 "tau_canopy_full_sun",
                                 "tau_canopy_mostly_shaded",
                                 "tau_canopy_open_grown",
                                 "tau_canopy_partially_shaded"))
#print(fit_summary)
                     

# write results to file
outfilename = paste(paste("parameter_values_plants_and_factors_",pheno,sep=""),".csv",sep="")
write.csv(as.data.frame(fit_summary),outfilename,quote=FALSE)



            
                 
