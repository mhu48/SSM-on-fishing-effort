rm(list = ls())
set.seed(141)
setwd("D:/Research/STAT262/JSM 2020/July meeting prep/codes")

#install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library("INLA")
library(fields)
library(dplyr)

library(rgdal)
library(rgeos)
library(sf)
library("pryr")
library(sp)
library(ggplot2)

# Reorganize the weekly data: find all occurences of locations, 
# and for weeks where the locations were not present, add in the location for the week with value 0
setwd("D:/Research/STAT262/JSM 2020")
data = read.csv("fe_2017_1_2017_12_west_weekly.csv", header = TRUE)  #Dropping the last couple of days exceeding the 52 weeks
# write.csv(tmp,"fe_2018_1_2018_12_east_weekly.csv")
data = data[,-1]

setwd("D:/Research/STAT262/JSM 2020/July meeting prep/codes")
locations = matrix(c(0,0), 1,2)
colnames(locations) = c("long", "lat")

for(i in c(1:length(unique(data$time)))) {
  tmp = filter(data, data$time == unique(data$time)[i])[,2:3]
  tmp_locations = do.call(rbind, mget(c("locations", "tmp")))
  locations = tmp_locations[!duplicated(tmp_locations),]
}

locations = locations[-1, ]
rownames(locations) = c(1:nrow(locations))
rm(tmp_locations, tmp)
newdata = data[1,]
for (i in c(1:length(unique(data$time)))) {
  weeklydata = filter(data, data$time == unique(data$time)[i])
  tmp = weeklydata[,2:3]
  tmp_locations = do.call(rbind, mget(c("tmp", "locations")))
  fill_locations = tmp_locations[!duplicated(tmp_locations),]
  fill_locations = fill_locations[(nrow(tmp)+1):nrow(fill_locations),]
  
  for (j in c(1:nrow(fill_locations))) {
      if (fill_locations[j, 2] < 30) {
        storage = matrix(cbind(as.character(unique(weeklydata$time)), fill_locations[j,1], fill_locations[j,2], 0,
          mean(filter(weeklydata, weeklydata$lat < 30)$sst) + rnorm(1, 0, sd(filter(weeklydata, weeklydata$lat < 30)$sst))),nrow = 1
          )
        colnames(storage) = colnames(newdata)
        newdata = rbind(newdata, storage)
      } else if (fill_locations[j, 2] >= 30 && fill_locations[j, 2] < 40) {
        storage = matrix(cbind(as.character(unique(weeklydata$time)), fill_locations[j,1], fill_locations[j,2], 0,
                               mean(filter(weeklydata, fill_locations[j, 2] >= 30 & fill_locations[j, 2] < 40)$sst) + rnorm(1, 0, sd(filter(weeklydata, fill_locations[j, 2] >= 30 && fill_locations[j, 2] < 40)$sst))),nrow = 1
        )
        colnames(storage) = colnames(newdata)
        newdata = rbind(newdata, storage)
        } else{
          storage = matrix(cbind(as.character(unique(weeklydata$time)), fill_locations[j,1], fill_locations[j,2], 0,
                                 mean(filter(weeklydata, weeklydata$lat > 40)$sst) + rnorm(1, 0, sd(filter(weeklydata, weeklydata$lat > 40)$sst))),nrow = 1
          )
          colnames(storage) = colnames(newdata)
          newdata = rbind(newdata, storage)
          }
  }
}
newdata = newdata[-1,]
newdata$long = as.numeric(newdata$long)
newdata$lat = as.numeric(newdata$lat)
newdata$fh = as.numeric(newdata$fh)
newdata$sst = as.numeric(newdata$sst)
str(newdata)
data = rbind(data, newdata)
data$time = as.Date(data$time, "%Y-%m-%d")
data = arrange(data, time)
str(data)  #2017east: 6307/53 = 119  #2017west: 3445/53 = 65
write.csv(data, "fe_2017_1_2017_12_west_weekly_filled.csv")
# save.image("D:/STAT262/JSM 2020/July meeting prep/codes/Starting_data_east_2017.RData")
rm(fill_locations, newdata, storage, tmp, tmp_locations,weeklydata)

data_saved = data
data$time = as.numeric(factor(data$time))
n_stations = nrow(data[data$time == unique(data$time)[1],])
n_data = nrow(data)
n_weeks = length(unique(data$time))
mean_covariates = mean(data$sst)
sd_covariates = sd(data$sst)
data$sst = (data$sst - mean_covariates)/(sd_covariates)

max.edge = (diff(range(data$long))/10)
bound.outer = (diff(range(data$lat))/3)

# system.time(mesh1 <- inla.mesh.create.helper(points=cbind(data$long,data$lat),
#                                             max.edge = c(1,5)*max.edge,
#                                             offset = c(max.edge, bound.outer),
#                                             min.angle=c(45, 60)))  #59s
system.time(mesh2 <- inla.mesh.2d(loc=cbind(data$long,data$lat),
                                 max.edge = c(1,5)*max.edge,
                                 cutoff = max.edge/5,
                                 offset = c(max.edge, bound.outer))) #0.61s works well for 2017 west
# plot(mesh1)
# points(data$long,data$lat,pch=20,cex=2, col=2)
plot(mesh2)
points(data$long,data$lat,pch=20,cex=1, col=2)

#What if we add boundary?
setwd("D:/Research/FrontiersMarine")
system.time(eez <- st_read('World_EEZ_v10/eez_v10.shp'))
setwd("D:/Research/STAT262/JSM 2020/July meeting prep/codes")
eez_usa <-  dplyr::filter(eez, ISO_Ter1 == "USA")
eez_usa_48 <-  dplyr::filter(eez, MRGID==8456)
system.time(eez_usa <- st_as_sf( rmapshaper::ms_simplify(input = as(eez_usa_48, 'Spatial')) ))
#3749.46
object_size(eez_usa_48) # ~4.19MB
eez_west = st_cast(eez_usa_48,"POLYGON")[2,]  #use 2 for west coast
# st_write(eez_usa_48, "eez_usa_48.shp")
# st_write(eez_east, "eez_east.shp")
eez_polygon = as_Spatial(eez_west)

system.time(mesh3 <- inla.mesh.2d(loc=cbind(data$long,data$lat),
                                  max.edge = c(3,15)*max.edge,
                                  cutoff = max.edge/5,
                                  offset = c(max.edge, bound.outer),
                                  boundary=eez_polygon))
plot(mesh3)
points(data$long,data$lat,pch=20,cex=0.8, col=2)
rm(eez, eez_usa, eez_polygon, mesh1)
# save.image("D:/STAT262/JSM 2020/July meeting prep/codes/Starting_data_east_2017_0723mesh.RData")

mesh = mesh3
spde = inla.spde2.matern(mesh=mesh, alpha=2)
A = inla.spde.make.A(mesh, loc=as.matrix(data[,c("long", "lat")]), group=data$time, n.group=n_weeks)
mesh3$n  #419  #676
field.indices = inla.spde.make.index("field", n.spde=spde$n.spde, n.group=n_weeks) 

stack = inla.stack(data=list(y = data$fh),A=list(A, 1),
                     effects=list(c(field.indices ,list(Intercept=1)),
                                  list(sst = data$sst)),
                     tag="est")
# stack <- inla.stack(data = list(y = coords.st$y), 
#                     A = list(A, 1), effects = list( s.index, list(Intercept = rep(1, nrow(coords.st))) ), tag = "est")
formula <- ( y~ -1 + Intercept + sst + f(field , model = spde, group = field.group, control.group = list(model="ar1")) )
system.time( result <- inla(formula, data = inla.stack.data(stack, spde=spde),
                           family = "gaussian", control.predictor = list(A=inla.stack.A(stack), compute=TRUE)) )
# 5548.71
# Analysis of results
# posterior estimates of covariate coefficients for sst
summary(result)

#Construct Table II
# marginal density of 1/(sigma2_epsilon), transformed to sigma2_epsilon for the Gaussian distribution of Y
sigma2eps_marg = inla.tmarginal(function(x) 1/x, result$marginals.hyperpar$"Precision for the Gaussian observations")
sigma2eps_m1 = inla.emarginal(function(x) x, sigma2eps_marg)  #expected value of sigma2_epsilon based on marginal density
# sigma2eps_m2 = inla.emarginal(function(x) x^2, sigma2eps_marg)  #second order moment of sigma2_epsilon based on marginal density
# sigma2eps_stdev = sqrt(sigma2eps_m2 - sigma2eps_m1^2) #std dev of sigma2_epsilon based on marginal density
# sigma2eps_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2eps_marg)  #quantiles of sigma2_epsilon based on marginal density

# marginal variance, sigma2gamma
result.field = inla.spde.result(result , "field", spde, do.transform=TRUE)
kappa_marg = inla.tmarginal(function(x) x, result.field$marginals.kappa[[1]])
kappa_m1 = inla.emarginal(function(x) x, kappa_marg)  #expected value of kappa based on marginal density
tau_marg = inla.tmarginal(function(x) x, result.field$marginals.tau[[1]])
tau_m1 = inla.emarginal(function(x) x, tau_marg)  #expected value of tau based on marginal density
sigma2gamma_m1 = 1/(4*pi*(kappa_m1^2)*tau_m1)

# empiricaly derived correlation range, rho
rho = sqrt(8)/kappa_m1
# posterior estimates of correlation of thi in state equation under AR(1), lambda
result$summary.hyperpar["GroupRho for field",]
lambda = result$summary.hyperpar["GroupRho for field",][1]

# Plot the posterior marginals
# result.field = inla.spde.result(result , "field", spde, do.transform=TRUE)
# inla.emarginal(function(x) x, result.field$marginals.range.nominal[[1]])
# inla.emarginal(function(x) x, result.field$marginals.variance.nominal[[1]])
par(mfrow = c(2, 2), mar = c(3, 3, 1, 0.1), mgp = 2:0)
plot(result$marginals.hyper$`GroupRho for field`, type = "l", xlab = "beta", ylab = "Density")
abline(v = lambda, col = 2)
plot(result.field$marginals.variance.nominal[[1]], type = "l", xlab = "variance nominal", ylab = "Density")
abline(v = 18200, col = 2)
plot(result.field$marginals.kappa[[1]], type = "l", xlab = "kappa", ylab = "Density")
abline(v = kappa_m1, col = 2)
plot(result.field$marginals.range.nominal[[1]], type = "l", xlab = "range nominal", ylab = "Density")
abline(v = rho, col = 2)
par(mfrow = c(1,1))

# Plot the posterior of the random field
str(idat <- inla.stack.index(stack, "est")$data)
cor(data$fh, result$summary.linear.predictor$mean[idat])
stepsize <- 4 * 1/111
coords.sp <- SpatialPolygons(list(Polygons(list(Polygon(cbind(x = locations$long, y = locations$lat))), ID = "1")))
coords.sp <- spsample(coords.sp, n = 120, type = "regular")
head(coords.sp)
class(coords.sp)
str(coords.sp)
projgrid <- inla.mesh.projector(mesh, loc = as.matrix(as.data.frame(coords.sp)))
# xmean <- list()
newdata <- NULL
for (j in c(1,5,9,13,17,21,25,29,33,37,41,45,49,53)){
  newdata = rbind(newdata, data.frame(time = j, coords.sp, mean = inla.mesh.project(mesh, loc = as.matrix(as.data.frame(coords.sp)),
                                                                                    result$summary.random$field$mean[field.indices$field.group == j] 
                                                                                    + result$summary.fixed$mean[1])))
}

head(newdata)
ggplot(newdata, aes(y = x2, x = x1)) + 
  geom_tile(aes(fill = mean)) + 
  facet_wrap(~time) + 
  scale_fill_gradientn(colours = terrain.colors(10)) +
  theme_bw() +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()  #remove y axis ticks
  )
  















