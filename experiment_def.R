#################################################################################################################
#                                                                                                               #
# R script to conduct analyses for manipulative experimental data                                               # 
#                                                                                                               #
# Fungal Biodiversity Mediates the Effects of Drying on Freshwater Ecosystem Functioning                        #
#                                                                                                               #
# Rebeca Arias-Real, Cayetano Gutiérrez-Cánovas, Isabel Muñoz, Cláudia Pascoal, Margarita Menéndez              #
#                                                                                                               #
# Code written by Rebeca Arias-Real and Cayetano Gutiérrez-Cánovas                                              #
# email for queries: rebeca.arias.real@gmail.com or cayeguti@um.es                                              #
#################################################################################################################

# Set your working directory
setwd("")

# Loading libraries
library(lme4)
library(car)
library(MuMIn)
library(multcomp)
library(sandwich)
library(lmerTest)
library(lattice)
library(variancePartition)
library(nlme)
library(gtools)
library(plyr)

# Loading data
dat<-read.table("experimental_dat.txt",h=T,sep="\t", dec=".")
dat_m<-read.table("experimental_dat.txt",h=T,sep="\t", dec=".")

# Removing NAs
na.omit(dat)->dat

# Calculating mean and SD of response variables

k_m<-ddply(dat_m, .(ric, treat), function(x) mean(1-x$k, na.rm =T))
k_sd<-ddply(dat_m, .(ric, treat), function(x) sd(1-x$k, na.rm =T))
k_n<-ddply(dat_m, .(ric, treat), function(x) length(na.omit(x$k)))

fb_m<-ddply(dat_m, .(ric, treat), function(x) mean(x$fb, na.rm =T))
fb_sd<-ddply(dat_m, .(ric, treat), function(x) sd(x$fb, na.rm =T))
fb_n<-ddply(dat_m, .(ric, treat), function(x) length(na.omit(x$fb)))

prod_m<-ddply(dat_m, .(ric, treat), function(x) mean(x$prod, na.rm =T))
prod_sd<-ddply(dat_m, .(ric, treat), function(x) sd(x$prod, na.rm =T))
prod_n<-ddply(dat_m, .(ric, treat), function(x) length(na.omit(x$prod)))

fun_var<-data.frame(k_m, k_sd=k_sd$V1, k_n=k_n$V1, 
                    fb_m=fb_m$V1, fb_sd=fb_sd$V1, fb_n=fb_n$V1,
                    prod_m=prod_m$V1, prod_sd=prod_sd$V1, prod_n=prod_n$V1)

write.table(fun_var,"fun_var.txt", sep="\t")

# Exploring data distributions

par(mfrow=c(1,2))
for (i in c(8,9)) hist(dat[,i], main=names(dat)[i])
par(mfrow=c(1,1))

# Transforming fungal biomass data
dat$fb<-log(dat$fb)

# Exploring responses of experimental manipulation of richness and 
# drying duration

###############################
# Organic matter decomposition
###############################

# Checking the necessity of random factor
mod1<-gls (k~Ric*treat*Time, method="REML", data=dat) 

mod2<-lme (k~Ric*treat*Time, random=~1|exp_unit, method="REML", data=dat)

AICc(mod1,mod2) # Random factor is not required

# Global model
mod<-lm(k~Ric*treat*Time, data=dat) 

# Checking residuals on global model
par(mfrow=c(1,2))
plot(fitted(mod),resid(mod),main="Residual homocedasticity")
abline(h=0,col="blue")
hist(resid(mod),main="Residual normality")
par(mfrow=c(1,1))

# Main test
Anova(mod, type=3)

options(na.action = "na.fail") # necessary to run dredge()

# runs all possible models 	and ranks the output according to the AIC, R2 is printed 
mod_d <- dredge (mod, rank = "AICc", extra=c("R^2")) 
mod_d

# Retaining models with AICc less or equal to 2
mod_set <- get.models (mod_d, subset=delta<= 2) # subset 

# Main test for best model
bmod<- mod_set[[1]]
summary(bmod)
anova(bmod)

# Checking residuals on the best-fitting model
par(mfrow=c(1,2))
plot(fitted(bmod),resid(bmod),main="Residual homocedasticity")
abline(h=0,col="blue")
hist(resid(bmod),main="Residual normality")
par(mfrow=c(1,1))

# Variance partitioning
round(variancePartition::calcVarPart(bmod),3)

# Post-hoc for final model
amod<-aov(k~Ric+treat+Time+Ric:Time+Time:treat,data=dat)

no_t1 <- which(dat$Time=="T2")

amod<-aov(dat$k[no_t1]~dat$Ric[no_t1]+dat$treat[no_t1] )
summary(amod)

# post-hoc anaysis
posthoc<-TukeyHSD(amod, conf.level=0.95)
posthoc

###############################
# Fungal biomass
###############################

# Checking the necessity of random factor
mod1<-gls (fb~Ric*treat*Time, method="REML", data=dat) 

mod2<-lme (fb~Ric*treat*Time, random=~1|exp_unit, method="REML", data=dat)

AICc(mod1,mod2) # Random factor is not required

# Global model
mod<-lm(fb~Ric*treat*Time, data=dat) 

# Checking residuals on global model
par(mfrow=c(1,2))
plot(fitted(mod),resid(mod),main="Residual homocedasticity")
abline(h=0,col="blue")
hist(resid(mod),main="Residual normality")
par(mfrow=c(1,1))

options(na.action = "na.fail") # necessary to run dredge()

# runs all possible models 	and ranks the output according to the AIC, R2 is printed 
mod_d <- dredge (mod, rank = "AICc", extra=c("R^2")) 
View(mod_d)

# Retaining models with AICc less or equal to 2
mod_set <- get.models (mod_d, subset=delta<= 2) # subset 

# Main test for best model
bmod<- mod_set[[1]]
summary(bmod)
anova(bmod)

# Checking residuals on the best-fitting model
par(mfrow=c(1,2))
plot(fitted(bmod),resid(bmod),main="Residual homocedasticity")
abline(h=0,col="blue")
hist(resid(bmod),main="Residual normality")
par(mfrow=c(1,1))

# Variance partitioning
round(variancePartition::calcVarPart(bmod),3)

# Final model for post-hoc
amod<-aov(fb~Ric*treat*Time,data=dat)

# post-hoc anaysis
posthoc<-TukeyHSD(amod, conf.level=0.95)
posthoc

