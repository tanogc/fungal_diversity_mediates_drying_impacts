#################################################################################################################
#                                                                                                               #
# R script to conduct analyses for field data                                                                   # 
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

# libraries
library(MuMIn)
library(variancePartition)
library(plyr)
library(vegan)
library(AER)
library(betapart)
library(ecodist)
library(ade4)
library(adespatial)
library(spdep)
library(usdm)
library(piecewiseSEM)

#loading data
dat<-read.table("fungi_dat.txt",h=T,sep="\t", dec=".")

# Arranging data
rownames(dat)<-dat$site
dat<-dat[,-c(1)]

q<-sapply(dat,class)=="numeric" | sapply(dat,class)=="integer"# 

par(mfrow=c(2,2))
for (i in which(q==T)) hist(dat[,i], main=names(dat)[i])
par(mfrow=c(1,1))

dat$ZFL<-(sqrt(dat$ZFL))
dat$k<-(sqrt(dat$k))
dat$srp<-log(dat$srp)
dat$din<-log(dat$din)
dat$area<-log(dat$area)
log(dat$FB)->dat$FB

par(mfrow=c(1,1))

# Identifying non-defined numbers
apply(dat,2,function(x) any(is.infinite(x)==T))

# Exploring which environmental variable best correlates with fungal richness

cor.res<-t(round(cor(dat$ric, dat[,8:ncol(dat)]),2))
write.table(cor.res,"cor_res.txt",sep="\t")
cor.test(dat$ric, dat$alt)

# Creating a standardized dataset

dat[, c("alt","ZFL","ZFP", "RE","ZFT")]<-scale(dat[, c("alt","ZFL","ZFP", "RE","ZFT")])

# Collinearity

# Pairwise correlations
round(as.dist(cor(dat[, c("alt","ZFL","ZFP", "RE","ZFT")] , method=c("pearson"))),2)

# VIF
vifstep(dat[, c("alt","ZFT")])
vifstep(dat[, c("alt","ZFP")])
vifstep(dat[, c("alt","ZFL")])

# Models

mod0 <- lm(ric ~ alt , data=dat)

mod1 <- lm(ric ~ alt + ZFT, data=dat)

mod2 <- lm(ric ~ alt + ZFP, data=dat)

mod3 <- lm(ric ~ alt + ZFL, data=dat)

mod4 <- lm(ric ~ ZFT, data=dat)

mod5 <- lm(ric ~ ZFP, data=dat)

mod6 <- lm(ric ~ ZFL, data=dat)

mod.list<-list(mod0,mod1, mod2, mod3, mod4, mod5, mod6)

model.sel(mod.list, extra=c(r2=function(x) r.squaredGLMM(x)[1]))

# Variance partitioning
round(calcVarPart(mod1),3)

#############
# Exploring changes in beta-diversity
#############

# Loading fungal species composition
spp<-read.table("fungal_spp.txt",h=T,sep="\t", dec=".")

## Calculating the contribution of species turnover and nestedness to explain
## community changes in response to drought intensity

# Calculating overall, turnover and nestedness-resultant dissimilarities
betapart.core(spp [,-c(1)])->spp.beta
beta.pair(spp.beta, index.family="sor")->spp.pair

beta.multi(spp[,-c(1)])->pi

# Euclidean distance for annual drying duration
dry.d<-dist(scale(dat$ZFT))

# Euclidean distance (altitude)
env.d<-dist(scale(dat$alt))

# Creating spatial structure
xy <- as.matrix(cbind(dat$x, dat$y))
nb1 <- graph2nb(gabrielneigh(xy), sym = T)
lw1 <- nb2listw(nb1)
w <- dist(xy)

f.sor <- spp.pair$beta.sor
f.sim <- spp.pair$beta.sim
f.sne <- spp.pair$beta.sne

# Not significant relationship
env.sor <- mantel.randtest(quasieuclid(f.sor), env.d)
env.sor

# Relationship between beta-diversity components and dry.d
m.sor <- mantel.randtest(quasieuclid(f.sor), dry.d)
m.sim <- mantel.randtest(quasieuclid(f.sim), dry.d)
m.sne <- mantel.randtest(quasieuclid(f.sne), dry.d)

# Removing potential spatial auto-correlation
msr.sor <- msr(m.sor, lw1, 999)  
msr.sim <- msr(m.sim, lw1, 999) 
msr.sne <- msr(m.sne, lw1, 999) 
 
# value of the statistic corrected for spatial autocorrelation
c(as.numeric(msr.sor$obs - msr.sor$expvar["Expectation"]),
  as.numeric(msr.sim$obs - msr.sim$expvar["Expectation"]),
  as.numeric(msr.sne$obs - msr.sne$expvar["Expectation"]))->r

# p-value corrected for spatial autocorrelation
c(msr.sor$pvalue, msr.sim$pvalue, msr.sne$pvalue)->pval

# results
beta_dry <- cbind(r, pval)
rownames(beta_dry) <- c("Sorensen","Turnover", "Nestedness")

# Plotting results
pdf(file="beta_resp.pdf",onefile=T,width=10,height=4)
par(mfrow=c(1,3),cex.lab=1.5,cex.axis=1.5,mar=c(5,5,5,1),pch=16)
plot(spp.pair$beta.sor~dry.d,ylab="Sorensen dissimilarity",xlab="hydrological distance",main="Overall")
abline(lm(as.numeric(f.sor)~as.numeric(dry.d)),col="blue",lwd=2)
plot(spp.pair$beta.sim~dry.d,ylab="Simpson dissimilarity",xlab="hydrological distance",main="Turnover")
abline(lm(as.numeric(f.sim)~as.numeric(dry.d)),col="blue",lwd=2)
plot(spp.pair$beta.sne~dry.d,ylab="Nestedness-resultant dissimilarity",xlab="hydrological distance",main="Nestedness")

dev.off()

# Exporting beta-diversity results
write.table(beta_dry, "beta_dry.txt", sep="\t")

# Structural Equation Modelling (SEM)

# Creating SEM dataframe

sem_set <- data.frame(k=scale(dat$k),
                      FB=scale(dat$FB),
                      ric=scale(dat$ric), 
                      ZFT=scale(dat$ZFT),
                      ZFP=scale(dat$ZFP),
                      ZFL=scale(dat$ZFL),
                      alt=scale(dat$alt))

# Direct relationships between fungal richness and ecosystem processes
summary(lm(k~ric, sem_set))
summary(lm(FB~ric, sem_set))

####################################
# SEM with multiple paths ####
####################################

# Creating a data frame to store SEM results
sem_res <- data.frame(process=rep(c("decomp", "fungal biomass"), each=12), matrix(NA, nrow=3*4*2, 12))

res_names <- c("Process","Response","Predictor","Std.Estimate","Std.Error","DF","Crit.Value","P.Value", 
  "r2","AICc", "Fisher.C", "df", "pval")

colnames(sem_res) <- res_names 

# Producing SEM results

for (i in 1:4) {
  
  a=i+3
  
  pred=sem_set[,a]

mod <- lm(k ~ pred, data.frame(sem_set,pred))

sem <- psem(lm(ric ~ pred, data.frame(sem_set,pred)), lm(k ~ ric, data.frame(sem_set,pred)))
sem_summary <-summary(sem, .progressBar = F)


(((i-1)*3)+1):(((i-1)*3)+2)

sem_summary$coefficients[-c(8,9)] -> sem_res[(((i-1)*3)+1):(((i-1)*3)+2),2:8]
c(rsquared(sem)$R.squared,summary(mod)$adj.r) -> sem_res[(((i-1)*3)+1):(i*3),9]
sem_summary$IC[2] -> sem_res[(((i-1)*3)+1),10]
sem_summary$Cstat -> sem_res[(((i-1)*3)+1),11:13]
c(coef(summary(mod))[2,1:2],13, coef(summary(mod))[2,3:4])-> sem_res[(i*3),4:8]

colnames(sem_set)[a] -> sem_res[(((i-1)*3)+1),3] -> sem_res[(i*3),3]
"k"-> sem_res[(i*3),2]
}

for (i in 5:8) {
  
  a=i-1
  
  pred=sem_set[,a]
  
  mod <- lm(FB ~ pred, data.frame(sem_set,pred))
  
  sem <- psem(lm(ric ~ pred, data.frame(sem_set,pred)), lm(FB ~ ric, data.frame(sem_set,pred)))
  sem_summary <-summary(sem, .progressBar = F)

  sem_summary$coefficients[-c(8,9)] -> sem_res[(((i-1)*3)+1):(((i-1)*3)+2),2:8]
  c(rsquared(sem)$R.squared,summary(mod)$adj.r) -> sem_res[(((i-1)*3)+1):(i*3),9]
  sem_summary$IC[2] -> sem_res[(((i-1)*3)+1),10]
  sem_summary$Cstat -> sem_res[(((i-1)*3)+1),11:13]
  c(coef(summary(mod))[2,1:2],13, coef(summary(mod))[2,3:4])-> sem_res[(i*3),4:8]
  
  colnames(sem_set)[a] -> sem_res[(((i-1)*3)+1),3] -> sem_res[(i*3),3]
  "FB"-> sem_res[(i*3),2]
}

# Saving SEM results
write.table(sem_res, "sem_res.txt", sep="\t")


# Effect sizes for organic matter decomposition SEMs

# Richness
ric_ef=sem_res[2,4]
ric_se=sem_res[2,5]

ric_ci1<-ric_ef - qnorm(0.975)*as.numeric(ric_se)
ric_ci2<-ric_ef + qnorm(0.975)*as.numeric(ric_se)

# ZFT
ZFT_ef=sem_res[3,4]
ZFT_se=sem_res[3,5]

ZFT_ci1<-ZFT_ef - qnorm(0.975)*as.numeric(ZFT_se)
ZFT_ci2<-ZFT_ef + qnorm(0.975)*as.numeric(ZFT_se)

# ZFP
ZFP_ef=sem_res[6,4]
ZFP_se=sem_res[6,5]

ZFP_ci1<-ZFP_ef - qnorm(0.975)*as.numeric(ZFP_se)
ZFP_ci2<-ZFP_ef + qnorm(0.975)*as.numeric(ZFP_se)

# ZFL
ZFL_ef=sem_res[9,4]
ZFL_se=sem_res[9,5]

ZFL_ci1<-ZFL_ef - qnorm(0.975)*as.numeric(ZFL_se)
ZFL_ci2<-ZFL_ef + qnorm(0.975)*as.numeric(ZFL_se)

# Alt
alt_ef=sem_res[12,4]
alt_se=sem_res[12,5]

alt_ci1<-alt_ef - qnorm(0.975)*as.numeric(alt_se)
alt_ci2<-alt_ef + qnorm(0.975)*as.numeric(alt_se)

# Aggregating data
m=c(ric_ef, ZFT_ef, ZFP_ef, ZFL_ef, alt_ef)
ci1=c(ric_ci1, ZFT_ci1, ZFP_ci1, ZFL_ci1, alt_ci1)
ci2=c(ric_ci2, ZFT_ci2, ZFP_ci2, ZFL_ci2, alt_ci2)

names(m)<-c("Ric", "ZFT", "ZFP", "ZFL", "Alt")

# Setting colours for each hydrologic metric
col.hyd<-c("purple","#56B4E9", "#7AD151FF","brown1", "gold3")

# Effect sizes for organic matter decomposition SEMs

pdf(file="decomp_sem.pdf",,useDingbats=FALSE,onefile=T,width=5.5,height=5.5)

dotplot(m[5:1], xlab="Standardized effect", xlim=c(-1.2, 1.2),
        
        par.settings = list(axis.text = list(cex = 2, font=1), 
                            par.xlab.text = list(cex = 2, font=1), 
                            par.ylab.text = list(cex = 2, font=1)),
        panel=function(x,y){
          panel.segments(ci1[5:1], as.numeric(y), ci2[5:1], as.numeric(y), lty=1, col=col.hyd[5:1])
          panel.xyplot(x, y, pch=15, cex=2,col=col.hyd[5:1])
          panel.abline(v=0, col="black", lty=2)
          
        })

dev.off()


## Effect sizes for fungal biomass SEMs

# Richness
ric_ef=sem_res[14,4]
ric_se=sem_res[14,5]

ric_ci1<-ric_ef - qnorm(0.975)*as.numeric(ric_se)
ric_ci2<-ric_ef + qnorm(0.975)*as.numeric(ric_se)

# ZFT
ZFT_ef=sem_res[15,4]
ZFT_se=sem_res[15,5]

ZFT_ci1<-ZFT_ef - qnorm(0.975)*as.numeric(ZFT_se)
ZFT_ci2<-ZFT_ef + qnorm(0.975)*as.numeric(ZFT_se)

# ZFP
ZFP_ef=sem_res[18,4]
ZFP_se=sem_res[18,5]

ZFP_ci1<-ZFP_ef - qnorm(0.975)*as.numeric(ZFP_se)
ZFP_ci2<-ZFP_ef + qnorm(0.975)*as.numeric(ZFP_se)

# ZFL
ZFL_ef=sem_res[21,4]
ZFL_se=sem_res[21,5]

ZFL_ci1<-ZFL_ef - qnorm(0.975)*as.numeric(ZFL_se)
ZFL_ci2<-ZFL_ef + qnorm(0.975)*as.numeric(ZFL_se)

# Alt
alt_ef=sem_res[24,4]
alt_se=sem_res[24,5]

alt_ci1<-alt_ef - qnorm(0.975)*as.numeric(alt_se)
alt_ci2<-alt_ef + qnorm(0.975)*as.numeric(alt_se)


# Aggregating data
m=c(ric_ef, ZFT_ef, ZFP_ef, ZFL_ef, alt_ef)
ci1=c(ric_ci1, ZFT_ci1, ZFP_ci1, ZFL_ci1, alt_ci1)
ci2=c(ric_ci2, ZFT_ci2, ZFP_ci2, ZFL_ci2, alt_ci2)

names(m)<-c("Ric", "ZFT", "ZFP", "ZFL", "Alt")

# Setting colours for each hydrologic metric
col.hyd<-c("purple","#56B4E9", "#7AD151FF","brown1", "gold3")

# Effect sizes for fungal biomass SEMs

pdf(file="fb_sem.pdf",,useDingbats=FALSE,onefile=T,width=5.5,height=5.5)

dotplot(m[5:1], xlab="Standardized effect", xlim=c(-1.2, 1.2),
        
        par.settings = list(axis.text = list(cex = 2, font=1), 
                            par.xlab.text = list(cex = 2, font=1), 
                            par.ylab.text = list(cex = 2, font=1)),
        panel=function(x,y){
          panel.segments(ci1[5:1], as.numeric(y), ci2[5:1], as.numeric(y), lty=1, col=col.hyd[5:1])
          panel.xyplot(x, y, pch=15, cex=2,col=col.hyd[5:1])
          panel.abline(v=0, col="black", lty=2)
          
        })

dev.off()

# Barplot showing explained variance for fungal richness and drying characteristics

# Preparing data
k.r2=c(sem_res[2,9], sem_res[3,9], sem_res[6,9], sem_res[9,9], sem_res[12,9])
fb.r2=c(sem_res[14,9], sem_res[15,9], sem_res[18,9], sem_res[21,9], 0)

names(k.r2)<-names(fb.r2)<-c("Ric", "ZFT", "ZFP", "ZFL", "Alt")

# Barplot

pdf(file="barplot_sem.pdf",,useDingbats=FALSE,onefile=T,width=12.5,height=5.5)

par(mfrow=c(1,2),cex.lab=1.9, cex.axis=1.75, mar=c(4,5,4,4))

barplot(k.r2*100, ylim=c(0,60), main="", col=col.hyd,ylab="Explained variance (%)",xlab="")
abline(h=sem_res[2,9]*100, lty=2)

barplot(fb.r2*100, ylim=c(0,60), main="", col=col.hyd,ylab="",xlab="")
abline(h=sem_res[14,9]*100, lty=2)

dev.off()

