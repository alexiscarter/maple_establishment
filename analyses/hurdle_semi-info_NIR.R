### Modelling of survival and biomass of sugar maple seedlings 
# Non-irradiated soil only
# semi-informativepriors for betas for Gamma and Bernoulli

## Libraries and data ####
library(plyr)
library(tidyverse)
library(rjags)
library(R2jags)
library(lattice)
library(ggplot2)
library(ggpubr)

load("data/traits_as.rdata")

# Manipulation ####
Use<-c("id" ,"species", "forest", "block" ,"treatment", "total.w")
trait_as<-traits_as[, names(traits_as) %in% Use]

# For modeling with survival data, create a column with O/1
trait_as$surv <- trait_as$total.w
trait_as$surv[which(trait_as$surv > 0)] <- 1

# Keep only untreated
trait_as_nir <- trait_as %>%
  filter(treatment == 'NIR') %>% 
  droplevels()

# change names of treatment covariables
trait_as_nir$forest <-revalue(trait_as_nir$forest, c("FT"="Temperate", "FM"="Mixed",'FB'= 'Boreal'))

## Hurdle model on acer data   #####
# Create a covariable matrix for the biomasse data
Xc<-model.matrix(~ forest, data= trait_as_nir)

# Create a covariable matrix for the survival
Xb<-model.matrix(~ forest, data= trait_as_nir)

# Define the number of parameter biomass part
Kc=ncol(Xc)

# Define the number of parameter survival part
Kb=ncol(Xb)

# Create a random effect id
block<- as.numeric(as.factor(trait_as_nir$block))

# Define the number of random effects
Nre<-length(unique(block))

# List all the required variables for JAGS
JAGS.data_full<-list(
  Y     =  trait_as_nir$total.w,
  Xc    =  Xc,
  Xb    =  Xb,
  Kc    =  Kc,
  Kb    =  Kb,
  N     = nrow(trait_as_nir),
  Block = block,
  Nre   = Nre,
  Zeros = rep(0, nrow(trait_as_nir)))

# Write down the model
load.module("glm")
sink("ZA_LASSO.txt")
cat("
    model{
    #1A Priors beta and Gamma
    for(i in 1:Kc)  { beta[i]  ~ dnorm(0, 0.333333)}
    for(i in 1:Kb)  { gamma[i] ~ dnorm(0, 0.333333)}
    
    #1B Priors for r parameter of gamma distribution
    r~ dunif(0,5)
    
    #1B Priors random effect
    for(i in 1:Nre){
    a1[i] ~ dnorm(0, tau1_Block)
    }
    
    #1C Priors random effect
    for(i in 1:Nre){
    a2[i] ~ dnorm(0, tau2_Block)
    }
    
    #1D Diffuse uniform prior for sigma_Block
    tau1_Block<-1/(sigma1_Block*sigma1_Block)
    sigma1_Block ~ dunif(0,100)
    
    #1E Diffuse uniform prior for sigma_Block
    tau2_Block<-1/(sigma2_Block*sigma2_Block)
    sigma2_Block ~ dunif(0,100)
    
    #2 likelihood ( zero trick)
    C <- 1000
    for( i in 1:N){
    Zeros[i] ~ dpois(-ll[i] + C)
    z[i] <- step(Y[i] - 0.0001)
    l1[i] <- (1 - z[i]) * log(1-Pi[i])
    l2[i] <- z[i] * ( log(Pi[i]) - loggam(r) + r*log(r/mu[i])+
    (r-1)*log(Y[i]) - (Y[i] * r)/mu[i] )
    
    ll[i]<- l1[i] + l2[i]
    
    log(mu[i]) <- inprod( beta[], Xc[i,]) + a1[Block[i]]
    logit(Pi[i]) <- inprod( gamma[], Xb[i,]) + a2[Block[i]]
    }
    
    #3 Discrepancy measures
    for(i in 1:N){
    ExpY[i] <- Pi[i] * mu[i]
    VarY[i] <- (Pi[i] * r + Pi[i] - Pi[i] *
    Pi[i]*r) * (mu[i]+mu[i]/r)
    
    PRes[i] <- (Y[i] - ExpY[i])/ sqrt(VarY[i])
    }
    }", fill = TRUE)
sink()

# Specify the initiale values
inits<- function() {
  list(beta  = rnorm(Kc, 0, 0.1),
       gamma = rnorm(Kb, 0, 0.1),
       r = runif(1,0,5),
       a1 = rnorm(Nre, 0,0.1),
       sigma1_Block = runif(1,0.001,5),
       a2 = rnorm(Nre, 0,0.1),
       sigma2_Block = runif(1,0.001,5))  
}

# Specify the parameter to store
params<-c("beta", "gamma", "sigma1_Block", "sigma2_Block",
          "PRes", "r", "ExpY")

# Run the hurdle model
nir.semi.info <- jags( data        = JAGS.data_full,
                       inits       = inits,
                       parameters  = params,
                       model       = "ZA_LASSO.txt",
                       n.thin      = 10,
                       n.chains    = 3,
                       n.burnin    = 4000,
                       n.iter      = 5000)

# /!\ Can be long to run
nir.semi.info2 <- update(nir.semi.info, n.iter = 500000, n.thin = 10)

# Extract info
out.nir <-nir.semi.info2$BUGSoutput

# check if the mixing is good ####
MyNames<- c(
  paste(c(colnames(Xc), "sigma1_Block") , "Gamma", sep = " "),
  paste(c(colnames(Xb), "sigma2_Block") , "Bern" , sep = " "),
  "r Gamma")

MyBUGSChains(out.nir, # run MCMCSupportHighstatV2.R before
             c(uNames("beta",Kc), "sigma1_Block",
               uNames("gamma",Kb), "sigma2_Block",
               "r"),
             PanelNames = MyNames)

# get the numerical output of JAGS  ####
OUT.nir <- MyBUGSOutput(out.nir,  c(uNames("beta",Kc), "sigma1_Block",
                                    uNames("gamma",Kb),"sigma2_Block",
                                    "r"),
                        VarNames = MyNames)
print(OUT.nir, digits = 5)

# get the posterior distribution ####
# Gamma part of the model
MyNamesG<- c(colnames(Xc), "sigma1_Block","r")

MyBUGSHist(out.nir,
           c(uNames("beta",Kc), "sigma1_Block" ,"r"),
           PanelNames = MyNamesG)

# Bernoulli part of the model
MyNamesB<- c(colnames(Xb), "sigma2_Block")

MyBUGSHist(out.nir,
           c(uNames("gamma",Kb), "sigma2_Block"),
           PanelNames = MyNamesB)

## Get Pearson residuals and Expected values of the hurdle model ####
E <-out.nir$mean$PRes
ExpY <-out.nir$mean$ExpY

## calculate de dispersion of the model ####
D<- sum(E^2)/(nrow(trait_as_nir) + (ncol(Xc)+ ncol(Xb)))

## Model validation ####
theme_set(theme_bw())

data1<-data.frame(E,ExpY)
Pres_ExpV<-ggplot(data= data1, aes( x= ExpY, y=E))+
  geom_point()+
  geom_hline(yintercept = 0, linetype="dotted")+
  ylab('Pearson residuals')+
  xlab('Expected values')
Pres_ExpV

data3<-data.frame(trait_as_nir$forest, E)
Pres_bio<-ggplot(data= data3, aes( x= trait_as_nir$forest  , y=E))+
  geom_boxplot()+
  geom_hline(yintercept = 0, linetype="dotted")+
  ylab('Pearson residuals')+
  theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1)) +
  xlab('Forest type')
Pres_bio

data4<-data.frame(trait_as_nir$total.w, ExpY)
ObvsExp<-ggplot(data= data4, aes(y = ExpY, x=trait_as_nir$total.w ))+
  geom_point()+
  labs(x ='Observed biomass (g)', y='Expected biomass (g)')
ObvsExp

model_val<-ggarrange(Pres_ExpV, ObvsExp, Pres_bio ,
                     ncol=2, nrow=2, labels=c("a", "b", "c"))

## get the R2 of obs.value ~ expected.value
tt<-lm(ExpY~trait_as_nir$total.w)
summary(tt)$adj.r.squared

## sketch the model fit
MyData <- expand.grid(Forest = levels(trait_as_nir$forest))

X<-model.matrix(~Forest, data= MyData)
beta.mcmc  <- out.nir$sims.list$beta
gamma.mcmc <- out.nir$sims.list$gamma

mu.mcmc    <- exp(X %*% t(beta.mcmc))
Pi.mcmc    <- exp(X %*% t(gamma.mcmc)) / (1 + exp(X %*% t(gamma.mcmc)))
ExpY.mcmc <- Pi.mcmc * mu.mcmc

## Violin plot ####
## Hurdle posterior
## format data for ggplot and select 90%
ExpY.mcmc.df <- ExpY.mcmc %>%
  t() %>%
  as.data.frame() %>%
  dplyr::rename('Temperate'='1', 'Mixed'='2', 'Boreal'='3') %>%
  gather(key = Soil, factor_key = TRUE) %>%
  group_by(Soil) %>%
  filter(value <= quantile(value, 0.95) & value >= quantile(value, 0.05))

Hurdle<- ggplot(data = ExpY.mcmc.df, aes(x = Soil,  y = value)) +
  geom_violin()+
  stat_summary(fun.y=median, geom="point", shape=23, size=2)+
  theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1)) +
  labs(x = "Forest type", y = "Dry mass (including survival, g)")

## Gamma posterior
## format data for ggplot and select 90%
mu.mcmc.df <- mu.mcmc %>%
  t() %>%
  as.data.frame() %>%
  dplyr::rename('Temperate'='1', 'Mixed'='2', 'Boreal'='3') %>%
  gather(key = Soil, factor_key = TRUE) %>%
  group_by(Soil) %>%
  filter(value <= quantile(value, 0.95) & value >= quantile(value, 0.05))

Gamm<- ggplot(data = mu.mcmc.df, aes(x = Soil,  y = value)) +
  geom_violin()+
  stat_summary(fun.y=median, geom="point", shape=23, size=2)+
  theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1)) +
  labs(x = "", y = "Dry mass (of survivors, g)")

## Bernoulli posterior
## format data for ggplot and select 90%
Pi.mcmc.df <- Pi.mcmc %>%
  t() %>%
  as.data.frame() %>%
  dplyr::rename('Temperate'='1', 'Mixed'='2', 'Boreal'='3') %>%
  gather(key = Soil, factor_key = TRUE) %>%
  group_by(Soil) %>%
  filter(value <= quantile(value, 0.95) & value >= quantile(value, 0.05))

Bern<- ggplot(data = Pi.mcmc.df, aes(x = Soil,  y = value)) +
  geom_violin()+
  stat_summary(fun.y=median, geom="point", shape=23, size=2)+
  theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1)) +
  labs(x = "", y = "Survival probability")

plot.all<-ggarrange(Bern, Gamm, Hurdle,
                    ncol=1, nrow=3, labels=c("a", "b", "c")) 
