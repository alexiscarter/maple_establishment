### Modelling of survival and biomass ofsugar maple seedlings 
# Irradiated soil only
# semi-informativepriors for betas for Gamma and Bernoulli

## Libraries and data ####
library(plyr)
library(tidyverse)
library(rjags)
library(R2jags)
library(lattice) # using in the MyBUGSChains function
library(ggplot2)
library(ggpubr)

load("data/traits_as.rdata")

# Manipulation #### 
Use<-c("id" ,"species", "forest", "block" ,"treatment", "total.w")
trait_as<-traits_as[, names(traits_as) %in% Use]

# For modeling with survival data, create a column with O/1
trait_as$surv <- trait_as$total.w
trait_as$surv[which(trait_as$surv > 0)] <- 1

# Remove untreated
trait_as_ir <- trait_as %>%
  filter(!treatment == 'NIR') %>% 
  droplevels()

## change names of treatment covariables
trait_as_ir$treatment<-revalue(trait_as_ir$treatment, c("IFT"="Temperate inoculant",
                                                  'IFM'= 'Mixed inoculant', 'IFB'= 'Boreal inoculant', 'SI'= 'Sterile'))
trait_as_ir$forest <-revalue(trait_as_ir$forest, c("FT"="Temperate", "FM"="Mixed",'FB'= 'Boreal'))

# Create a covariable matrix for the biomasse data
Xc<-model.matrix(~ (forest + treatment)^2 , data= trait_as_ir)

# Create a covariable matrix for the survival
Xb<-model.matrix(~ (forest + treatment)^2 , data= trait_as_ir)

# Define the number of parameter biomass part
Kc=ncol(Xc)

# Define the number of parameter survival part
Kb=ncol(Xb)

# Create a random effect id
block<- as.numeric(as.factor(trait_as_ir$block))

# Define the number of random effects
Nre<-length(unique(block))

# List all the required variables for JAGS
JAGS.data_full<-list(
  Y     =  trait_as_ir$total.w,
  Xc    =  Xc,
  Xb    =  Xb,
  Kc    =  Kc,
  Kb    =  Kb,
  N     = nrow(trait_as_ir),
  Block = block,
  Nre   = Nre,
  Zeros = rep(0, nrow(trait_as_ir)))

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

## Run the hurdle model
ir.semi.info <- jags( data        = JAGS.data_full,
                       inits       = inits,
                       parameters  = params,
                       model       = "ZA_LASSO.txt",
                       n.thin      = 10,
                       n.chains    = 3,
                       n.burnin    = 4000,
                       n.iter      = 5000)

# /!\ Can be long to run
ir.semi.info2 <- update(ir.semi.info, n.iter = 500000, n.thin = 10) # do 500,000 or 100,000,000 if possible

# Extract info
out.ir <-ir.semi.info2$BUGSoutput

# Check if the mixing is good ####
MyNames<- c(
  paste(c(colnames(Xc), "sigma1_Block") , "Gamma", sep = " "),
  paste(c(colnames(Xb), "sigma2_Block") , "Bern" , sep = " "),
  "r Gamma")

MyBUGSChains(out.ir, # run MCMCSupportHighstatV2.R
             c(uNames("beta",Kc), "sigma1_Block",
               uNames("gamma",Kb), "sigma2_Block",
               "r"),
             PanelNames = MyNames)

# get the numerical output of JAGS  ####
OUT.nir <- MyBUGSOutput(out.ir,  c(uNames("beta",Kc), "sigma1_Block",
                                    uNames("gamma",Kb),"sigma2_Block",
                                    "r"),
                        VarNames = MyNames)
print(OUT.nir, digits = 5)

# get the posterior distribution ####
# Gamma part of the model
MyNamesG<- c(colnames(Xc), "sigma1_Block","r")

MyBUGSHist(out.ir,
           c(uNames("beta",Kc), "sigma1_Block" ,"r"),
           PanelNames = MyNamesG)

# Bernoulli part of the model
MyNamesB<- c(colnames(Xb), "sigma2_Block")

MyBUGSHist(out.ir,
           c(uNames("gamma",Kb), "sigma2_Block"),
           PanelNames = MyNamesB)

# Get Pearson residuals and Expected values of the hurdle model ####
E <-out.ir$mean$PRes
ExpY <-out.ir$mean$ExpY

# Calculate de dispersion of the model ####
D<- sum(E^2)/(nrow(trait_as_ir) + (ncol(Xc)+ ncol(Xb)))

# Model validation of the Hurdle model ####
theme_set(theme_bw())

data1<-data.frame(E,ExpY)
Pres_ExpV<-ggplot(data= data1, aes( x= ExpY, y=E))+
  geom_point()+
  geom_hline(yintercept = 0, linetype="dotted")+
  ylab('Pearson residuals')+
  xlab('Expected values')

data2<-data.frame(trait_as_ir$treatment, E)
Pres_bio<-ggplot(data= data2, aes( x=trait_as_ir$treatment , y=E))+
  geom_boxplot()+
  geom_hline(yintercept = 0, linetype="dotted")+
  scale_x_discrete(labels=c("Temperate inoculant"="Temperate","Mixed inoculant"="Mixed","Boreal inoculant"="Boreal")) +
  theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1)) +
  ylab('Pearson residuals')+
  xlab('Inoculum source')

data3<-data.frame(trait_as_ir$forest, E)
Pres_abio<-ggplot(data= data3, aes( x= trait_as_ir$forest  , y=E))+
  geom_boxplot()+
  geom_hline(yintercept = 0, linetype="dotted")+
  theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1)) +
  ylab('Pearson residuals')+
  xlab('Forest type')

data4<-data.frame(trait_as_ir$total.w, ExpY)
ObvsExp<-ggplot(data= data4, aes(y = ExpY, x=trait_as_ir$total.w ))+
  geom_point()+
  labs(x ='Observed biomass (g)', y='Expected biomass (g)')

model_val<-ggarrange(Pres_ExpV,ObvsExp, Pres_abio ,Pres_bio,
                     ncol=2, nrow=2, labels=c("a", "b", "c", "d"))

# get the R2 of obs.value ~ expected.value
tt<-lm(ExpY~trait_as_ir$total.w)
summary(tt)$adj.r.squared

# sketch the model fit
MyData <- expand.grid(Treatment = levels(trait_as_ir$treatment), Forest = levels(trait_as_ir$forest))

X<-model.matrix(~ (Forest + Treatment)^2 , data= MyData)
beta.mcmc  <- out.ir$sims.list$beta
gamma.mcmc <- out.ir$sims.list$gamma

mu.mcmc    <- exp(X %*% t(beta.mcmc))
Pi.mcmc    <- exp(X %*% t(gamma.mcmc)) / (1 + exp(X %*% t(gamma.mcmc)))
ExpY.mcmc <- Pi.mcmc * mu.mcmc

MyData$key <- factor(1:12)

# Violin plots ####
soil_names <- c(Temperate = "Temperate soil origin", Mixed = "Mixed soil origin", Boreal = "Boreal soil origin")
theme_set(theme_bw())

# Hurdle posterior ####
# format data for ggplot and select 90% by forest
ExpY.mcmc.F <- ExpY.mcmc %>%
  t() %>%
  as.data.frame() %>%
  gather(factor_key = TRUE) %>%
  left_join(MyData, by = 'key') %>%
  group_by(Forest) %>%
  filter(value <= quantile(value, 0.95) & value >= quantile(value, 0.05))
  
# By forest
Hurdle_F<- ggplot(data = ExpY.mcmc.F, aes(x = Forest,  y = value)) +
  geom_violin()+
  stat_summary(fun.y=median, geom="point", shape=23, size=2)+
  theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1)) +
  ylim(0,1.6) +
  labs(x = "Soil origin", y = "Performance (g)")

## By treatment
ExpY.mcmc.T <- ExpY.mcmc %>%
  t() %>%
  as.data.frame() %>%
  gather(factor_key = TRUE) %>%
  left_join(MyData, by = 'key') %>%
  group_by(Treatment) %>%
  filter(value <= quantile(value, 0.95) & value >= quantile(value, 0.05))

Hurdle_T<- ggplot(data = ExpY.mcmc.T, aes(x = Treatment,  y = value)) +
  geom_violin()+
  stat_summary(fun.y=median, geom="point", shape=23, size=2)+
  theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1)) +
  scale_x_discrete(labels=c("Temperate inoculant"="Temperate","Mixed inoculant"="Mixed","Boreal inoculant"="Boreal")) +
  ylim(0,1.6) +
  labs(x = "Inoculum source", y = "")

## Overall
ExpY.mcmc.O <- ExpY.mcmc %>%
  t() %>%
  as.data.frame() %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  filter(value <= quantile(value, 0.95) & value >= quantile(value, 0.05)) %>%
  left_join(MyData, by = 'key')

Hurdle_O <- ggplot(data = ExpY.mcmc.O, aes(x = Treatment,  y = value)) +
  geom_violin()+
  facet_grid(~Forest, labeller = as_labeller(soil_names)) +
  stat_summary(fun.y=median, geom="point", shape=23, size=2)+
  theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1)) +
  scale_x_discrete(labels=c("Temperate inoculant"="Temperate","Mixed inoculant"="Mixed","Boreal inoculant"="Boreal")) +
  ylim(0,1.6) +
  labs(x = "Inoculum source", y = "Performance (g)")

# Bernoulli posterior ####
# format data for ggplot and select 90% by forest
Pi.mcmc.F <- Pi.mcmc %>%
  t() %>%
  as.data.frame() %>%
  gather(factor_key = TRUE) %>%
  left_join(MyData, by = 'key') %>%
  group_by(Forest) %>%
  filter(value <= quantile(value, 0.95) & value >= quantile(value, 0.05))

# By forest
Bern_F<- ggplot(data = Pi.mcmc.F, aes(x = Forest,  y = value)) +
  geom_violin()+
  stat_summary(fun.y=median, geom="point", shape=23, size=2)+
  theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1)) +
  ylim(0,1) +
  labs(x = "", y = "Survival probability")

# By treatment
Pi.mcmc.T <- Pi.mcmc %>%
  t() %>%
  as.data.frame() %>%
  gather(factor_key = TRUE) %>%
  left_join(MyData, by = 'key') %>%
  group_by(Treatment) %>%
  filter(value <= quantile(value, 0.95) & value >= quantile(value, 0.05))

Bern_T<- ggplot(data = Pi.mcmc.T, aes(x = Treatment,  y = value)) +
  geom_violin()+
  stat_summary(fun.y=median, geom="point", shape=23, size=2)+
  ylim(0,1) +
  theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1)) +
  scale_x_discrete(labels=c("Temperate inoculant"="Temperate","Mixed inoculant"="Mixed","Boreal inoculant"="Boreal")) +
  labs(x = "", y = "")

# Overall
Pi.mcmc.O <- Pi.mcmc %>%
  t() %>%
  as.data.frame() %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  filter(value <= quantile(value, 0.95) & value >= quantile(value, 0.05)) %>%
  left_join(MyData, by = 'key')

Bern_O<- ggplot(data = Pi.mcmc.O, aes(x = Treatment,  y = value)) +
  geom_violin()+
  facet_grid(~Forest, labeller = as_labeller(soil_names)) +
  stat_summary(fun.y=median, geom="point", shape=23, size=2)+
  ylim(0,1) +
  theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1)) +
  scale_x_discrete(labels=c("Temperate inoculant"="Temperate","Mixed inoculant"="Mixed","Boreal inoculant"="Boreal")) +
  labs(x = "", y = "Survival probability")

## Gamma posterior ####
## format data for ggplot and select CI 90% by forest
mu.mcmc.F <- mu.mcmc %>%
  t() %>%
  as.data.frame() %>%
  gather(factor_key = TRUE) %>%
  left_join(MyData, by = 'key') %>%
  group_by(Forest) %>%
  filter(value <= quantile(value, 0.95) & value >= quantile(value, 0.05))

## By forest
Gamm_F<- ggplot(data = mu.mcmc.F, aes(x = Forest,  y = value)) +
  geom_violin()+
  stat_summary(fun.y=median, geom="point", shape=23, size=2)+
  theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1)) +
  ylim(0,4) +
  labs(x = "", y = "Dry mass (of survivors, g)")

## By treatment
mu.mcmc.T <- mu.mcmc %>%
  t() %>%
  as.data.frame() %>%
  gather(factor_key = TRUE) %>%
  left_join(MyData, by = 'key') %>%
  group_by(Treatment) %>%
  filter(value <= quantile(value, 0.95) & value >= quantile(value, 0.05))

Gamm_T<- ggplot(data = mu.mcmc.T, aes(x = Treatment,  y = value)) +
  geom_violin()+
  stat_summary(fun.y=median, geom="point", shape=23, size=2)+
  theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1)) +
  scale_x_discrete(labels=c("Temperate inoculant"="Temperate","Mixed inoculant"="Mixed","Boreal inoculant"="Boreal")) +
  ylim(0,4) +
  labs(x = "", y = "Dry mass (g)")

## Overall
mu.mcmc.O <- mu.mcmc %>%
  t() %>%
  as.data.frame() %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  filter(value <= quantile(value, 0.95) & value >= quantile(value, 0.05)) %>%
  left_join(MyData, by = 'key')

Gamm_O<- ggplot(data = mu.mcmc.O, aes(x = Treatment,  y = value)) +
  geom_violin()+
  facet_grid(~Forest, labeller = as_labeller(soil_names)) +
  stat_summary(fun.y=median, geom="point", shape=23, size=2)+
  theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1)) +
  scale_x_discrete(labels=c("Temperate inoculant"="Temperate","Mixed inoculant"="Mixed","Boreal inoculant"="Boreal")) +
  ylim(0,4) +
  labs(x = "Inoculum source", y = "Dry mass of survivors (g)")

# Coefficient of Variation ####
# Abiotic
sd(c(expH[4,4], expH[8,4], expH[12,4]))/mean(c(expH[4,4], expH[8,4], expH[12,4]))*100

sd(c(
  mean(c(expH[1,4], expH[2,4], expH[3,4])),
  mean(c(expH[5,4], expH[6,4], expH[5,4])),
  mean(c(expH[9,4], expH[10,4], expH[11,4])))
)/mean(c(
  mean(c(expH[1,4], expH[2,4], expH[3,4])),
  mean(c(expH[5,4], expH[6,4], expH[5,4])),
  mean(c(expH[9,4], expH[10,4], expH[11,4])))
  )

# Biotic temperate
sd(c(expH[1,4], expH[2,4], expH[3,4]))/mean(c(expH[1,4], expH[2,4], expH[3,4]))*100

# Biotic Mixed
sd(c(expH[5,4], expH[6,4], expH[7,4]))/mean(c(expH[5,4], expH[6,4], expH[5,4]))*100

# Biotic boreal
sd(c(expH[9,4], expH[10,4], expH[11,4]))/mean(c(expH[9,4], expH[10,4], expH[11,4]))*100
