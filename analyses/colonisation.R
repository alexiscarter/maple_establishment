# Analysis of the mycorrhizal colonisation on sugar maple seedlings

# Load packages and data ####
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(ggpubr)

load("data/colo.rdata")

## Create length columns per categories ####
## create Length vector per observation (m)
len<-rep(0.00114,nrow(colo))

## fill new column with vector
colo$length<-len

## HYPHAE ####
## calculate length per idividual and per type (and leave out other numerical columns)
tot.length<-ddply(colo[,-c(8,9,10,11,13)], .(id, type), numcolwise(sum))
colnames(tot.length)[4]<-"total.length"

## calculate length per indiviual, type AND class for one structure
class.length<-ddply(colo[,-c(8,9,10,11,13)], .(id, type, hyp.amf), numcolwise(sum)) #add forest??
colnames(class.length)[4]<-"class.length"

## merge class and total length according to id and type
hyp_amf_data<-merge(class.length, tot.length, by=c("id", "type"))

#### create numerical class categories 
## create new class column
hyp_amf_data$class<-NA

## attribute numerical median to classes
hyp_amf_data[regexpr("0", hyp_amf_data[,3])>0,7]=0
hyp_amf_data[regexpr("1", hyp_amf_data[,3])>0,7]=0.125
hyp_amf_data[regexpr("2", hyp_amf_data[,3])>0,7]=0.375
hyp_amf_data[regexpr("3", hyp_amf_data[,3])>0,7]=0.650

hyp_amf_data$class<-as.numeric(hyp_amf_data$class)

#### calculate weight of class per length (no. of observations)
## relative length analysed corresponding to class values (Ltot.ci x ci)
weight<-function(x){
  x$weighted_lenght<-x$class.length*x$class
  return(x)
}

hyp_data_1<-weight(hyp_amf_data)

## clean data
NO_use<-c("hyp.amf.x", "hyp.amf.y","class")
hyp_data_1<-hyp_data_1[,!names(hyp_data_1) %in% NO_use]

## sum of all weighted data per type per id
sum_weighted_length<-ddply(hyp_data_1, .(id, type), numcolwise(sum))
colnames(sum_weighted_length)[5]<-"sum_weighted_length"

## add total.length to data.frame
hyp_data_2<-merge(sum_weighted_length, tot.length, by=c("id", "type") )

## clean data
hyp_data_2<-hyp_data_2[,-c(3,4,6) ]
colnames(hyp_data_2)[4]<- "total_length"
colnames(hyp_data_2)[3]<-"w_colo_length"

#### calculate the degree of colonization     
hyp_data_2$colonisation<-(hyp_data_2$w_colo_length/hyp_data_2$total_length)*100

## COIL ####
## calculate length per idividual and per type (and leave out other numerical columns)
tot.length<-ddply(colo[,-c(8,9,10,12,13)], .(id, type), numcolwise(sum))
colnames(tot.length)[4]<-"total.length"

## calculate length per indiviual, type AND class for one structure
class.length<-ddply(colo[,-c(8,9,10,12,13)], .(id, type, coil), numcolwise(sum)) #add forest??
colnames(class.length)[4]<-"class.length"

## merge class and total length according to id and type
coil_data<-merge(class.length, tot.length, by=c("id", "type"))

##### create numerical class categories
## create new class column
coil_data$class<-NA

## attribute numerical median to classes
coil_data[regexpr("0", coil_data[,3])>0,7]=0
coil_data[regexpr("1", coil_data[,3])>0,7]=0.125
coil_data[regexpr("2", coil_data[,3])>0,7]=0.375
coil_data[regexpr("3", coil_data[,3])>0,7]=0.650

coil_data$class<-as.numeric(coil_data$class)

#### calculate weight of class per length (no. of observations)
## relative length analysed corresponding to class values (Ltot.ci x ci)
weight<-function(x){
  x$weighted_lenght<-x$class.length*x$class
  return(x)
}
coil_1<-weight(coil_data)

## clean data
NO_use<-c("coil.x", "coil.y","class")
coil_1<-coil_1[,!names(coil_1) %in% NO_use]

## sum of all weighted data per type per id
sum_weighted_length_c<-ddply(coil_1, .(id, type), numcolwise(sum))
colnames(sum_weighted_length_c)[5]<-"sum_weighted_length"

## add total.length to data.frame
coil_2<-merge(sum_weighted_length_c, tot.length, by=c("id", "type"))

## clean data
coil_2<-coil_2[,-c(3,4,6) ]
colnames(coil_2)[4]<- "total_length"
colnames(coil_2)[3]<-"w_colo_length"

## calculate the degree of colonization     
coil_2$colonisation<-(coil_2$w_colo_length/coil_2$total_length)*100

## ARBUSCULES ####
## calculate length per idividual and per type (and leave out other numerical columns)
tot.length<-ddply(colo[,-c(8,10,11,12,13)], .(id, type), numcolwise(sum))
colnames(tot.length)[4]<-"total.length"

## calculate length per indiviual, type AND class for one structure
class.length<-ddply(colo[,-c(8,10,11,12,13)], .(id, type, arb), numcolwise(sum)) #add forest??
colnames(class.length)[4]<-"class.length"

## merge class and total length according to id and type
arb_data<-merge(class.length, tot.length, by=c("id", "type"))

#### create numerical class categories
## create new class column
arb_data$class<-NA

## attribute numerical median to classes
arb_data[regexpr("0", arb_data[,3])>0,7]=0
arb_data[regexpr("1", arb_data[,3])>0,7]=0.125
arb_data[regexpr("2", arb_data[,3])>0,7]=0.375
arb_data[regexpr("3", arb_data[,3])>0,7]=0.650

arb_data$class<-as.numeric(arb_data$class)

#### calculate weight of class per length (no. of observations)
## relative length analysed corresponding to class values (Ltot.ci x ci)
weight<-function(x){
  x$weighted_length<-x$class.length*x$class
  return(x)
}
arb_1<-weight(arb_data)

## clean data
NO_use<-c("arb.x", "arb.y","class")
arb_1<-arb_1[,!names(arb_1) %in% NO_use]

## sum of all weighted data per type per id
sum_weighted_length_a<-ddply(arb_1, .(id, type), numcolwise(sum))
colnames(sum_weighted_length_a)[5]<-"sum_weighted_length"

## add total.length to data.frame
arb_2<-merge(sum_weighted_length_a, tot.length, by=c("id", "type"))

## clean data
arb_2<-arb_2[,-c(3,4,6) ]

colnames(arb_2)[4]<- "total_length"
colnames(arb_2)[3]<-"w_colo_length"

## calculate the degree of colonization     
arb_2$colonisation<-(arb_2$w_colo_length/arb_2$total_length)*100


## VESICLES ####
## calculate length per idividual and per type (and leave out other numerical columns)
tot.length<-ddply(colo[,-c(8,9,11,12,13)], .(id, type), numcolwise(sum))
colnames(tot.length)[4]<-"total.length"

## calculate length per indiviual, type AND class for one structure
class.length<-ddply(colo[,-c(8,9,11,12,13)], .(id, type, ves), numcolwise(sum))
colnames(class.length)[4]<-"class.length"

## merge class and total length according to id and type
ves_data<-merge(class.length, tot.length, by=c("id", "type"))

#### create numerical class categories
## create new class column
ves_data$class<-NA

## attribute numerical median to classes
ves_data[regexpr("0", ves_data[,3])>0,7]=0
ves_data[regexpr("1", ves_data[,3])>0,7]=0.125
ves_data[regexpr("2", ves_data[,3])>0,7]=0.375
ves_data[regexpr("3", ves_data[,3])>0,7]=0.650

ves_data$class<-as.numeric(ves_data$class)

#### calculate weight of class per length (no. of observations)
## relative length analysed corresponding to class values (Ltot.ci x ci)
weight<-function(x){
  x$weighted_length<-x$class.length*x$class
  return(x)
}
ves_1<-weight(ves_data)

## clean data
NO_use<-c("ves.x", "ves.y","class")
ves_1<-ves_1[,!names(ves_1) %in% NO_use]

## sum of all weighted data per type per id
sum_weighted_length_a<-ddply(ves_1, .(id, type), numcolwise(sum))
colnames(sum_weighted_length_a)[5]<-"sum_weighted_length"

## add total.length to data.frame
ves_2<-merge(sum_weighted_length_a, tot.length, by=c("id", "type"))

## clean data
ves_2<-ves_2[,-c(3,4,6) ]

colnames(ves_2)[4]<- "total_length"
colnames(ves_2)[3]<-"w_colo_length"

## calculate the degree of colonization     
ves_2$colonisation<-(ves_2$w_colo_length/ves_2$total_length)*100

## ENDOPHYTES ####
## calculate length per idividual and per type (and leave out other numerical columns)
tot.length<-ddply(colo[,-c(8,9,10,11,12)], .(id, type), numcolwise(sum))
colnames(tot.length)[4]<-"total.length"

## calculate length per indiviual, type AND class for one structure
class.length<-ddply(colo[,-c(8,9,10,11,12)], .(id, type, endo), numcolwise(sum)) #add forest??
colnames(class.length)[4]<-"class.length"

## merge class and total length according to id and type
endo_data<-merge(class.length, tot.length, by=c("id", "type"))

#### create numerical class categories
## create new class column
endo_data$class<-NA

## attribute numerical median to classes
endo_data[regexpr("0", endo_data[,3])>0,7]=0
endo_data[regexpr("1", endo_data[,3])>0,7]=0.125
endo_data[regexpr("2", endo_data[,3])>0,7]=0.375
endo_data[regexpr("3", endo_data[,3])>0,7]=0.650

endo_data$class<-as.numeric(endo_data$class)

#### calculate weight of class per length (no. of observations)
## relative length analysed corresponding to class values (Ltot.ci x ci)
weight<-function(x){
  x$weighted_lenght<-x$class.length*x$class
  return(x)
}
endo_1<-weight(endo_data)

## clean data
NO_use<-c("endo.x", "endo.y","class")
endo_1<-endo_1[,!names(endo_1) %in% NO_use]

## sum of all weighted data per type per id
sum_weighted_length_e<-ddply(endo_1, .(id, type), numcolwise(sum))
colnames(sum_weighted_length_e)[5]<-"sum_weighted_length"

## add total.length to data.frame
endo_2<-merge(sum_weighted_length_e, tot.length, by=c("id", "type"))

## clean data
endo_2<-endo_2[,-c(3,4,6) ]

colnames(endo_2)[4]<- "total_length"
colnames(endo_2)[3]<-"w_colo_length"

## calculate the degree of colonization     
endo_2$colonisation<-(endo_2$w_colo_length/endo_2$total_length)*100

## Standardization of data taking under account mortality per forest ####
## separation infos retenues
colo2<-colo[,c(2,3,5,7)]
colo2<-colo2[!duplicated(colo2),]

## count number of survivors per forest per treatment (on 10 individuals)
colo3<-colo[,c(2,3,5,6)]
colo3<-colo3[!duplicated(colo3),]

Count_id<-colo3 %>% 
  dplyr::count(forest, treatment)

#### HYP.AMF
## merge retained info
hyp_data_3<-merge(hyp_data_2, colo2, by=c("id", "type")) 
## merge count with colonization data
hyp_data_4<-merge(hyp_data_3, Count_id, by=c("forest", "treatment")) 
## obtain ratio of colonization normalized by number of survivors
hyp_data_4$norm_colo<-(hyp_data_4$colonisation*((hyp_data_4$n)/10))

#### Coil 
## merge retained info
coil_3<-merge(coil_2, colo2, by=c("id", "type")) 
## merge count with colonization data
coil_4<-merge(coil_3, Count_id, by=c("forest", "treatment")) 
## obtain ratio of colonization normalized by number of survivors
coil_4$norm_colo<-(coil_4$colonisation*((coil_4$n)/10))

#### ARB 
## merge retained info
arb_3<-merge(arb_2, colo2, by=c("id", "type")) 
## merge count with colonization data
arb_4<-merge(arb_3, Count_id, by=c("forest", "treatment")) 
## obtain ratio of colonization normalized by number of survivors
arb_4$norm_colo<-(arb_4$colonisation*((arb_4$n)/10))

#### VES
## merge retained info
ves_3<-merge(ves_2, colo2, by=c("id", "type")) 
## merge count with colonization data
ves_4<-merge(ves_3, Count_id, by=c("forest", "treatment")) 
## obtain ratio of colonization normalized by number of survivors
ves_4$norm_colo<-(ves_4$colonisation*((ves_4$n)/10))

#### ENDO
## merge retained info
endo_3<-merge(endo_2, colo2, by=c("id", "type")) 
## merge count with colonization data
endo_4<-merge(endo_3, Count_id, by=c("forest", "treatment")) 
## obtain ratio of colonization normalized by number of survivors
endo_4$norm_colo<-(endo_4$colonisation*((endo_4$n)/10))


# Figure with treated and untreated samples
colo_hyp_irn<-ggplot(hyp_data_4, mapping=aes(x= treatment, y=colonisation)) +
  facet_grid(~forest, labeller = as_labeller(soil_names)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  labs(x = "Inoculum source", y = "Hyphal colonization intensity (%)") +
  scale_x_discrete(labels=c("Temperate inoculant"="Temperate","Mixed inoculant"="Mixed","Boreal inoculant"="Boreal")) +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

colo_hyp_irnn<-ggplot(hyp_data_4, mapping=aes(x= treatment, y=colonisation)) +
  facet_grid(~forest, labeller = as_labeller(soil_names)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  labs(x = "Inoculum source", y = "Hyphae") +
  scale_x_discrete(labels=c("Temperate inoculant"="Temperate","Mixed inoculant"="Mixed","Boreal inoculant"="Boreal")) +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

colo_coil_irn<-ggplot(coil_4, mapping=aes(x= treatment, y=colonisation, cex.main=0.5 )) +
  facet_grid(~forest, labeller = as_labeller(soil_names)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  labs(x = "Inoculum source", y = "Coil") +
  scale_x_discrete(labels=c("Temperate inoculant"="Temperate","Mixed inoculant"="Mixed","Boreal inoculant"="Boreal")) +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

colo_arb_irn<-ggplot(arb_4, mapping=aes(x= treatment, y=colonisation, cex.main=0.5 )) +
  facet_grid(~forest, labeller = as_labeller(soil_names)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  labs(x = "Inoculum source", y = "Arbuscules") +
  scale_x_discrete(labels=c("Temperate inoculant"="Temperate","Mixed inoculant"="Mixed","Boreal inoculant"="Boreal")) +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

colo_ves_irn<-ggplot(ves_4, mapping=aes(x= treatment, y=colonisation, cex.main=0.5 )) +
  facet_grid(~forest, labeller = as_labeller(soil_names)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  labs(x = "Inoculum source", y = "Vesicles") +
  scale_x_discrete(labels=c("Temperate inoculant"="Temperate","Mixed inoculant"="Mixed","Boreal inoculant"="Boreal")) +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

colo_endo_irn<-ggplot(endo_4, mapping=aes(x= treatment, y=colonisation, cex.main=0.5 )) +
  facet_grid(~forest, labeller = as_labeller(soil_names)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  labs(x = "Inoculum source", y = "Endophytes") +
  scale_x_discrete(labels=c("Temperate inoculant"="Temperate","Mixed inoculant"="Mixed","Boreal inoculant"="Boreal")) +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

all_irn <- ggarrange(colo_arb_irn, colo_coil_irn, colo_hyp_irnn, colo_ves_irn, colo_endo_irn,
                     ncol=1, nrow=5, labels=c("a", "b", "c", "d", "e"))  

## Keep only untreated samples (experiment 1)
hyp_nir <- hyp_data_4 %>% subset(treatment == "Untreated")
coil_nir <- coil_4 %>% subset(treatment == "Untreated")
arb_nir <- arb_4 %>% subset(treatment == "Untreated")
ves_nir <- ves_4 %>% subset(treatment == "Untreated")
endo_nir <- endo_4 %>% subset(treatment == "Untreated")

## Remove untreated samples (experiment 2)
hyp_data_4 <- hyp_data_4 %>% subset(!treatment == "Untreated")
coil_4 <- coil_4 %>% subset(!treatment == "Untreated")
arb_4 <- arb_4 %>% subset(!treatment == "Untreated")
ves_4 <- ves_4 %>% subset(!treatment == "Untreated")
endo_4 <- endo_4 %>% subset(!treatment == "Untreated")


## Figures on observed data (Experiment 2) #### 
soil_names <- c(Temperate = "Temperate soil origin", Mixed = "Mixed soil origin", Boreal = "Boreal soil origin")
#### Hyphae
colo_hyp<-ggplot(hyp_data_4, mapping=aes(x= treatment, y=colonisation)) +
  facet_grid(~forest, labeller = as_labeller(soil_names)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  labs(x = "Inoculum source", y = "Hyphal colonization intensity (%)") +
  scale_x_discrete(labels=c("Temperate inoculant"="Temperate","Mixed inoculant"="Mixed","Boreal inoculant"="Boreal")) +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

colo_hyp_treatment<-ggplot(hyp_data_4, mapping=aes(x= treatment, y=colonisation)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  labs(x = "Inoculum source", y = "") +
  scale_x_discrete(labels=c("Temperate inoculant"="Temperate","Mixed inoculant"="Mixed","Boreal inoculant"="Boreal")) +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

colo_hyp_forest<-ggplot(hyp_data_4, mapping=aes(x= forest, y=colonisation)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  labs(x = "Soil origin", y = "Hyphal colonization\nintensity (%)") +
  scale_x_discrete(labels=c("Temperate forest"="Temperate","Mixed forest"="Mixed","Boreal forest"="Boreal")) +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

#### CoiL
colo_coil<-ggplot(coil_4, mapping=aes(x= treatment, y=colonisation, cex.main=0.5 )) +
  facet_grid(~forest, labeller = as_labeller(soil_names)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  labs(x = "", y = "") +
  scale_x_discrete(labels=c("Temperate inoculant"="Temperate","Mixed inoculant"="Mixed","Boreal inoculant"="Boreal")) +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

colo_coil_treatment<-ggplot(coil_4, mapping=aes(x= treatment, y=colonisation, cex.main=0.5 )) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  labs(x = "", y = "") +
  scale_x_discrete(labels=c("Temperate inoculant"="Temperate","Mixed inoculant"="Mixed","Boreal inoculant"="Boreal")) +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

colo_coil_forest<-ggplot(coil_4, mapping=aes(x= forest, y=colonisation)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  labs(x = "", y = "Colonization intensity\nby coils (%)") +
  scale_x_discrete(labels=c("Temperate forest"="Temperate","Mixed forest"="Mixed","Boreal forest"="Boreal")) +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

#### Arbuscules
colo_arb<-ggplot(arb_4, mapping=aes(x= treatment, y=colonisation, cex.main=0.5 )) +
  facet_grid(~forest, labeller = as_labeller(soil_names)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  labs(x = "", y = "") +
  scale_x_discrete(labels=c("Temperate inoculant"="Temperate","Mixed inoculant"="Mixed","Boreal inoculant"="Boreal")) +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

colo_arb_treatment<-ggplot(arb_4, mapping=aes(x= treatment, y=colonisation, cex.main=0.5 )) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  labs(x = "", y = "") +
  scale_x_discrete(labels=c("Temperate inoculant"="Temperate","Mixed inoculant"="Mixed","Boreal inoculant"="Boreal")) +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

colo_arb_forest<-ggplot(arb_4, mapping=aes(x= forest, y=colonisation)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  labs(x = "", y = "Colonization intensity\nby arbuscules (%)") +
  scale_x_discrete(labels=c("Temperate forest"="Temperate","Mixed forest"="Mixed","Boreal forest"="Boreal")) +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

#### Vesicles
colo_ves<-ggplot(ves_4, mapping=aes(x= treatment, y=colonisation, cex.main=0.5 )) +
  facet_grid(~forest, labeller = as_labeller(soil_names)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  labs(x = "", y = "") +
  scale_x_discrete(labels=c("Temperate inoculant"="Temperate","Mixed inoculant"="Mixed","Boreal inoculant"="Boreal")) +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

colo_ves_treatment<-ggplot(ves_4, mapping=aes(x= treatment, y=colonisation, cex.main=0.5 )) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  labs(x = "", y = "") +
  scale_x_discrete(labels=c("Temperate inoculant"="Temperate","Mixed inoculant"="Mixed","Boreal inoculant"="Boreal")) +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

colo_ves_forest<-ggplot(ves_4, mapping=aes(x= forest, y=colonisation)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  labs(x = "", y = "Colonization intensity\nby vesicles (%)") +
  scale_x_discrete(labels=c("Temperate forest"="Temperate","Mixed forest"="Mixed","Boreal forest"="Boreal")) +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

#### Endophytes
colo_endo<-ggplot(endo_4, mapping=aes(x= treatment, y=colonisation, cex.main=0.5 )) +
  facet_grid(~forest, labeller = as_labeller(soil_names)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  labs(x = "Inoculum source", y = ")") +
  scale_x_discrete(labels=c("Temperate inoculant"="Temperate","Mixed inoculant"="Mixed","Boreal inoculant"="Boreal")) +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

colo_endo_treatment<-ggplot(endo_4, mapping=aes(x= treatment, y=colonisation, cex.main=0.5 )) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  labs(x = "Inoculum source", y = "") +
  scale_x_discrete(labels=c("Temperate inoculant"="Temperate","Mixed inoculant"="Mixed","Boreal inoculant"="Boreal")) +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

colo_endo_forest<-ggplot(endo_4, mapping=aes(x= forest, y=colonisation)) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  labs(x = "Soil origin", y = "Colonization intensity\nby endophytes (%)") +
  scale_x_discrete(labels=c("Temperate forest"="Temperate","Mixed forest"="Mixed","Boreal forest"="Boreal")) +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

## Figures, only untreated data (experiment 1) ####
colo_arb_nir <-ggplot(arb_nir, mapping=aes(x= forest, y=colonisation, cex.main=0.5 )) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  labs(x = "", y = "Colonization intensity\nby arbuscules (%)") +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

colo_hyp_nir <- ggplot(hyp_nir, mapping=aes(x= forest, y=colonisation, cex.main=0.5 )) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  labs(x = "", y = "Hyphal colonization\nintensity (%)") +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

colo_coil_nir <- ggplot(coil_nir, mapping=aes(x= forest, y=colonisation, cex.main=0.5 )) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  labs(x = "", y = "Colonization intensity\nby coils (%)") +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

colo_ves_nir <- ggplot(ves_nir, mapping=aes(x= forest, y=colonisation, cex.main=0.5 )) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  labs(x = "", y = "Colonization intensity\nby vesicles (%)") +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))

colo_endo_nir <- ggplot(endo_nir, mapping=aes(x= forest, y=colonisation, cex.main=0.5 )) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position="dodge", width = .3, fun.args=list(conf.int=.68)) +
  stat_summary(fun.data = mean_cl_boot, geom = "bar", position="dodge", width = .8, fun.args=list(conf.int=0)) +
  labs(x = "Soil", y = "Colonization intensity\nby endophytes (%)") +
  theme_bw() + theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1), axis.title = element_text(size=10), strip.text.x = element_text(size=10))


## figures for the article
colo_hyp

colo_all <- ggarrange(colo_arb_forest, colo_arb_treatment, colo_arb,
                      colo_coil_forest, colo_coil_treatment, colo_coil,
                      colo_ves_forest, colo_ves_treatment, colo_ves,
                      colo_endo_forest, colo_endo_treatment, colo_endo,
                      ncol=3, nrow=4, widths =c(1,1,2),
                      labels=c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"))

colo_nir_all <- ggarrange(colo_arb_nir, colo_hyp_nir, colo_coil_nir, colo_ves_nir, colo_endo_nir,
                      ncol=1, nrow=5,
                      labels=c("a", "b", "c", "d", "e"))
