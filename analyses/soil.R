#### Analysis of soil data from Mégantic

## Load packages and data ####
library(ggplot2)
library(nlme)
library(emmeans)
library(ggpubr)

load('data/soil.rdata')

# plot settings
soil_names <- c(FT = "Temperate forest soil", FM = "Mixed forest soil", FB = "Boreal forest soil")
theme_set(theme_bw())

## Keep only H, Ae, and B for cations analysis (L and F not analyzed)
soil.cat <- soil[soil$horizon %in% c('H', 'Ae', 'B'), ]

## Full models ####

## pH ####
mod.ph1 <- lme(pH.CaCl2 ~ forest*horizon, random = ~ 1|block,
              method = 'REML',
              weights = NULL, data = soil)
anova(mod.ph1)
plot(mod.ph1)
qqnorm(resid(mod.ph1)); qqline(resid(mod.ph1))

## forest*horizon
ph.forest.horiz <- emmeans(mod.ph1, pairwise ~ forest*horizon, adjust = "tukey")
mult.ph.forest.horiz <- CLD(ph.forest.horiz, alpha=0.05, Letters=letters, adjust="tukey", reversed = TRUE)
plot.mult.ph.forest.horiz <- ggplot(mult.ph.forest.horiz, aes(x = horizon, y = emmean)) + 
  geom_point(stat = "identity", position=position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.5), width = 0.2) +
  facet_grid(~forest, labeller = as_labeller(soil_names))+
  labs(x="Horizons", y= expression(paste("pH (in CaCl"[2],")")))

## forest
ph.forest <- emmeans(mod.ph1, pairwise ~ forest, adjust = "tukey")
mult.ph.forest <- CLD(ph.forest, alpha=0.05, Letters=letters, adjust="tukey", reversed = TRUE)
plot.mult.ph.forest <- ggplot(mult.ph.forest, aes(x = forest, y = emmean)) + 
  geom_point(stat = "identity", position=position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.5), width = 0.2) +
  scale_x_discrete(labels=c("Temperate", "Mixed", "Boreal")) +
  theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1)) +
  labs(x="", y= expression(paste("pH (in CaCl"[2],")")))

## ECEC ####
mod.ecec1 <- lme(ECEC ~ forest*horizon, random = ~ 1|block,
               method = 'REML',
               weights = NULL, data = soil.cat)
anova(mod.ecec1)
plot(mod.ecec1)
qqnorm(resid(mod.ecec1)); qqline(resid(mod.ecec1))

## forest*horizon
ecec.forest.horiz <- emmeans(mod.ecec1, pairwise ~ forest*horizon, adjust = "tukey")
mult.ecec.forest.horiz <- CLD(ecec.forest.horiz, alpha=0.05, Letters=letters, adjust="tukey", reversed = TRUE)
plot.mult.ecec.forest.horiz <- ggplot(mult.ecec.forest.horiz, aes(x = horizon, y = emmean)) + 
  geom_point(stat = "identity", position=position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.5), width = 0.2) +
  facet_grid(~forest, labeller = as_labeller(soil_names))+
  labs(x="Horizons", y= expression(paste("ECEC (cmol"[c]," kg"^"-1",")")))

## forest
ecec.forest <- emmeans(mod.ecec1, pairwise ~ forest, adjust = "tukey")
mult.ecec.forest <- CLD(ecec.forest, alpha=0.05, Letters=letters, adjust="tukey", reversed = TRUE)
plot.mult.ecec.forest <- ggplot(mult.ecec.forest, aes(x = forest, y = emmean)) + 
  geom_point(stat = "identity", position=position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.5), width = 0.2) +
  scale_x_discrete(labels=c("Temperate", "Mixed", "Boreal")) +
  theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1)) +
  labs(x="Forest type", y= expression(paste("ECEC (cmol"[c]," kg"^"-1",")")))

## C:N ratio ####
mod.CN3 <- lme(totalC/totalN ~ forest*horizon, random = ~ 1|block,
               method = 'REML',
               weights = varExp(), data = soil)
anova(mod.CN3)
plot(mod.CN3)
qqnorm(resid(mod.CN3)); qqline(resid(mod.CN3))

## forest*horizon
cn.forest.horiz <- emmeans(mod.CN3, pairwise ~ forest*horizon, adjust = "tukey")
mult.cn.forest.horiz <- CLD(cn.forest.horiz, alpha=0.05, Letters=letters, adjust="tukey", reversed = TRUE)
plot.mult.cn.forest.horiz <- ggplot(mult.cn.forest.horiz, aes(x = horizon, y = emmean)) + 
  geom_point(stat = "identity", position=position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.5), width = 0.2) +
  facet_grid(~forest, labeller = as_labeller(soil_names))+
  labs(x="Horizons", y= "C:N ratio")

## forest 
cn.forest <- emmeans(mod.CN3, pairwise ~ forest, adjust = "tukey")
mult.cn.forest <- CLD(cn.forest, alpha=0.05, Letters=letters, adjust="tukey", reversed = TRUE)
plot.mult.cn.forest <- ggplot(mult.cn.forest, aes(x = forest, y = emmean)) + 
  geom_point(stat = "identity", position=position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.5), width = 0.2) +
  scale_x_discrete(labels=c("Temperate", "Mixed", "Boreal")) +
  theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1)) +
  labs(x="", y= "C:N ratio")

## TOTAL P ####
mod.tp2 <- lme(totalP ~ forest*horizon, random = ~ 1|block,
               method = 'REML',
               weights = varPower(), data = soil)
anova(mod.tp2)
plot(mod.tp2)
qqnorm(resid(mod.tp2)); qqline(resid(mod.tp2))

## forest*horizon
tp.forest.horiz <- emmeans(mod.tp2, pairwise ~ forest*horizon, adjust = "tukey")
mult.tp.forest.horiz <- CLD(tp.forest.horiz, alpha=0.05, Letters=letters, adjust="tukey", reversed = TRUE)
plot.mult.tp.forest.horiz <- ggplot(mult.tp.forest.horiz, aes(x = horizon, y = emmean)) + 
  geom_point(stat = "identity", position=position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.5), width = 0.2) +
  facet_grid(~forest, labeller = as_labeller(soil_names))+
  labs(x="Horizons", y= expression(paste("Total P (mg kg"^"-1",")")))

## forest
tp.forest <- emmeans(mod.tp2, pairwise ~ forest, adjust = "tukey")
mult.tp.forest <- CLD(tp.forest, alpha=0.05, Letters=letters, adjust="tukey", reversed = TRUE)
plot.mult.tp.forest <- ggplot(mult.tp.forest, aes(x = forest, y = emmean)) + 
  geom_point(stat = "identity", position=position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.5), width = 0.2) +
  scale_x_discrete(labels=c("Temperate", "Mixed", "Boreal")) +
  theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1)) +
  labs(x="", y= expression(paste("Total P (mg kg"^"-1",")")))

## LABILE P ####
mod.bp2 <- lme(BrayP ~ forest*horizon, random = ~ 1|block,
               method = 'REML',
               weights = varPower(), data = soil)
anova(mod.bp2)
plot(mod.bp2)
qqnorm(resid(mod.bp2)); qqline(resid(mod.bp2))

## forest*horizon
bp.forest.horiz <- emmeans(mod.bp2, pairwise ~ forest*horizon, adjust = "tukey")
mult.bp.forest.horiz <- CLD(bp.forest.horiz, alpha=0.05, Letters=letters, adjust="tukey", reversed = TRUE)
plot.bp.mult.forest.horiz <- ggplot(mult.bp.forest.horiz, aes(x = horizon, y = emmean)) + 
  geom_point(stat = "identity", position=position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.5), width = 0.2) +
  facet_grid(~forest, labeller = as_labeller(soil_names))+
  labs(x="Horizons", y= expression(paste("Labile P (mg kg"^"-1",")")))

## forest
bp.forest <- emmeans(mod.bp2, pairwise ~ forest, adjust = "tukey")
mult.bp.forest <- CLD(bp.forest, alpha=0.05, Letters=letters, adjust="tukey", reversed = TRUE)
plot.mult.bp.forest <- ggplot(mult.bp.forest, aes(x = forest, y = emmean)) + 
  geom_point(stat = "identity", position=position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.5), width = 0.2) +
  scale_x_discrete(labels=c("Temperate", "Mixed", "Boreal")) +
  theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1)) +
  labs(x="", y= expression(paste("Labile P (mg kg"^"-1",")")))

## BASE SATURATION ####
mod.bs2 <- lme(BS ~ forest*horizon, random = ~ 1|block,
               method = 'REML',
               weights = varPower(), data = soil.cat)
anova(mod.bs2)
plot(mod.bs2)
qqnorm(resid(mod.bs2)); qqline(resid(mod.bs2))

## forest*horizon
bs.forest.horiz <- emmeans(mod.bs2, pairwise ~ forest*horizon, adjust = "tukey")
mult.bs.forest.horiz <- CLD(bs.forest.horiz, alpha=0.05, Letters=letters, adjust="tukey", reversed = TRUE)
plot.mult.bs.forest.horiz <- ggplot(mult.bs.forest.horiz, aes(x = horizon, y = emmean)) + 
  geom_point(stat = "identity", position=position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.5), width = 0.2) +
  facet_grid(~forest, labeller = as_labeller(soil_names))+
  labs(x="Horizons", y= 'Base Sat. (%)')

## forest
bs.forest <- emmeans(mod.bs2, pairwise ~ forest, adjust = "tukey")
mult.bs.forest <- CLD(bs.forest, alpha=0.05, Letters=letters, adjust="tukey", reversed = TRUE)
plot.mult.bs.forest <- ggplot(mult.bs.forest, aes(x = forest, y = emmean)) + 
  geom_point(stat = "identity", position=position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), position=position_dodge(width=0.5), width = 0.2) +
  scale_x_discrete(labels=c("Temperate", "Mixed", "Boreal")) +
  theme(axis.text.x = element_text(angle=45, vjust= 1, size = 10, hjust=1)) +
  labs(x="Forest type", y= 'Base Sat. (%)')

## Arrange FOREST in one page
## Note: Tukey post hoc analysis not included in any plots

soil_forest<-ggarrange(plot.mult.ph.forest, plot.mult.cn.forest, 
                       plot.mult.tp.forest, plot.mult.bp.forest, 
                       plot.mult.ecec.forest, plot.mult.bs.forest,
                       ncol=2, nrow=3, labels=c("a", "b", "c", "d", "e", "f"))

## Arrange FOREST*HORIZONS in one page
soil_forest_horiz<-ggarrange(plot.mult.ph.forest.horiz, plot.mult.cn.forest.horiz, 
                       plot.mult.tp.forest.horiz, plot.bp.mult.forest.horiz, 
                       plot.mult.ecec.forest.horiz, plot.mult.bs.forest.horiz,
                       ncol=2, nrow=3, labels=c("a", "b", "c", "d", "e", "f"))


### DEPTH by horizons ####
library(tidyr)
library(dplyr)
library(plyr)

## load data
soil_thickness <- read.csv("data/soil_thick.csv", sep = ";")

## manipulate
thick <- gather(soil_thickness, Horizon, Thickness, -Forest, -Block, factor_key=TRUE)
thick$Forest <- recode(thick$Forest, FT="Temperate", FM="Mixed",FB= "Boreal")
thick$Forest <- factor(thick$Forest, c("Temperate", "Mixed", "Boreal"))
thick$Horizon <- factor(thick$Horizon, c("B", "Ae", "H", "F", "L"))

## Mean by forest and block
thick_mean <- ddply(thick ,.(Forest, Block, Horizon), numcolwise(mean))

ggplot(thick_mean, aes(x = Forest, y = Thickness, fill = Horizon)) +
  geom_bar(stat = "summary", fun.y = "mean")+
  scale_y_reverse() +
  guides(fill = guide_legend(reverse=T)) +
  labs(y = "Soil depth (cm)") +
  scale_fill_manual(values = c('darkorange2', 'grey60', 'grey1', 'sienna4', 'darkgreen')) +
  theme(panel.background = element_blank(), panel.border = element_blank(),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "black"), 
  panel.grid.minor = element_line(size = 0.5, linetype = 'solid',colour = "black"),
  axis.line = element_blank())
