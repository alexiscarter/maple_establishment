## Analyses with root hyphal colonization as an explanatory variable of dry mass, 
## and with soil characteristics as explanatory variables of performance
## using generalized linear mixed-effect models.

library(tidyverse); theme_set(theme_bw())
library(brms)
library(brmstools)
library(tidybayes)
library(ggpubr)

colo_as <- read.csv('data/colo_as.csv')
traits_as <- read.csv('data/traits_as.csv')
soil <- read.csv('data/soil.csv')

## Colonization variable ####
colo_w_mod <- colo_as %>% 
  mutate(colo_tot = rowSums(colo_as[, c('mean_arb','mean_hyp','mean_ves','mean_coil')])) %>%
  left_join(traits_as, by = 'id') %>%
  subset(!treatment == "NIR") %>% 
  subset(total.w > 0) %>% 
  subset(!colo_tot > 80) %>% # remove extreme values, id 124 32
  select("id", "total.w", "mean_arb", "mean_hyp", "mean_ves", "mean_endo", "mean_coil","colo_tot", "forest", "block", "treatment")

## Model
mod.colo.gamm <- brm(formula = total.w ~ mean_hyp + (1|block),
                data = colo_w_mod, family = Gamma(link="log"),
                warmup = 1000, iter = 10000, chains = 4, thin = 10,
                control = list(adapt_delta = 0.90))
mod.colo.gamm.1 <- update(mod.colo.gamm, iter = 50000) 
summary(mod.colo.gamm.1, prob = 0.90, mc_se = TRUE, priors = TRUE)

## Plot marginal effect
plot.colo <- plot(marginal_effects(mod.colo.gamm.1),points=T)
plot.colo$mean_hyp +
  labs(x= 'Root hyphal colonization (%)', y = "Dry mass (of survivors, g)")

## Soil chemistry variables ##### 
soil.mean <- soil %>%
  select("forest", "block", "pH.CaCl2", "totalC", "totalN", "totalP", "BrayP", "ECEC", "BS") %>%
  group_by(forest, block) %>%
  summarise_at(vars(pH.CaCl2:BS), mean, na.rm = TRUE) %>% 
  mutate(CN = totalC/totalN)

soil_w_ir_hgamm <- traits_as %>% 
  subset(!treatment == "NIR") %>%
  select(id, forest, block, total.w) %>% 
  left_join(soil.mean, by = c("forest", "block")) %>% 
  mutate(CN = totalC/totalN) %>% 
  rename(pH = pH.CaCl2, biomass = total.w)

mod.ir.hgamm.soil <- brm(formula = biomass ~ pH + CN + totalP + BrayP + ECEC + BS + (1|block),
                           data = soil_w_ir_hgamm, family = hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"),
                           warmup = 1000, iter = 10000, chains = 4, thin = 10,
                           control = list(adapt_delta = 0.95))
mod.ir.hgamm.soil.1 <- update(mod.ir.hgamm.soil, iter = 50000) 
summary(mod.ir.hgamm.soil.1, prob = 0.90, mc_se = TRUE, priors = TRUE)

## Plot marginal effects
fit.hgamm.soil.gg <- plot(marginal_effects(mod.ir.hgamm.soil.1),points=T, prob = 0.90, plot = FALSE)

plot.tp <- fit.hgamm.soil.gg$totalP +
  labs(y = "Performance (g)", x = expression(paste("Total P (mg kg"^"-1",")")))

plot.bp <- fit.hgamm.soil.gg$BrayP +
  labs(y = "Performance (g)", x = expression(paste("Labile P (mg kg"^"-1",")")))

plot.ecec <- fit.hgamm.soil.gg$ECEC +
  labs(y ="Performance (g)", x = expression(paste("ECEC (cmol"[c]," kg"^"-1",")")))

plot.bs <- fit.hgamm.soil.gg$BS +
  labs(y ="Performance (g)", x = 'Base Sat. (%)')

plot.ph <- fit.hgamm.soil.gg$pH +
  labs(y ="Performance (g)", x = expression(paste("pH (in CaCl"[2],")")))

plot.cn <- fit.hgamm.soil.gg$CN +
  labs(y ="Performance (g)", x = 'C:N ratio')

fit.hgamm.soil.coef$data$Parameter <- c('Labile P', 'Base Sat.', 'C:N', 'ECEC', 'pH', 'Total P')
plot.coef<- fit.hgamm.soil.coef +
  labs(y = "", x = 'Estimate')
