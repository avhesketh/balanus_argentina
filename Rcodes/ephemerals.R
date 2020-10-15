# May 2020 after Chris meeting to finalize stats

## clean time-adjusted model for algal cover
## after much discussion, should use 3 dates for each response
## and also use a generalized linear model with a Tweedie distribution

library(tidyverse)
library(glmmTMB)
library(car)
library(DHARMa)
library(ggcorrplot)

# starting with ephemeral algae

ephemerals <- read_csv("../data/ephemerals.csv") %>% 
  select(-1)
glimpse(ephemerals)

# now visualize the response variable

hist(ephemerals$percent_cover, breaks = 100)

# overdispersed with most of the weight at or close to zero
# suited to a Tweedie distribution without transformation
# first, restrict to desired time points

eph_arg <- ephemerals %>% 
  filter(location == "GN")

eph_bc <- ephemerals %>% 
  filter(location == "BS")

# argentina data are missing initial time points since ephemerals counted then
# are not real data ... hiding under barnacles, so trmt specific
levels(as.factor(eph_arg$date))
eph_arg_keep <- eph_arg %>% 
  filter(date == "2006-04-26" | date == "2006-10-07" | date == "2007-02-19")
levels(as.factor(eph_bc$date))
eph_bc_keep <- eph_bc %>%
  filter(date == "2006-10-06" | date == "2007-04-06" | date == "2007-07-17")

ephemerals_keep <- eph_bc_keep %>% 
  full_join(eph_arg_keep)

# have kept the point closest to fall eq., then spring eq., then end of experiment

plot(ephemerals_keep$percent_cover ~ ephemerals_keep$timediff, col = c("blue","red")[as.factor(ephemerals_keep$location)])
plot(ephemerals_keep$percent_cover ~ ephemerals_keep$timediff, col = c("blue","red", "green")[as.factor(ephemerals_keep$limpets)])
plot(ephemerals_keep$percent_cover ~ ephemerals_keep$timediff, col = c("blue","red")[as.factor(ephemerals_keep$barnacles)])

# more variation in argentina plots than in bc, where limpets are excluded. Barnacle effect harder to see.

eph_exclude <- ephemerals_keep %>% 
  filter(limpets == "out")

# try keeping all the limpet treatments

glmm.eph.all.limp <- glmmTMB(percent_cover ~ (timediff + location + barnacles + limpets)^2 
                             + (1|block/location),
                             dispformula = ~limpets + location,
                      data = ephemerals_keep,
                      family = tweedie())

plot(simulateResiduals(glmm.eph.all.limp))
AIC(glmm.eph.all.limp)
summary(glmm.eph.all.limp)

res.glmm.eph.all <- residuals(glmm.eph.all.limp)
plot(res.glmm.eph.all ~ as.factor(ephemerals_keep$location))
plot(res.glmm.eph.all ~ as.factor(ephemerals_keep$barnacles))
plot(res.glmm.eph.all ~ as.factor(ephemerals_keep$limpets))
plot(res.glmm.eph.all ~ ephemerals_keep$timediff)

drop1(glmm.eph.all.limp)

glmm.eph.all.a <- update(glmm.eph.all.limp, ~. - barnacles:limpets)

drop1(glmm.eph.all.a, test = "Chisq")

glmm.eph.all.b <- update(glmm.eph.all.a, ~. - timediff:location)

drop1(glmm.eph.all.b, test = "Chisq")

## stop here using p values

summary(glmm.eph.all.b)
Anova(glmm.eph.all.b)

acf(residuals(glmm.eph.all.b))
plot(simulateResiduals(glmm.eph.all.a))

# can't do higher order interactions with limpet treatments

glmm.eph.1 <- glmmTMB(percent_cover ~ (timediff + location + barnacles)^3 +
                        (1|block/location),
                      data = eph_exclude,
                      family = tweedie())

resdh.glm.eph.1 <- simulateResiduals(glmm.eph.1)
plot(resdh.glm.eph.1)

# QQ plot and residuals are actually pretty good, but not too perfect that it looks bad

res.glm.eph.1 <- residuals(glmm.eph.1)
acf(res.glm.eph.1)

# autocorrelation is under control, so no need to add that element

plot(res.glm.eph.1 ~ as.factor(eph_exclude$location))
plot(res.glm.eph.1 ~ eph_exclude$timediff)
plot(res.glm.eph.1 ~ as.factor(eph_exclude$barnacles))

# do we need block term?

glmm.eph.2 <- glmmTMB(percent_cover ~ (timediff + location + barnacles)^3,
                      data = eph_exclude,
                      family = tweedie())

AIC(glmm.eph.1, glmm.eph.2)
summary(glmm.eph.2)
# AIC definitely lower without block term

# variance between locations is highest for any variable
# include that too in dispformula? what does AIC say

glmm.eph.3 <- glmmTMB(percent_cover ~ (timediff + location + barnacles)^3,
                      data = eph_exclude,
                      dispformula = ~location,
                      family = tweedie())

AIC(glmm.eph.2, glmm.eph.3)
# AIC is somewhat reduced by adding in dispersion formula, but not really

resdh.glmm.eph.2 <- simulateResiduals(glmm.eph.2)
plot(resdh.glmm.eph.1)

# the QQ plot is improved (I think) by this inclusion

summary(glmm.eph.2)

# glmm #2 wins when grazers excluded. congrats!

