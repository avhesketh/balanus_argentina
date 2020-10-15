## May 2020 post-Chris meeting

perennials <- read_csv("../data/perennials.csv") %>% 
  select(-1)

# visualize the response variable

hist(perennials$percent_cover, breaks = 100)

# definitely more of a spiky distribution than ephemerals
# but should still be fitt-able with a Tweedie family model

per_arg <- perennials %>% 
  filter(location == "GN")

per_bc <- perennials %>% 
  filter(location == "BS")

# argentina data are missing initial time points since ephemerals counted then
# are not real data ... hiding under barnacles, so trmt specific
levels(as.factor(per_arg$date))
per_arg_keep <- per_arg %>% 
  filter(date == "2006-04-26" | date == "2006-10-07" | date == "2007-02-19")
levels(as.factor(per_bc$date))
per_bc_keep <- per_bc %>%
  filter(date == "2006-10-06" | date == "2007-04-06" | date == "2007-07-17")

perennials_keep <- per_bc_keep %>% 
  full_join(per_arg_keep)

## look at the histogram again
hist(perennials_keep$percent_cover, breaks = 100)

# now standardize the model variables
# need to first encode limpet as three separate columns

# now I can create a model

plot(perennials_keep$percent_cover ~ as.factor(perennials_keep$location))
plot(perennials_keep$percent_cover ~ as.factor(perennials_keep$barnacles))
plot(perennials_keep$percent_cover ~ as.factor(perennials_keep$limpets))
plot(perennials_keep$percent_cover ~ perennials_keep$timediff)

# it looks like BC affects the number of zeroes more than anything else

per.glmm.1 <- glmmTMB(percent_cover ~ (timediff + barnacles + location + limpets)^3 +
                        (1|block/location),
                      data = perennials_keep)
summary(per.glmm.1)


resdh.glmm.per.1 <- simulateResiduals(per.glmm.1)
plot(resdh.glmm.per.1)
# four-way interactions cannot be included, sadly, due to overparameterization
# overdispersion is not terrible, but let's try adding ziformula
# also RESIDUALS ARE TERRIBLE for this initial model!!!

# first tried ziformula for location, then tried dispformula for time since 
# variance increases with time (as one would expect).

# it looks like Tweedie is not a good fit for this distribution, 
# probably because the weight at 1 is quite high
# try beta distribution with transformed proportions?

perennials_tr <- perennials_keep %>% 
  mutate(tr_pc = (percent_cover/100*(length(percent_cover)-1)+0.5)/length(percent_cover))
range(perennials_tr$tr_pc)

per.glmm.2 <- glmmTMB(tr_pc ~ (timediff + barnacles + location + limpets)^3 +
                        (1|block/location),
                      dispformula = ~timediff,
                      family = beta_family(link = "logit"),
                      data = perennials_tr)


resdh.glmm.per.2 <- simulateResiduals(per.glmm.2)
plot(resdh.glmm.per.2)

# it's OK, but not amazing. KS test shows that data barely fit this distribution

# can't fit ziformula without error message, but even so, doesn't help things

AIC(per.glmm.1, per.glmm.2)

# try dropping block term

per.glmm.3 <- glmmTMB(percent_cover ~ (timediff + barnacles + location + limpets)^3,
                      family = tweedie(),
                      data = perennials_keep)

AIC(per.glmm.1, per.glmm.3)


resdh.glm.per.3 <- simulateResiduals(per.glmm.3)
plot(resdh.glm.per.3)

# block term really does seem to matter, so keep it in. BUT residuals look way better...
# but maybe we can fix the overdispersion by adding in time again


per.glmm.4 <-  glmmTMB(percent_cover ~ (timediff + barnacles + location + limpets)^3,
                       family = tweedie(),
                       dispformula = ~(block/location),
                       data = perennials_keep)

resdh.glmm.per.4 <- simulateResiduals(per.glmm.4)
plot(resdh.glmm.per.4)

# glmm #2 has the lowest AIC

drop1(per.glmm.2)

# no large improvements occur when things are dropped, so keep in all the three-way interactions

# now move on to transforming the data for normal dist

perennials_keep$logit_pc <- logit(perennials_keep$percent_cover)

per.glmm.5 <- glmmTMB(logit_pc ~ (timediff + barnacles + location + limpets)^3
                      + (1|block/location),
                      data = perennials_keep)

# what do the residuals look like?

res.per.5 <- residuals(per.glmm.5)
plot(res.per.5 ~ as.factor(perennials_keep$barnacles))
plot(res.per.5 ~ as.factor(perennials_keep$location))
plot(res.per.5 ~ as.factor(perennials_keep$limpets))
plot(res.per.5 ~ perennials_keep$timediff)

per.glmm.6 <- glmmTMB(logit_pc ~ (timediff + barnacles + location + limpets)^3
                      + (1|block/location),
                     dispformula = ~location*barnacles,
                      data = perennials_keep)

resdh.glmm.per.6 <- simulateResiduals(per.glmm.6)
plot(resdh.glmm.per.6)

AIC(per.glmm.5, per.glmm.6)
summary(per.glmm.6)

res.glmm.per.6 <- residuals(per.glmm.6)
plot(res.glmm.per.6 ~ as.factor(perennials_keep$barnacles))
plot(res.glmm.per.6 ~ as.factor(perennials_keep$limpets))
plot(res.glmm.per.6 ~ as.factor(perennials_keep$location))
plot(res.glmm.per.6 ~ jitter(perennials_keep$timediff, 4))

Anova(per.glmm.6)

## what if we keep deleting insignificant terms?

drop1(per.glmm.6, test = "Chisq")

per.glmm.6a <- update(per.glmm.6, ~. - barnacles:location:limpets)

drop1(per.glmm.6a, test = "Chisq")

per.glmm.6b <- update(per.glmm.6a, ~. -timediff:location:limpets)

drop1(per.glmm.6b, test = "Chisq")

per.glmm.6c <- update(per.glmm.6b, ~. - timediff:barnacles:limpets)

drop1(per.glmm.6c, test = "Chisq")

per.glmm.6d <- update(per.glmm.6c, ~. - barnacles:limpets)

drop1(per.glmm.6d, test = "Chisq")

per.glmm.6e <- update(per.glmm.6d, ~. - timediff:limpets)

drop1(per.glmm.6e, test = "Chisq")

summary(per.glmm.6e)
Anova(per.glmm.6e)

acf(residuals(per.glmm.6e))
plot(simulateResiduals(per.glmm.6e))

# now just fucus in BC

#fucus_nontrace <- fucus %>% 
 # filter(prop_cover == 0 | prop_cover >= 0.013888889)

#b <- length(fucus$prop_cover)
#fucus_nontrace$tr_pc <- (fucus_nontrace$prop_cover*(b-1)+0.5)/b

# what does the response variable look like? Lots of zeroes!
hist(fucus$prop_cover, breaks = 100)
fucus$cover <- fucus$prop_cover*100
# try for an initial full model
glmm.fuc.1 <- glmmTMB(cover ~(barnacles+limpets+timediff)^3
                    + (1|block),
                    fucus, family = tweedie())

resdh.fucus.1 <- simulateResiduals(glmm.fuc.1)
plot(resdh.fucus.1)

# looks good... try taking away block term??

glmm.fuc.2 <- glmmTMB(cover ~(barnacles+limpets+timediff)^3,
                      fucus, family = tweedie())

AIC(glmm.fuc.1, glmm.fuc.2)

# NOPE. leave in that block term

drop1(glmm.fuc.1, test = "Chisq")
# and leave in all the terms currently there.

summary(glmm.fuc.1)
Anova(glmm.fuc.1)
