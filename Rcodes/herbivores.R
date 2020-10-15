# 2: examining herbivore abundance between sites

herbs <- read_csv("../data/herbivores2.csv") %>% 
  select(-1) %>% 
  unite("plot_location", c(plot, location), remove = FALSE)

## negbin dist glmmtmb

glmm.herb.1 <- glmmTMB(count ~ (timediff + location + barnacles)^3 
                    + (1|block/location) + ar1(factor(timediff) + 0|plot_location),
                       data = herbs,
                       family = nbinom1())

resdh.glmm.herb.1 <- simulateResiduals(glmm.herb.1)
plot(resdh.glmm.herb.1)

acf(residuals(glmm.herb.1))

# residual plots ok, but try dropping block effect


glmm.herb.2 <- glmmTMB(count ~ (timediff + location + barnacles)^3 
                       + ar1(factor(timediff) + 0|plot_location),
                       data = herbs,
                       family = nbinom1())

resdh.glmm.herb.2 <- simulateResiduals(glmm.herb.2)
plot(resdh.glmm.herb.2)

acf(residuals(glmm.herb.2))

AIC(glmm.herb.1, glmm.herb.2)

# removing block effect decreases AIC, but still some dispersion problems...

res.glmm.herb.2 <- residuals(glmm.herb.2)

plot(res.glmm.herb.2 ~ herbs$timediff)
plot(res.glmm.herb.2 ~ as.factor(herbs$location))
plot(res.glmm.herb.2 ~ as.factor(herbs$barnacles))

# i'd say the most deviation occurs with time - variance highest initially, and declines with time
# but not enough problem that I want to change anything...

summary(glmm.herb.2)
drop1(glmm.herb.2, test = "Chisq")

Anova(glmm.herb.2)

# don't need to drop anything


# now for BC specific results of litt and lott abundance

lotts <- limpets_barnacles %>% 
  filter(species == "Lottia") %>% 
  filter(limpets == "con")

hist(lotts$count)

# looks like negative binomial glmmtmb

glmm.lott.1 <- glmmTMB(count ~ barnacles*timediff + (1|block),
                       data = lotts,
                       family = nbinom1())
plot(simulateResiduals(glmm.lott.1))
acf(residuals(glmm.lott.1))

drop1(glmm.lott.1, test = "Chisq") # keep everything

Anova(glmm.lott.1)

summary(glmm.lott.1)


# great! done. now littorines

litts <- limpets_barnacles %>% 
  filter(species == "Littorina")

hist(litts$count)

glmm.litt.1 <- glmmTMB(count ~ barnacles*timediff*limpets
                       + (1|block),
                       data = litts,
                       family = nbinom1)

plot(simulateResiduals(glmm.litt.1))
acf(residuals(glmm.litt.1))

drop1(glmm.litt.1, test = "Chisq") # keep everything

Anova(glmm.litt.1)
summary(glmm.litt.1)

## and now siphonaria in argentina

siphonaria <- herb_counts %>% 
  filter(species == "Siphonaria")

hist(siphonaria$count)

# maybe a negative binomial? or zero-inflated poisson? start with negbin

glmm.siph.1 <- glmmTMB(count ~ barnacles*timediff,
                       family = nbinom1,
                       data = siphonaria)

# all checking out with assumptions
res.dh.s1 <- plot(simulateResiduals(glmm.siph.1))
acf(residuals(glmm.siph.1))

summary(glmm.siph.1)

Anova(glmm.siph.1)
