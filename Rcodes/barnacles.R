## balanus and chthamalus models

balanus <- read_csv("../data/balanus_recruitment.csv") %>% 
  select(-1)

hist(balanus$count, breaks = 100)

## negative binomial, and VERY zero inflated (or else overdispersed)

plot(count ~ location, balanus) # BC has way more dispersion going on
plot(count ~ barnacles, balanus) # look similar
plot(count ~ limpets, balanus)
plot(count ~ timediff, balanus)

# one timepoint in BC carries all the dispersion (>50 counts). Might need to omit?

glmm.bal.1 <- glmmTMB(count ~ (location + barnacles + limpets + timediff)^3 +
                        (1|block/location),
                      dispformula = ~location*timediff,
                      data = balanus,
                      family = nbinom1())
resdh.glmm.bal.1 <- simulateResiduals(glmm.bal.1)
plot(resdh.glmm.bal.1)

summary(glmm.bal.1)

# removing problematic timepoint from barnacle data doesn't help
# removing block term lowers AIC

glmm.bal.2 <- glmmTMB(count ~ (location + barnacles + limpets + timediff)^3,
                      dispformula = ~location*timediff,
                      data = balanus,
                      family = nbinom1())

resdh.glmm.bal.2 <- simulateResiduals(glmm.bal.2)
plot(resdh.glmm.bal.2)

# getting rid of the dispersion formula is a BAD idea though based on AIC
glmm.bal.3 <- glmmTMB(count ~ (location + barnacles + limpets + timediff)^3,
                      data = balanus,
                      dispformula = ~timediff + location,
                      family = nbinom2())
resdh.glmm.bal.3 <- simulateResiduals(glmm.bal.3)
plot(resdh.glmm.bal.3)

# finally the plots look good! 
drop1(glmm.bal.3, test = "Chisq")

glmm.bal.3a <- update(glmm.bal.3, ~. - barnacles:limpets:timediff)

drop1(glmm.bal.3a, test = "Chisq")

glmm.bal.3b <- update(glmm.bal.3a, ~. - location:barnacles:limpets)

drop1(glmm.bal.3b, test = "Chisq")

glmm.bal.3c <- update(glmm.bal.3b, ~. - barnacles:limpets)

drop1(glmm.bal.3c, test = "Chisq")

summary(glmm.bal.3c)

plot(simulateResiduals(glmm.bal.3c))
Anova(glmm.bal.3c)


## chthamalus now in BC

hist(chthamalus$count, breaks = 100)

glmm.ch.1 <- glmmTMB(count ~ (barnacles + limpets + timediff)^3 
                     + (1|block),
                     family = nbinom2(),
                     data = chthamalus)

resdh.glmm.ch.1 <- simulateResiduals(glmm.ch.1)
plot(resdh.glmm.ch.1)

res.glmm.ch.1 <- residuals(glmm.ch.1)

plot(res.glmm.ch.1 ~ chthamalus$barnacles)
plot(res.glmm.ch.1 ~ chthamalus$limpets)
plot(res.glmm.ch.1 ~ chthamalus$timediff)

# there's a dispersion problem, but everything seems to be working ok
# keep the block term though

glmm.ch.2 <- glmmTMB(count ~ (barnacles + limpets + timediff)^3 
                     + (1|block),
                     family = nbinom2(),
                     disp = ~timediff,
                     data = chthamalus)

AIC(glmm.ch.1, glmm.ch.2)
# adding this dispersion helps

resdh.glmm.ch.2 <- simulateResiduals(glmm.ch.2)
plot(resdh.glmm.ch.2)

acf(residuals(glmm.ch.2))

## ok it's all GOOD! Except positive autocorrelation, so let's try adding that.

glmm.ch.3 <- glmmTMB(count ~ (barnacles + limpets + timediff)^3
                     + (1|block) + diag(timediff + 0 |plot),
                     family = nbinom2(),
                     disp = ~timediff,
                     data = chthamalus)

res <- residuals(glmm.ch.3)
pacf(res)
plot(simulateResiduals(glmm.ch.3))

# not really a way to take out autocorrelation, unfortunately. it's a weird pattern
#now drop terms...

drop1(glmm.ch.2, test = "Chisq")

glmm.ch.2a <- update(glmm.ch.2, ~. - barnacles:limpets:timediff)

drop1(glmm.ch.2a, test = "Chisq")

glmm.ch.2b <- update(glmm.ch.2a, ~. - barnacles:timediff)

drop1(glmm.ch.2b, test = "Chisq")

summary(glmm.ch.2b)
Anova(glmm.ch.2b)

acf(residuals(glmm.ch.2b))
# definitely some autocorrelation ... how to fix??