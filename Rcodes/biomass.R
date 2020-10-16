# 2: analysis and visualization of herbivore abundance & biomass data

# load packages required for script
packages <- c("tidyverse", "glmmTMB","car", "dplyr", "DHARMa")

for (i in 1:length(packages)){
if (!require(packages[i], character.only = TRUE)) {
  install.packages(packages[i], dependencies = TRUE)
  library(packages[i], character.only = TRUE)
}
}

# first, examine the herbivore count data
herb_counts <- read_csv("./clean_data/herbivores.csv", index_col = FALSE) 
# 1: Lottia spp. abundance model
lottia <- herb_counts %>% 
  filter(species == "Lottia")

# response variable has a negative binomial distribution, so use a generalized linear model
# random effect of block within location
# random effect for temporal autocorrelation

glmm.lott.1 <- glmmTMB(count ~ barnacles*timediff + (1|block),
                         data = lottia,
                         family = nbinom1())

# plot the residuals to check that assumptions are met; they are
plot(simulateResiduals(glmm.lott.1))

# and check the residuals for temporal autocorrelation; there is none
par(mfrow = c(1,1))
acf(residuals(glmm.lott.1))

# check to see if dropping the random effects term reduces AIC
glmm.lott.2 <- glmmTMB(count ~ barnacles*timediff,
                       data = lottia,
                       family = nbinom1())

AIC(glmm.lott.1, glmm.lott.2) # keep in the random effect

# determine if the full model is the most parsimonious
drop1(glmm.lott.1, test = "Chisq")

# we can drop the second-order interaction to reduce AIC
glmm.lott.3 <- update(glmm.lott.1, ~. -barnacles:timediff)

drop1(glmm.lott.3, test = "Chisq") # this model is the most parsimonious

plot(simulateResiduals(glmm.lott.3)) 
acf(residuals(glmm.lott.3)) # assumptions met, autocorrelation absent

summary(glmm.lott.3)
Anova(glmm.lott.3)

# 2: Littorina spp. abundance model
littorina <- herb_counts %>% 
  filter(species == "Littorina")

glmm.litt.1 <- glmmTMB(count ~ barnacles*timediff*limpets
                       + (1|block),
                       data = littorina,
                       family = nbinom1)

plot(simulateResiduals(glmm.litt.1))
acf(residuals(glmm.litt.1))

drop1(glmm.litt.1, test = "Chisq") # keep everything

Anova(glmm.litt.1)
summary(glmm.litt.1)

##

herb_counts$total_dw <- NA

lengths <- read_csv("../data/herbivores_sl.csv")

lengths$sp <- gsub("siphonaria", "Siphonaria", lengths$sp)
lengths$sp <- gsub("lsc", "Littorina", lengths$sp)
lengths$sp <- gsub("ldig", "Lottia", lengths$sp)

# get the herbivore dataframe up in here

herb_counts <- read_csv("../data/herbivore_counts.csv") %>% 
  select(-1)

herb_counts <- herb_counts %>% 
  na.omit()

# define relationships of biomass to organism length

lottia_dw <- function(len) {
  log_dw <- 2.45328*log(len) + 0.57312
}

litt_dw <- function(len) {
  log_dw = 2.45328*log(len) - 9.67482
}

siph_dw <- function(len) {
  log_dw = 2.8318*log(len) - 10.6002
}



set.seed(26)

for (x in 1:length(herb_counts$count)){
  sub <- subset(lengths$size, lengths$sp == herb_counts$species[x])
  y <- sample(sub, size = herb_counts$count[x], replace = TRUE)
  if (herb_counts$species[x] == "Siphonaria"){
    herb_counts$total_dw[x] <- sum(exp(siph_dw(y)))
  }
  if (herb_counts$species[x] == "Littorina"){
    herb_counts$total_dw[x] <- sum(exp(litt_dw(y)))
  }
  if (herb_counts$species[x] == "Lottia") {
    herb_counts$total_dw[x] <- sum(exp(siph_dw(y)))
  }
}

# biomass is all figured out! all above zero, all in grams
# amalgamate litt and lott for BC
# somehow there was a duplicate of siphonaria for Arg which was inflating the initial timept
# manually removed with slice for now

herb_biomass <- herb_counts %>% 
  slice(-(372:383)) %>% 
  group_by(location, timediff, plot, block, barnacles) %>% 
  summarize(biomass = sum(total_dw)) %>% 
  unite(plot_location, c(plot,location), remove = FALSE)

write.csv(herb_biomass, "../data/herb_biomass.csv")

# stats!
hist(herb_biomass$biomass)

mod.herb.dw <- glmmTMB(biomass ~ location*barnacles*timediff +
                     (1 | block/location) + 
                     ar1(factor(timediff) + 0 | plot_location ), 
                   data = herb_biomass,
                   family = tweedie)

summary(mod.herb.dw)
Anova(mod.herb.dw)
