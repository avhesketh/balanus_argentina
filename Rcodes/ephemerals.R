# 4: Analyzing and visualizing algae cover datasets

## clean time-adjusted model for algal cover
## after much discussion, should use 3 dates for each response
## and also use a generalized linear model with a Tweedie distribution

packages <- c("tidyverse", "glmmTMB", "DHARMa", "car", "metaDigitise")

for (i in 1:length(packages)){
  if (!require(packages[i], character.only = TRUE)) {
    install.packages(packages[i], dependencies = TRUE)
    library(packages[i], character.only = TRUE)
  }
}

se <- function(x){
  sd(x)/sqrt(length(x))
}


# load in appropriate data

ephemerals <- read_csv("./clean_data/bio_responses.csv") %>% 
  mutate(species = str_replace_all(species, c("Ephemeral Algae" = "Ephemerals"))) %>% 
  filter(species == "Ephemerals") %>% 
  # note that argentina data should be missing initial time points because of treatment artifact
  mutate(remove = if_else(location == "PA" & timediff == 0, TRUE, FALSE)) %>%
  filter(remove == FALSE) %>% 
  select(-remove)

# keep only three timepoints for each site: fall and spring equinoxes

eph_arg_keep <- ephemerals %>%
  filter(location == "PA") %>% 
  filter(date == "2006-04-26" | date == "2006-10-07" | date == "2007-02-19")
  
eph_bc_keep <- ephemerals %>%
  filter(location == "BP") %>% 
  filter(date == "2006-10-06" | date == "2007-04-06" | date == "2007-07-17")

ephemerals_keep <- eph_bc_keep %>% 
  full_join(eph_arg_keep)

# Model 5: ephemeral cover
eph.cover.1 <- glmmTMB(percent_cover ~ (timediff + location + barnacles + limpets)^2 
                             + (1|block/location),
                      data = ephemerals_keep,
                      family = tweedie())

sim.res.5 <- simulateResiduals(eph.cover.1)
plot(sim.res.5)
testTemporalAutocorrelation(sim.res.5)

# try adding in a dispersion formula

plot(residuals(eph.cover.2) ~ as.factor(ephemerals_keep$barnacles)) # some variation
plot(residuals(eph.cover.2) ~ as.factor(ephemerals_keep$location)) # some variation
plot(residuals(eph.cover.1) ~ as.factor(ephemerals_keep$limpets)) #major source of dispersion
plot(residuals(eph.cover.2) ~ ephemerals_keep$timediff) # not an issue

# add in an additive dispersion formula for barnacles and limpets
eph.cover.2 <- glmmTMB(percent_cover ~ (timediff + location + barnacles + limpets)^2 
                       + (1|block/location),
                       dispformula = ~limpets + barnacles,
                       data = ephemerals_keep,
                       family = tweedie())

sim.res.5a <- simulateResiduals(eph.cover.2)
plot(sim.res.5a)
testTemporalAutocorrelation(sim.res.5a)
# assumptions of model are satisfied

drop1(eph.cover.2, test = "Chisq")
eph.cover.3 <- update(eph.cover.2, ~. - timediff:location)

drop1(eph.cover.3, test = "Chisq")
eph.cover.4 <- update(eph.cover.3, ~. -barnacles:limpets)

drop1(eph.cover.4, test = "Chisq")

## the model is finished
summary(eph.cover.4)
Anova(eph.cover.4)
