# 3: analyzing and modeling herbivore abundance and biomass data

##

pkgs <- c("tidyverse", "glmmTMB", "DHARMa", "car", "metaDigitise")
#lapply(pkgs, install.packages, character.only = TRUE)
lapply(pkgs, library, character.only = TRUE)

rm(pkgs)

se <- function(x){
  sd(x)/sqrt(length(x))
}

##

# create a useful subset of data
herbivores <- read_csv("./clean_data/bio_responses.csv") %>% 
  # all three of our grazer species
  filter(species == "Lottia" | species == "Siphonaria" | species == "Littorina") %>%
  # for Lottia and Siphonaria, we only want control plots where limpet densities reflect natural responses
  mutate(remove = ifelse(species == "Lottia" & limpets != "con", TRUE, FALSE)) %>% 
  filter(remove == FALSE) %>% 
  mutate(remove = ifelse(species == "Siphonaria" & limpets != "con", TRUE, FALSE)) %>% 
  filter(remove == FALSE) %>% 
  dplyr::select(-percent_cover, -remove) %>% 
  # we need a column for unique plot identifiers for modeling temporal autocorrelation later
  unite("plot_location", c(plot, location), remove = FALSE)

#write_csv(herbivores, "./clean_data/herbivore_abundance.csv")

# to compare between sites, we need to convert to biomass

# read in shell length data for each species measured from photos
shell_length <- read_csv("./raw_data/herbivores_sl.csv") %>% 
  select(-limpets, -barnacles, -date) %>% 
  mutate(sp = str_replace_all(sp, c("siphonaria" = "Siphonaria", 
                                    "lsc" = "Littorina",
                                    "ldig" = "Lottia"))) %>% 
  rename(size_mm = size)

# building models for relationship of shell length and dry weight for each species
siph <- read_csv("./raw_data/siphonaria_dwsl.csv") %>% 
  rename(dw_mg = dw, sl_mm = sl)
# these data are extracted from Toblado & Gappa 2001
siph_dw_sl <- lm(log(dw_g) ~ log(sl_mm), data = siph)
summary(siph_dw_sl)

# create a function from the output
siph_dw <- function(len) {
  log_dw = 2.8318*log(len) - 10.6002
}

# need to extract the data from existing figures
#lott_litt <- metaDigitise("./raw_data", summary = FALSE)
# littorine data from North 1954
#litt <- do.call(rbind, lott_litt$scatterplot[1]) %>% 
  # get length in terms of mm, not cm
  #mutate(x = x*10) %>% 
  #rename(sl_mm = x, dw_g = y)

litt_dw_sl <- lm(log(dw_g) ~ log(sl_mm), data = litt)
summary(litt_dw_sl)

litt_dw <- function(len) {
  log_dw = 2.45328*log(len) - 9.67482
}

# lottia data from Frank 1965
#lott <- do.call(rbind, lott_litt$scatterplot[2]) %>% 
  # convert length to mm, convert volume to dry weight 
  #mutate(x = x*10, y = y*0.35) %>% 
  #rename(sl_mm = x, dw_g = y)

lott_dw_sl <- lm(log(dw_g) ~ log(sl_mm), data = lott)
summary(lott_dw_sl)

# create function from model outputs
lottia_dw <- function(len) {
  log_dw <- 1.55329*log(len) - 6.51138
}

# now for biomass estimation
set.seed(26)
# create an empty column to receive dry weights
herbivores$total_dw_g <- NA

for (x in 1:length(herbivores$species)){
  # subset shell lengths of the appropriate species for sampling
  sub <- subset(shell_length$size_mm, shell_length$sp == herbivores$species[x])
  # and sample this subset
  y <- sample(sub, size = herbivores$count[x], replace = TRUE)
  # if the species is siphonaria, perform the appropriate function on the sample
  if (herbivores$species[x] == "Siphonaria"){
    herbivores$total_dw_g[x] <- sum(exp(siph_dw(y)))
  }
  # and same for littorines
  if (herbivores$species[x] == "Littorina"){
    herbivores$total_dw_g[x] <- sum(exp(litt_dw(y)))
  }
  # and for lottia
  if (herbivores$species[x] == "Lottia") {
    herbivores$total_dw_g[x] <- sum(exp(siph_dw(y)))
  }
}

# create a new dataframe for total biomass
herb_biomass <- herbivores %>% 
  group_by(location, timediff, plot, block, barnacles) %>% 
  summarize(biomass_g = sum(total_dw_g)) %>% 
  unite(plot_location, c(plot,location), remove = FALSE)

#write_csv(herb_biomass, "./clean_data/herbivore_biomass.csv")
#herb_biomass$timediff <- numFactor(herb_biomass$timediff)

# Model 1: biomass between locations

#herb_biomass <- read_csv("./clean_data/herbivore_biomass.csv")

herb.bio.1 <- glmmTMB(biomass_g ~ (location+barnacles+timediff)^3 +
                        (1|block/location),
                   data = herb_biomass,
                   family = tweedie()
                   )

sim.res.1 <- simulateResiduals(herb.bio.1)
testTemporalAutocorrelation(simulationOutput = sim.res.1, time = unique(herb_biomass$timediff))
# based on test of temporal autocorrelation, no additional term is needed
plot(sim.res.1)
# need to continue tweaking the model - there are some dispersion issues with the residuals

plot(residuals(herb.bio.1) ~ as.factor(herb_biomass$barnacles)) # a lot of variation
plot(residuals(herb.bio.1) ~ as.factor(herb_biomass$location)) # some variation here too
plot(residuals(herb.bio.1) ~ (herb_biomass$timediff)) 

herb.bio.2 <- glmmTMB(biomass_g ~ (location+barnacles+timediff)^3 +
                      (1|block/location),
                    data = herb_biomass,
                    dispformula = ~barnacles,
                    family = tweedie())

sim.res.1a <- simulateResiduals(herb.bio.2)
testTemporalAutocorrelation(simulationOutput = sim.res.1a, time = unique(herb_biomass$timediff))
# based on test of temporal autocorrelation, no additional term is needed
plot(sim.res.1a)
# looks good! adding the dispersion formula solved the issues with the model.
# slight pattern in residuals, but looks relatively good

summary(herb.bio.2)
Anova(herb.bio.2)

drop1(herb.bio.2, test = "Chisq")
herb.bio.3 <- update(herb.bio.2, ~. -location:barnacles:timediff)
drop1(herb.bio.3, test = "Chisq")
herb.bio.4 <- update(herb.bio.3, ~. -location:barnacles)
drop1(herb.bio.4, test = "Chisq")
herb.bio.5 <- update(herb.bio.4, ~. -barnacles:timediff)
drop1(herb.bio.5, test = "Chisq")

sim.res.1b <- simulateResiduals(herb.bio.5)
plot(sim.res.1b)
summary(herb.bio.5)
Anova(herb.bio.5)

# all the terms have been dropped that needed to be

## Model 2: Count data for Lottia at BP

lott_abund <- herbivores %>% 
  filter(species == "Lottia")

# Lottia counts looks like they have a negative binomial distribution
# random effect of block within location
glimpse(lott_abund)

lott.abund.1 <- glmmTMB(count ~ barnacles*timediff + (1|block),
                       data = lott_abund,
                       family = nbinom1())
# dropping the block effect does not help AIC at all, so it stays in
sim.res.2 <- simulateResiduals(lott.abund.1)
# plot the residuals to check that assumptions are met; they aren't completely
plot(sim.res.2)
# temporal autocorrelation check; there is none
testTemporalAutocorrelation(simulationOutput = sim.res.2)
# the residuals have an odd pattern to them that may indicate dispersion formula is needed

plot(residuals(lott.abund.1)~ lott_abund$timediff)
plot(residuals(lott.abund.1)~ as.factor(lott_abund$barnacles))
# based on these residuals plots, there looks to be greater variance at certain timepoints
# add in a dispersion formula

lott.abund.2 <- glmmTMB(count ~ barnacles*timediff + (1|block),
                        data = lott_abund,
                        dispformula = ~timediff,
                        family = nbinom1())
sim.res.2a <- simulateResiduals(lott.abund.2)
plot(sim.res.2a)
# this looks much better now that dispersion through time is accounted for. keep this model.

# reduce model to minimize AIC

drop1(lott.abund.2, test = "Chisq") 

lott.abund.3 <- update(lott.abund.2, ~. -barnacles:timediff)

drop1(lott.abund.3, test = "Chisq")

# AIC is minimized.Double-check the assumptions.

sim.res.2b <- simulateResiduals(lott.abund.3)
plot(sim.res.2b)
testTemporalAutocorrelation(sim.res.2b)

# no assumptions appear to be violated. The final model:
print(lott.abund.3)
summary(lott.abund.3)
Anova(lott.abund.3)

## 3: Count data for Littorina at BP
litt_abund <- herbivores %>% 
  filter(species == "Littorina")

# response data are negative binomial. Start with a full model.
litt.abund.1 <- glmmTMB(count ~ barnacles*timediff*limpets
                       + (1|block),
                       data = litt_abund,
                       family = nbinom1)
sim.res.3 <- simulateResiduals(litt.abund.1)
plot(sim.res.3)
testTemporalAutocorrelation(sim.res.3)

# no temporal autocorrelation, but we may need a dispersion formula

plot(residuals(litt.abund.1)~ litt_abund$timediff)
plot(residuals(litt.abund.1)~ as.factor(litt_abund$barnacles))

# dispersion is less a factor of time, and more a factor of barnacles. introduce a dispersion formula

litt.abund.2 <- glmmTMB(count ~ barnacles*timediff*limpets
                        + (1|block),
                        data = litt_abund,
                        dispformula = ~barnacles,
                        family = nbinom1)
sim.res.3a <- simulateResiduals(litt.abund.2)
plot(sim.res.3a)
testTemporalAutocorrelation(sim.res.3a)

# adding a dispersion formula has solved the pattern in the residuals plot

drop1(litt.abund.2, test = "Chisq") # drop the three-way interaction
litt.abund.3 <- update(litt.abund.2, ~. -barnacles:timediff:limpets)
drop1(litt.abund.3, test = "Chisq") # don't drop anything else

summary(litt.abund.3)
Anova(litt.abund.3)

# Model 4: Count data for Siphonaria at PA

siph_abund <- herbivores %>% 
  filter(species == "Siphonaria")

# negative binomial distribution, no block effect needed

siph.abund.1 <- glmmTMB(count ~ barnacles*timediff,
                       family = nbinom1,
                       data = siph_abund)

# all checking out with assumptions
sim.res.4 <- simulateResiduals(siph.abund.1)
plot(sim.res.4)
testTemporalAutocorrelation(sim.res.4)
# all the assumptions of the tests check out

drop1(siph.abund.1, test = "Chisq")

summary(siph.abund.1)
Anova(siph.abund.1)

# script for creating figures for herbivory responses

# herbivore biomass by location and barnacle treatment (limpet control plots only)
herb_summary <- read_csv("./clean_data/herbivore_biomass.csv") %>% 
  group_by(location, barnacles, timediff) %>% 
  summarize(av_bi = mean(biomass_g), se_bi = se(biomass_g)) %>% 
  mutate(barnacles = str_replace_all(barnacles, c("no" = "-B",
                                                  "yes" = "+B")),
         Treatment = paste(barnacles,location))

herb_summary$Treatment <- factor(herb_summary$Treatment, 
                                  levels = c(
                                    "+B BP",
                                    "-B BP",
                                    "+B PA",
                                    "-B PA"))

herbs <- ggplot(aes(x = timediff, y = av_bi, shape = Treatment,
                    color = Treatment), 
                data = herb_summary) +
  geom_point(size = 3) +
  theme_classic() +
  geom_line(aes(lty = Treatment), size = 1) +
  ylab("Herbivore biomass (g dry tissue)") +
  xlab("Time (weeks)") +
  scale_shape_manual(values = c(16,1,16,1)) +
  scale_linetype_manual(values = c(1,6,1,6))+
  scale_colour_manual(values = c("steelblue3", "steelblue3","indianred3", "indianred3")) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  theme(axis.title.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18)) +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.title = element_text(size = 18)) +
  theme(legend.position = "top") +
  geom_errorbar(aes(ymin = av_bi - se_bi, ymax = av_bi + se_bi), width = 1.5)
herbs

#ggsave("./figures/Figure_2c.tiff", plot = herbs, 
#width = 9, height = 5, units = "in",
#dpi = 600)

# littorine and limpet abundance plots

litt_summ <- read_csv("./clean_data/herbivore_abundance.csv") %>% 
  filter(species == "Littorina") %>% 
  group_by(barnacles, limpets, timediff) %>% 
  summarize(av_ab = mean(count), se_ab = se(count)) %>% 
  rename("Barnacles" = barnacles) %>% 
  mutate(Barnacles = str_replace_all(Barnacles, c("no" = "-B",
                                                  "yes" = "+B")),
         limpets = str_replace_all(limpets, c("con" = "Control",
                                              "in" = "Inclusion",
                                              "out" = "Exclusion"))) %>% 
  mutate(limpets = factor(limpets, levels = c("Control","Inclusion","Exclusion")))

labels.litt <- data.frame(label = c("Control", "Inclusion", "Exclusion"),
                          limpets = as.factor(c("Control", "Inclusion", "Exclusion")))


litt_plot <- ggplot(aes(x = timediff, y = av_ab), 
                    data = litt_summ) +
  geom_point(size = 3, colour = "steelblue3", aes(shape = Barnacles)) +
  geom_line(colour = "steelblue3", aes(lty = Barnacles)) +
  theme_classic() +
  scale_linetype_manual(values = c(6,1)) +
  ylab(expression("Abundance of"~italic("Littorina")~ "spp.")) +
  xlab("Time (weeks)") +
  scale_shape_manual(values = c(1,16)) +
  theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18)) +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 20, hjust = 0.5)) +
  theme(strip.text.y = element_text(size = 18)) +
  facet_wrap(~limpets, nrow = 3, scales = "free_x") +
  theme(strip.text = element_blank(), strip.background = element_blank()) +
  geom_errorbar(aes(ymin = av_ab - se_ab, ymax = av_ab + se_ab), colour = "steelblue3", width = 1.5) +
  geom_text(aes(x = 45, y = 25, label = label), data = labels.litt, size = 7,
            fontface = "bold")
litt_plot

#ggsave("./figures/Figure_2a.tiff", plot = litt_plot, 
       #width = 6.5, height = 5, units = "in",
       #dpi = 600)

# limpets and "limpets", control treatment only

limpet_summ <- read_csv("./clean_data/herbivore_abundance.csv") %>% 
  filter(species == "Lottia" | species == "Siphonaria") %>% 
  group_by(barnacles, location, species, timediff) %>% 
  summarize(av_ab = mean(count), se_ab = se(count)) %>% 
  rename("Barnacles" = barnacles) %>%  
  mutate(Barnacles = str_replace_all(Barnacles, c("no" = "-B",
                                                  "yes" = "+B")),
         Treatment = paste(Barnacles, location))

limpet_summ$Treatment <- factor(limpet_summ$Treatment,
                              levels = c("+B PA",
                                         "-B PA",
                                         "+B BP",
                                         "-B BP"))

limpet_plot <- ggplot(aes(x = timediff, y = av_ab, colour = Treatment, shape = Treatment), data = limpet_summ) +
  geom_point(size = 3) +
  geom_line(aes(lty = Treatment)) +
  scale_shape_manual(values = c(16,1,16,1)) +
  scale_linetype_manual(values = c(1,6,1,6)) +
  scale_colour_manual(values = c("indianred3", "indianred3", "steelblue3", "steelblue3")) +
  theme_classic() +
  ylab(expression("Abundance of"~italic("Lottia")~"spp./"~italic("S. lessonii"))) +
  xlab("Time (weeks)") +
  theme(axis.text.x = element_text(size = 18)) +
  theme(axis.text.y = element_text(size = 18)) +
  theme(axis.title.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 19.5)) +
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 20, hjust = 0.5)) +
  geom_errorbar(aes(ymin = av_ab - se_ab, ymax = av_ab + se_ab), width = 1.5)
limpet_plot

#ggsave("./figures/Figure_2b.tiff", plot = limpet_plot, 
#width = 6.5, height = 5, units = "in",
#dpi = 600)