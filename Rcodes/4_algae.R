# 4: Analyzing and visualizing algae cover datasets

## clean time-adjusted model for algal cover
## after much discussion, should use 3 dates for each response
## and also use a generalized linear model with a Tweedie distribution

pkgs <- c("tidyverse", "glmmTMB", "DHARMa", "car")
lapply(pkgs, install.packages, character.only = TRUE)
lapply(pkgs, library, character.only = TRUE)

se <- function(x){
  sd(x)/sqrt(length(x))
}


# load in appropriate data

ephemerals <- read_csv("./clean_data/bio_responses.csv", col_types=
                         cols("D", "f", "f", "f", "f", "f", "n", "n", "f", "n")) %>% 
  filter(species == "Ephemerals") %>% 
  # note that argentina data should be missing initial time points because of treatment artifact
  mutate(remove = if_else(location == "PA" & timediff == 0, TRUE, FALSE)) %>%
  filter(remove == FALSE) %>% 
  select(-remove) %>% 
  mutate(herbivores = factor(herbivores, levels = c("con","in","out")),
          barnacles = factor(barnacles, levels = c("no","yes")))

# keep only three timepoints for each site: fall and spring equinoxes

eph_arg_keep <- ephemerals %>%
  filter(location == "PA") %>% 
  filter(date == "2006-04-26" | date == "2006-10-07" | date == "2007-02-19")
  
eph_bc_keep <- ephemerals %>%
  filter(location == "BP") %>% 
  filter(date == "2006-10-06" | date == "2007-04-06" | date == "2007-07-17")

ephemerals_eq<- eph_bc_keep %>% 
  full_join(eph_arg_keep) 

# Model 5: ephemeral cover

eph.cover.1 <- glmmTMB(percent_cover ~ (herbivores + barnacles + location + timediff)^2
                             + (1|block/location),
                      data = ephemerals_eq,
                      family = tweedie())

sim.res.5 <- simulateResiduals(eph.cover.1)
plot(sim.res.5)
testTemporalAutocorrelation(sim.res.5)

# the assumptions are met, and temporal autocorrelation looks to be a non-issue

drop1(eph.cover.1, test = "Chisq")
eph.cover.2 <- update(eph.cover.1, ~. -location:timediff)
drop1(eph.cover.2, test = "Chisq")

# stop here. deleting further terms does not substantially benefit AIC or models do not converge

sim.res.5a <- simulateResiduals(eph.cover.2)
plot(sim.res.5a)

# model assumptions are still met

summary(eph.cover.2)
Anova(eph.cover.2)

# Model 6:perennial cover

perennials <- read_csv("./clean_data/bio_responses.csv", col_types=
  cols("D", "f", "f", "f", "f", "f", "n", "n", "f", "n")) %>% 
  filter(species == "Perennials") %>% 
  mutate(herbivores = factor(herbivores, levels = c("con","in","out")))


# need to subset dates because the data are overly complex and non-linear with respect to time
per_arg_keep <- perennials %>%
  filter(location == "PA") %>% 
  filter(date == "2006-04-26" | date == "2006-10-07" | date == "2007-02-19")

per_bc_keep <- perennials %>%
  filter(location == "BP") %>% 
  filter(date == "2006-10-06" | date == "2007-04-06" | date == "2007-07-17")

perennials_eq <- per_arg_keep %>% 
  full_join(per_bc_keep)  %>% 
  mutate(barnacles = factor(barnacles, c("no","yes")))

per.cover.1 <- glmmTMB(percent_cover ~ (timediff + barnacles + location + herbivores)^3
                        ,
                       dispformula = ~location,
                      family = tweedie(),
                      data = perennials)

sim.res.6 <- simulateResiduals(per.cover.1)
plot(sim.res.6)
testTemporalAutocorrelation(sim.res.6)

# the distribution seems potentially wrong! try working with a different error distribution

# logit transform and work with gaussian distribution

per.cover.2 <- glmmTMB(logit(percent_cover) ~ (timediff + barnacles + location + herbivores)^3
                       + (1|block/location),
                       # add in dispersion formula based on residuals plots
                       dispformula = ~location*barnacles,
                       family = gaussian(),
                       data = perennials_eq)

sim.res.6a <- simulateResiduals(per.cover.2)
plot(sim.res.6a)
testTemporalAutocorrelation(sim.res.6a)
# assumptions are satisfied

# iteratively drop terms to minimize AIC
drop1(per.cover.2, test="Chisq")
per.cover.3 <- update(per.cover.2, ~. -barnacles:location:herbivores)
drop1(per.cover.3, test = "Chisq")
per.cover.4 <- update(per.cover.3, ~. -timediff:location:herbivores)
drop1(per.cover.4, test = "Chisq")
per.cover.5 <- update(per.cover.4, ~. - timediff:barnacles:herbivores)
drop1(per.cover.5, test = "Chisq")
per.cover.6 <- update(per.cover.5, ~. -barnacles:herbivores)
drop1(per.cover.6, test = "Chisq")
per.cover.7 <- update(per.cover.6, ~. -timediff:herbivores)
drop1(per.cover.7, test = "Chisq")
# finished dropping terms

sim.res.6b <- simulateResiduals(per.cover.7)
plot(sim.res.6b)
testTemporalAutocorrelation(sim.res.6b)
# updated model is simplified and still meets all assumptions

summary(per.cover.7)
Anova(per.cover.7)

# Model 7: Fucus in BC

fucus <- read_csv("./clean_data/fucus_clean.csv")

fuc.cover.1 <- glmmTMB(percent_cover ~(barnacles+herbivores+timediff)^3
                      + (1|block),
                      fucus, family = tweedie())

sim.res.7 <- simulateResiduals(fuc.cover.1)
plot(sim.res.7)
testTemporalAutocorrelation(fuc.cover.1)
# satisfies the assumptions of the test, no temporal autocorrelation
# dropping the block term increases AIC

drop1(fuc.cover.1, test = "Chisq")
# dropping the three-way interaction raises AIC, so leave it as a full model

summary(fuc.cover.1)
Anova(fuc.cover.1)

# FIGURES

# main plotting theme:

algae_theme <- theme_classic + 
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 14)) +

# ephemeral algae

ealgae_summary <- ephemerals %>% 
  group_by(location, barnacles, herbivores, timediff) %>% 
  summarize(av_cover = mean(percent_cover), se_cover = se(percent_cover)) %>%
  # rename treatments
  mutate(herbivores = str_replace_all(herbivores, c("con" = "Control",
                                              "in" = "Inclusion",
                                              "out" = "Exclusion")),
         barnacles = str_replace_all(barnacles, c("no" = "-B",
                                                  "yes" = "+B")),
         herbivores = factor(herbivores, levels = c("Control","Inclusion","Exclusion")),
        Treatment = paste(barnacles, location))

# separate analysed data from nonanalysed to make visual distinction

levels <- levels(as.factor(ephemerals_eq$timediff))
# separate out the analyzed timepoints from the dataset
ealgae_analysed <- ealgae_summary %>% 
  filter(timediff %in% levels)
# and create a new dataset where these same timepoints are removed
ealgae_nonan <- ealgae_summary %>% 
  anti_join(ealgae_analysed)
  
ealgae_summary$Treatment <- factor(ealgae_summary$Treatment, 
                                   levels = c("+B PA",
                                              "-B PA",
                                              "+B BP",
                                              "-B BP"))
col.names.eph <- colnames(ealgae_summary)
# add in dummy point to bump up y axis to 5 for control limpet plot
df.eph <- data.frame( values = c("BP","-B","Control",0,5,0,"dummy"), names  = col.names.eph) %>% 
  pivot_wider(names_from = names, values_from = values) %>% 
  mutate(timediff = as.numeric(timediff), av_cover = as.numeric(av_cover),
         se_cover = as.numeric(se_cover))
ealgae_summary <- ealgae_summary %>% 
  full_join(df.eph) %>% 
  mutate(Treatment = factor(Treatment, levels = c("+B PA","-B PA","+B BP","-B BP","dummy")),
         herbivores = factor(herbivores, levels = c("Control","Inclusion","Exclusion")))

ealgae <- ggplot(aes(x = timediff, y = av_cover, shape = Treatment, color = Treatment), 
                 data = ealgae_nonan) +
  geom_line(aes(lty=Treatment), data = ealgae_summary) +
  geom_point(size = 3, colour = "grey80") +
  geom_point(size = 3, data = ealgae_analysed) +
  scale_shape_manual(values = c(16,1,16,1, 2)) +
  scale_linetype_manual(values = c(1,6,1,6, 2))+
  scale_colour_manual(values = c("indianred3", "indianred3", "steelblue3", "steelblue3", "white")) +
  labs(y = "Ephemeral algal cover (%)", x = "Time (weeks)") +
  xlim(c(0,63)) +
  algae_theme +
  theme(panel.spacing = unit(1.2, "lines")) +
  scale_y_continuous(breaks = c(0,5,10,20,30,40,50,60,70,80)) +
  facet_wrap(~herbivores, scales = "free", strip.position = "right", nrow = 3)+
  geom_errorbar(colour = "grey80", aes(ymin = av_cover - se_cover, ymax = av_cover + se_cover), width = 1.5) +
  geom_errorbar(aes(ymin = av_cover - se_cover, ymax = av_cover + se_cover), data = ealgae_analysed, width = 1.5) 
ealgae

# adjust length of y axes in panels
eph.table <- ggplotGrob(ealgae)
eph.table$heights[[7]] <- unit(0.21, "null")
eph.table$heights[[11]] <- unit(0.5, "null")
eph.scaled <- as.ggplot(eph.table)
eph.scaled

#ggsave("./figures/Figure_3.tiff", plot = eph.scaled, 
#width = 6, height = 8, units = "in", compression = "lzw",
#dpi = 800)

# Perennial algae

palgae_summary <- perennials %>%
  mutate(logit_pc = car::logit(percent_cover)) %>% 
  group_by(location, barnacles, timediff, herbivores) %>% 
  summarize(log_av_cover = mean(logit_pc), log_se_cover = se(logit_pc),
            av_cover = mean(percent_cover), se_cover = se(percent_cover)) %>% 
  mutate(herbivores = str_replace_all(herbivores, c("con" = "Control",
                                              "in" = "Inclusion",
                                              "out" = "Exclusion")),
         barnacles = str_replace_all(barnacles, c("no"="-B",
                                                  "yes"="+B")),
         herbivores = factor(herbivores, levels = c("Control","Inclusion","Exclusion")),
         Treatment = paste(barnacles, location))

palgae_summary$Treatment <- factor(palgae_summary$Treatment, 
                                   levels = c("+B PA",
                                              "-B PA",
                                              "+B BP",
                                              "-B BP"))

# separate analysed data from nonanalysed

levels <- levels(as.factor(perennials_eq$timediff))

palgae_analysed <- palgae_summary %>% 
  filter(timediff %in% levels)

palgae_nonan <- palgae_summary %>% 
  anti_join(palgae_analysed)

palgae <- ggplot(aes(x = timediff, y = av_cover, shape = Treatment,
                     color = Treatment), 
                 data = palgae_nonan) +
  geom_line(aes(lty=Treatment), data = palgae_summary) +
  geom_point(size = 3, colour = "grey80") +
  geom_point(size = 3, data = palgae_analysed) +
  scale_colour_manual(values = c("indianred3", "indianred3", "steelblue3", "steelblue3")) +
  labs(y = "Perennial algal cover (%)", x = "Time (weeks)") +
  xlim(c(0,63)) +
  ylim(c(0,100)) +
  scale_shape_manual(values = c(16,1,16,1)) +
  scale_linetype_manual(values = c(1,6,1,6))+
  facet_wrap(~herbivores, nrow = 3, strip.position = "right", scales = "free_x") +
  algae_theme +
  theme(panel.spacing = unit(1.2, "lines")) +
  geom_errorbar(colour = "grey80", aes(ymin = av_cover - se_cover
                                       , ymax = av_cover + se_cover), width = 1.5) +
  geom_errorbar(aes(ymin = av_cover - se_cover, 
                    ymax = av_cover + se_cover), data = palgae_analysed, width = 1.5) 
palgae

#ggsave("./figures/Figure_4.tiff", plot = palgae, 
#width = 6, height = 8, units = "in", compression = "lzw",
#dpi = 800)

# Figure S5

lpalgae <- ggplot(aes(x = timediff, y = log_av_cover, shape = Treatment,
                      color = Treatment), 
                  data = palgae_nonan) +
  geom_line(aes(lty=Treatment), data = palgae_summary) +
  geom_point(size = 3, colour = "grey80") +
  geom_point(size = 3, data = palgae_analysed) +
  scale_colour_manual(values = c("indianred3", "indianred3", "steelblue3", "steelblue3")) +
  labs(y = "Logit [Perennial algal cover (%)]", x = "Time (weeks)") +
  xlim(c(0,70)) +
  scale_shape_manual(values = c(16,1,16,1)) +
  scale_linetype_manual(values = c(1,6,1,6))+
  facet_wrap(~herbivores, nrow = 3, strip.position = "right", scales = "free_x") +
  algae_theme +
  theme(panel.spacing = unit(1.2, "lines")) +
  geom_errorbar(colour = "grey80", aes(ymin = log_av_cover - log_se_cover
                                       , ymax = log_av_cover + log_se_cover), width = 1.5) +
  geom_errorbar(aes(ymin = log_av_cover - log_se_cover, 
                    ymax = log_av_cover + log_se_cover), data = palgae_analysed, width = 1.5) 
lpalgae

#ggsave("./figures/Figure_S5.tiff", plot = lpalgae, 
#width = 7, height = 5, units = "in",
#dpi = 600)

# fucus at BP
glimpse(fucus)
fucus_summ <- fucus %>% 
  group_by(barnacles, herbivores, timediff) %>% 
  summarize(av_cover = mean(percent_cover), se_cover = sd(percent_cover)/sqrt(length(percent_cover))) %>% 
  mutate(herbivores = str_replace_all(herbivores, c("con" = "Control",
                                                      "in" = "Inclusion",
                                                      "out" = "Exclusion")),
         barnacles = str_replace_all(barnacles, c("no"="-B",
                                                  "yes"="+B")))

fucus_summ$herbivores <- factor(fucus_summ$herbivores, levels = c("Control", "Inclusion", "Exclusion"))

fucus_summ$Barnacles <- factor(fucus_summ$barnacles,
                               levels = c("+B","-B"))

fucus_plot <- ggplot(aes(x = timediff, y = av_cover, shape = Barnacles), 
                     data = fucus_summ) +
  geom_point(size = 3, colour = "steelblue3") +
  geom_line(aes(lty = Barnacles), colour = "steelblue3") +
  scale_linetype_manual(values = c(1,6)) +
  labs(x = "Time (weeks)", y = expression(~italic("Fucus distichus")~ "cover (%)")) +
  scale_shape_manual(values = c(16,1)) +
  algae_theme +
  theme(plot.title = element_text(size = 16, hjust = 0.5)) +
  facet_wrap(~herbivores, nrow = 3, strip.position = "right", scales = "free_x") +
  geom_errorbar(aes(ymin = av_cover - se_cover, ymax = av_cover + se_cover), 
                width = 1.5, colour = "steelblue3")
fucus_plot

#ggsave("./figures/Figure_S6.tiff", plot = fucus_plot,
      #width = 7, height = 5, units = "in",
      #dpi = 600)
