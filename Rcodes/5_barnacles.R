## balanus and chthamalus models

pkgs <- c("tidyverse", "glmmTMB", "DHARMa", "car", "grid", "ggplotify")
lapply(pkgs, install.packages, character.only = TRUE)
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)
se <- function(x){
  sd(x)/sqrt(length(x))
}

# read in dataframes
balanus <- read_csv("./clean_data/bio_responses.csv", col_types=
                      cols("D", "f", "f", "f", "f", "f", "n", "n", "f", "n")) %>% 
  filter(species == "Balanus_recruits") %>% 
  # reorder factors for later visualization
  mutate(herbivores = factor(herbivores, levels = c("con","in","out")),
         barnacles = factor(barnacles, levels = c("no","yes")),
         location = factor(location, levels = c("BP","PA")))


# repeat for chthamalus barnacles
chthamalus <- read_csv("./clean_data/bio_responses.csv", col_types=
                         cols("D", "f", "f", "f", "f", "f", "n", "n", "f", "n")) %>% 
  filter(species == "Chthamalus_recruits") %>% 
  mutate(herbivores = factor(herbivores, levels = c("con","in","out")),
         barnacles = factor(barnacles, levels = c("no","yes")),
         location = factor(location, levels = c("BP","PA")))

# Model 8: Balanus recruitment

# start with 3-way interactions and block factor
bal.count.1 <- glmmTMB(count ~ (location + barnacles + herbivores + timediff)^3 +
                        (1|block/location),
                      data = balanus,
                      family = nbinom2())

sim.res.8 <- simulateResiduals(bal.count.1)
plot(sim.res.8)
testTemporalAutocorrelation(sim.res.8)

plot(residuals(bal.count.1) ~ balanus$barnacles)
plot(residuals(bal.count.1) ~ balanus$herbivores) # exclusion treatment has more dispersion
plot(residuals(bal.count.1) ~ balanus$location) # one pulse of recruitment in BC means residuals different between location
plot(residuals(bal.count.1) ~ balanus$timediff) # which also pops up in time

# removing problematic timepoint from barnacle data doesn't help
# removing block term raises AIC
# location*herbivores is the best dispersion formula in terms of AIC

bal.count.2 <- glmmTMB(count ~ (location + barnacles + herbivores + timediff)^3
                       + (1|block/location),
                      dispformula = ~location*herbivores,
                      data = balanus,
                      family = nbinom2())

sim.res.8a <- simulateResiduals(bal.count.2)
plot(sim.res.8a)

# the diagnostic plots all check out

drop1(bal.count.2, test = "Chisq")

bal.count.3 <- update(bal.count.2, ~. - barnacles:herbivores:timediff)

drop1(bal.count.3, test = "Chisq")

bal.count.4 <- update(bal.count.3, ~. - location:barnacles:herbivores)

drop1(bal.count.4, test = "Chisq")

# p value not significant for barnacles:herbivores, but AIC is not really improved. Retain term.

sim.res.8b <- simulateResiduals(bal.count.4)
plot(sim.res.8b)

summary(bal.count.4)
Anova(bal.count.4)

## Model 9: Chthamalus

chtham.count.1 <- glmmTMB(count ~ (barnacles + herbivores + timediff)^3 
                     + (1|block),
                     family = nbinom2(),
                     data = chthamalus)

sim.res.9 <- simulateResiduals(chtham.count.1)
plot(sim.res.9)

plot(residuals(chtham.count.1) ~ chthamalus$barnacles)
plot(residuals(chtham.count.1) ~ chthamalus$herbivores)
plot(residuals(chtham.count.1) ~ chthamalus$timediff) # time needs a dispersion term

chtham.count.2 <- glmmTMB(count ~ (barnacles + herbivores + timediff)^3 
                     + (1|block),
                     family = nbinom2(),
                     disp = ~timediff,
                     data = chthamalus)

sim.res.9a <- simulateResiduals(chtham.count.2)
plot(sim.res.9a)
testTemporalAutocorrelation(sim.res.9a)
# all assumptions now met

drop1(chtham.count.2, test = "Chisq")

chtham.count.3 <- update(chtham.count.2, ~. - barnacles:herbivores:timediff)

drop1(chtham.count.3, test = "Chisq")

chtham.count.4 <- update(chtham.count.3, ~. - barnacles:timediff)

drop1(chtham.count.4, test = "Chisq")

summary(chtham.count.4)
Anova(chtham.count.4)

sim.res.9b <- simulateResiduals(chtham.count.4)
testTemporalAutocorrelation(sim.res.9b)
plot(sim.res.9b)

# final model looks good - meets all the assumptions.

## FIGURES

barnacle_theme <- theme_classic +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(plot.title = element_text(size = 16, hjust = 0.5)) +
  theme(strip.text.y = element_text(size = 14))

# Balanus glandula

bal_summary <- balanus %>% 
  group_by(location, barnacles, herbivores, timediff) %>% 
  summarize(av_abund = mean(count), se_abund = sd(count)/sqrt(length(count))) %>% 
  mutate(barnacles = str_replace_all(barnacles, c("no" = "-B", "yes" = "+B")),
         herbivores = str_replace_all(herbivores, c("con" = "Control",
                                              "in" = "Inclusion",
                                              "out" = "Exclusion")),
         Treatment = paste(barnacles, location),
         Treatment = factor(Treatment, levels = c("+B PA",
                                                 "-B PA",
                                                 "+B BP",
                                                 "-B BP")),
         herbivores = factor(herbivores, levels = c("Control","Inclusion","Exclusion")))

bal <- ggplot(aes(x = timediff, y = av_abund, shape = Treatment, color = Treatment), 
              data = bal_summary) +
  geom_point(size = 3) +
  geom_line(aes(lty = Treatment)) +
  labs(x = "Time (weeks)", y = expression(~italic("Balanus glandula")~ "recruit abundance")) +
  scale_shape_manual(values = c(16,1,16,1)) +
  scale_linetype_manual(values = c(1,6,1,6))+
  scale_colour_manual(values = c("indianred3", "indianred3", "steelblue3", "steelblue3")) +
  facet_wrap(~herbivores, nrow = 3, scales = "free", strip.position = "right") +
  theme(panel.spacing = unit(1.2, "lines")) +
  scale_y_continuous(breaks = c(0,30,60,90,120,150)) +
  geom_errorbar(aes(ymin = av_abund - se_abund, ymax = av_abund + se_abund), width = 1.5) +
  barnacle_theme
bal

bal.table <- ggplotGrob(bal)
bal.table$heights[[7]] <- unit(0.533, "null")
bal.table$heights[[11]] <- unit(0.467, "null")
bal.scaled <- as.ggplot(bal.table)

#ggsave("./figures/Figure_5.tiff", plot = bal.scaled,
#width = 6, height = 8, units = "in", compression = "lzw",
#dpi = 800)

# Chthamalus

chtham_summ <- chthamalus %>% 
  group_by(barnacles, herbivores, timediff) %>% 
  summarize(av_abund = mean(count), se_abund = sd(count)/sqrt(length(count))) %>% 
  mutate(barnacles = str_replace_all(barnacles, c("no" = "-B", "yes" = "+B")),
         herbivores = str_replace_all(herbivores, c("con" = "Control",
                                              "in" = "Inclusion",
                                              "out" = "Exclusion"))) %>% 
  rename("Barnacles" = barnacles)

chtham <- ggplot(aes(x = timediff, y = av_abund, shape = Barnacles), 
                 data = chtham_summ) +
  geom_point(size = 3, colour = "steelblue3") +
  geom_line(aes(lty = Barnacles), colour = "steelblue3") +
  labs(y = expression(~italic("Chthamalus dalli")~ "recruit abundance"), x = "Time (weeks)") +
  scale_shape_manual(values = c(16,1)) +
  facet_wrap(~herbivores, nrow = 3, scales = "free_x", strip.position = "right")+
  geom_errorbar(aes(ymin = av_abund - se_abund, ymax = av_abund + se_abund),
                width = 1.5, colour = "steelblue3") +
  barnacles_theme

chtham

#ggsave("./figures/Figure_S7.tiff", plot = chtham,
#width = 7, height = 5, units = "in",
#dpi = 600)
