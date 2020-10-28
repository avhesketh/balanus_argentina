## balanus and chthamalus models

pkgs <- c("tidyverse", "glmmTMB", "DHARMa", "car")
lapply(pkgs, install.packages, character.only = TRUE)
lapply(pkgs, library, character.only = TRUE)

se <- function(x){
  sd(x)/sqrt(length(x))
}

# read in dataframes
balanus <- read_csv("./clean_data/bio_responses.csv", col_types=
                      cols("D", "f", "f", "f", "f", "f", "n", "n", "f", "n")) %>% 
  filter(species == "Balanus_recruits")

# reorder factors for later visualization
balanus$limpets <- factor(balanus$limpets, levels = c("con","in","out"))
balanus$barnacles <- factor(balanus$barnacles, levels = c("no","yes"))
balanus$location <- factor(balanus$location, levels = c("PA","BP"))

# repeat for chthamalus barnacles
chthamalus <- read_csv("./clean_data/bio_responses.csv", col_types=
                         cols("D", "f", "f", "f", "f", "f", "n", "n", "f", "n")) %>% 
  filter(species == "Chthamalus_recruits")

chthamalus$limpets <- factor(chthamalus$limpets, levels = c("con","in","out"))
chthamalus$barnacles <- factor(chthamalus$barnacles, levels = c("no","yes"))
chthamalus$location <- factor(chthamalus$location, levels = c("PA","BP"))

# Model 8: Balanus recruitment

# start with 3-way interactions and block factor
bal.count.1 <- glmmTMB(count ~ (location + barnacles + limpets + timediff)^3 +
                        (1|block/location),
                      data = balanus,
                      family = nbinom1())

sim.res.8 <- simulateResiduals(bal.count.1)
plot(sim.res.8)
testTemporalAutocorrelation(sim.res.8)

plot(residuals(bal.count.1) ~ balanus$barnacles)
plot(residuals(bal.count.1) ~ balanus$limpets)
plot(residuals(bal.count.1) ~ balanus$location) # one pulse of recruitment in BC
plot(residuals(bal.count.1) ~ balanus$timediff) # which also pops up in time

# removing problematic timepoint from barnacle data doesn't help
# removing block term lowers AIC

bal.count.2 <- glmmTMB(count ~ (location + barnacles + limpets + timediff)^3 +
                         (1|block/location),
                      dispformula = ~location*timediff,
                      data = balanus,
                      family = nbinom1())

sim.res.8a <- simulateResiduals(bal.count.2)
plot(sim.res.8a)

# getting rid of the dispersion formula is a BAD idea though based on AIC
bal.count.3 <- glmmTMB(count ~ (location + barnacles + limpets + timediff)^3 +
                        (1|block/location),
                      data = balanus,
                      dispformula = ~timediff + location,
                      family = nbinom2())
sim.res.8b <- simulateResiduals(bal.count.3)
plot(bal.count.3)

# finally the plots look good! 
drop1(bal.count.3, test = "Chisq")

bal.count.4 <- update(bal.count.3, ~. - barnacles:limpets:timediff)

drop1(bal.count.4, test = "Chisq")

bal.count.5 <- update(bal.count.4, ~. - location:barnacles:limpets)

drop1(bal.count.5, test = "Chisq")

bal.count.6 <- update(bal.count.5, ~. - barnacles:limpets)

drop1(bal.count.6, test = "Chisq")


sim.res.8b <- simulateResiduals(bal.count.6)
plot(sim.res.8b)

summary(bal.count.6)
Anova(bal.count.6)

## Model 9: Chthamalus

chtham.count.1 <- glmmTMB(count ~ (barnacles + limpets + timediff)^3 
                     + (1|block),
                     family = nbinom2(),
                     data = chthamalus)

sim.res.9 <- simulateResiduals(chtham.count.1)
plot(sim.res.9)

plot(residuals(chtham.count.1) ~ chthamalus$barnacles)
plot(residuals(chtham.count.1) ~ chthamalus$limpets)
plot(residuals(chtham.count.1) ~ chthamalus$timediff) # time needs a dispersion term

chtham.count.2 <- glmmTMB(count ~ (barnacles + limpets + timediff)^3 
                     + (1|block),
                     family = nbinom2(),
                     disp = ~timediff,
                     data = chthamalus)

sim.res.9a <- simulateResiduals(chtham.count.2)
plot(sim.res.9a)
testTemporalAutocorrelation(sim.res.9a)
# all assumptions now met

drop1(chtham.count.2, test = "Chisq")

chtham.count.3 <- update(chtham.count.2, ~. - barnacles:limpets:timediff)

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

# Balanus glandula

bal_summary <- balanus %>% 
  group_by(location, barnacles, limpets, timediff) %>% 
  summarize(av_abund = mean(count), se_abund = sd(count)/sqrt(length(count))) %>% 
  mutate(barnacles = str_replace_all(barnacles, c("no" = "-B", "yes" = "+B")),
         limpets = str_replace_all(limpets, c("con" = "Control",
                                              "in" = "Inclusion",
                                              "out" = "Exclusion")),
         Treatment = paste(barnacles, location),
         Treatment = factor(Treatment, levels = c("+B PA",
                                                 "-B PA",
                                                 "+B BP",
                                                 "-B BP")))

bal <- ggplot(aes(x = timediff, y = av_abund, shape = Treatment, color = Treatment), 
              data = bal_summary) +
  geom_point(size = 3) +
  geom_line(aes(lty = Treatment)) +
  theme_classic() +
  ylab(expression(~italic("Balanus glandula")~ "recruit abundance")) +
  xlab("Time (weeks)") +
  scale_shape_manual(values = c(16,1,16,1)) +
  scale_linetype_manual(values = c(1,6,1,6))+
  scale_colour_manual(values = c("indianred3", "indianred3", "steelblue3", "steelblue3")) +
  facet_wrap(~limpets, nrow = 3, scales = "free", strip.position = "right") +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(plot.title = element_text(size = 16, hjust = 0.5)) +
  theme(strip.text.y = element_text(size = 14)) +
  theme(panel.spacing = unit(1.2, "lines")) +
  geom_errorbar(aes(ymin = av_abund - se_abund, ymax = av_abund + se_abund), width = 1.5)
bal

#ggsave("./figures/Figure_5.tiff", plot = bal,
#width = 7, height = 5, units = "in",
#dpi = 600)

# Chthamalus

chtham_summ <- chthamalus %>% 
  group_by(barnacles, limpets, timediff) %>% 
  summarize(av_abund = mean(count), se_abund = sd(count)/sqrt(length(count))) %>% 
  mutate(barnacles = str_replace_all(barnacles, c("no" = "-B", "yes" = "+B")),
         limpets = str_replace_all(limpets, c("con" = "Control",
                                              "in" = "Inclusion",
                                              "out" = "Exclusion"))) %>% 
  rename("Barnacles" = barnacles)

chtham <- ggplot(aes(x = timediff, y = av_abund, shape = Barnacles), 
                 data = chtham_summ) +
  geom_point(size = 3, colour = "steelblue3") +
  geom_line(aes(lty = Barnacles), colour = "steelblue3") +
  theme_classic() +
  ylab(expression(~italic("Chthamalus dalli")~ "recruit abundance")) +
  xlab("Time (weeks)") +
  scale_shape_manual(values = c(16,1)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(plot.title = element_text(size = 16, hjust = 0.5)) +
  facet_wrap(~limpets, nrow = 3, scales = "free_x", strip.position = "right")+
  theme(strip.text.y = element_text(size = 14)) +
  geom_errorbar(aes(ymin = av_abund - se_abund, ymax = av_abund + se_abund),
                width = 1.5, colour = "steelblue3")
chtham

#ggsave("./figures/Figure_S2_S3.tiff", plot = chtham,
#width = 7, height = 5, units = "in",
#dpi = 600)