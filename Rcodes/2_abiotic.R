# 1: Analysis and visualization of temperature and humidity data


## First load packages of interest

pkgs <- c("tidyverse")
lapply(pkgs, install.packages, character.only = TRUE)
lapply(pkgs, library, character.only = TRUE)

# define standard error function for eventual plotting
se <- function(x){
  sd(x)/sqrt(length(x))
}

# First, read in dataframes and calculate the time since/until summer solstice 
# for each measurement in each hemisphere

bc <- read_csv("./raw_data/temperature_BP_200609.csv") %>% 
  mutate(site = "BP") %>% 
  mutate(date_ibutton = as.Date(date_ibutton, format = "%m/%d/%Y")) %>% 
  # here, we calculate dss as the 'days since solstice', in N. America, June 21
  mutate(dss = difftime(date_ibutton, "2006-06-21", units = "days")) %>% 
  na.omit()

arg <- read_csv("./raw_data/temperature_PA_200603.csv") %>% 
  mutate(site = "PA") %>% 
  # and in the southern hemisphere, solstice is instead Sept. 21
  mutate(dss = difftime(date_ibutton, "2009-12-21", units = "days"))

# join data from the two locations before starting calculations and write data to csv
temperature_clean <- bc %>% 
  full_join(arg)
#write_csv(temperature_clean, "./clean_data/temperature_clean.csv")

# calculate the average daily maximum temperature

temp <- read_csv("./clean_data/temperature_clean.csv")

adm <- temp %>% 
  group_by(site, dss, number) %>% 
  summarize(dm = max(temp, na.rm = TRUE)) %>% 
  group_by(site, number) %>% 
  summarize(adm = mean(dm))

# we primarily care about modeling the average daily maximum between sites,
# not trends through time, so a simple linear model will do

lm.adm <- lm(adm ~ site, data = adm)

par(mfrow=c(2,2))

plot(lm.adm)

# inspecting the diagnostic plots, residuals look fairly good
# there may be a bit more variability in temperature at one of the sites
# some non-normality is present from the Q-Q plot, but at acceptable levels

summary(lm.adm)
anova(lm.adm)

# according to statistical tests, site explains 40% of the error in these data

# now for the upper 99th percentile of temperature

nnpercentile <- temp %>% 
  group_by(site, number) %>% 
  summarize(nnquant = quantile(temp , prob = 0.99))

lm.nnquant <- lm(nnquant ~ site, data = nnpercentile)

plot(lm.nnquant)
summary(lm.nnquant)
anova(lm.nnquant)

# Now for humidity data

# read in environmental data for BC and isolate only monthly humidity data
bc_hum <-  read_csv("./raw_data/humidity_BP.csv") %>% 
  group_by(year, month) %>% 
  summarize(humidity = mean(humidity, na.rm=TRUE)) %>% 
  na.omit() %>% 
  # only have data from Jan 2015-May 2020 for Nuevo Gulf, so filter out other BC data
  filter(year < 2020)

# read in data for Argentina
arg_hum <- read_csv("./raw_data/humidity_PA.csv") %>%
  filter(year >= 2015 & year < 2020)

# join together and write combined clean dataframe
hum <- bc_hum %>% 
  full_join(arg_hum)
#write_csv(hum, "./clean_data/humidity_clean.csv")

# model humidity between sites since we only really care about the overall mean conditions
lm.hum <- lm(humidity ~ site, data = hum)

plot(lm.hum)

# model meets assumptions - no glaring violations of heterogeneity or normality

summary(lm.hum)
anova(lm.hum)

# site explains 91% of variance in humidity data

# and finally, the code for creating figures

# first, code for the ribbon plot of average daily maximum temperature through time
temp_time_fig <- read_csv("./clean_data/temperature_clean.csv") %>% 
  # calculate daily maximum for each ibutton at each site
  group_by(site, dss, number) %>% 
  summarize(dm = max(temp)) %>% 
  # calculate the overall average daily maximum for each day at each site
  group_by(site, dss) %>% 
  summarize(adm = mean(dm), se_adm = se(dm))

temp_time <- ggplot(aes(x = dss, y = adm, colour = site, fill = site), data = temp_time_fig) +
  geom_ribbon(aes(ymin = adm - se_adm, ymax = adm + se_adm), alpha = 0.3, lwd = 0.7) +
  theme_classic() +
  theme(legend.position = c(0.5,0.15)) +
  scale_color_manual(values = c("steelblue3", "indianred3")) +
  scale_fill_manual(values = c("steelblue3", "indianred3")) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  labs(x = "Time until/since summer solstice (days)", y = "Daily maximum temperature (˚C)",
       colour = "Site", fill = "Site")

# save files to figure folder
#ggsave("./figures/Figure_1a.tiff", plot = temp_time, 
      #width = 11, height = 5, units = "in",
      #dpi = 600)

## Figures 1 b-d plotting theme

fig1bd_theme <- theme_classic() +
  theme(axis.text.x = element_text(size = 22)) +
  theme(axis.text.y = element_text(size = 22)) +
  theme(axis.title.x = element_text(size = 24)) +
  theme(axis.title.y = element_text(size = 24)) +
  theme(legend.text = element_text(size = 22)) +
  theme(legend.title = element_text(size = 24)) +
  theme(legend.position = "none")

# summarize the average daily maximum temperature average through time
adm_fig <- read_csv("./clean_data/temperature_clean.csv") %>% 
  # calculate maximum temp each day for each ibutton at each site
  group_by(site, dss, number) %>% 
  summarize(dm = max(temp)) %>% 
  # then average over time, leaving ibutton as the plotting replicate
  group_by(site, number) %>% 
  summarize(adm = mean(dm))

# also want to summarize the upper 99th quantile of temperature
temp_site_summary <- read_csv("./clean_data/temperature_clean.csv") %>% 
  # calculate the 99th quantile of temperature for each site with ibutton as replicate
  group_by(site, number) %>% 
  summarize(nnquant = quantile(temp , prob = 0.99)) %>% 
  full_join(adm_fig)

# making and saving figures for 1) average daily maximum temperature
adm_summary <- ggplot(aes(x = site, y = adm, colour = site), data = temp_site_summary) +
  geom_boxplot(width = 0.8, lwd = 1) +
  geom_jitter(width = 0.25) +
  scale_color_manual(values = c("steelblue3", "indianred3")) +
  labs(x = "Study site", y = "Mean daily max. temp. (˚C)", colour = "Site") +
  fig1bd_theme

#ggsave("./figures/Figure_1b.tiff", plot = adm_summary, 
      #width = 5, height = 6, units = "in",
      #dpi = 600)

# and a figure for upper 99th quantile of temperature
nnq_summary <- ggplot(aes(x = Site, y = nnquant, colour = site), data = temp_site_summary) +
  geom_boxplot(width = 0.8, lwd = 1) +
  geom_jitter(width = 0.25) +
  scale_color_manual(values = c("steelblue3", "indianred3")) +
  labs(x = "Study site", y = "Upper 99th quantile temp. (˚C)") +
  fig1bd_theme

#ggsave("./figures/Figure_1c.tiff", plot = nnq_summary, 
       #width = 5, height = 6, units = "in",
       #dpi = 600)

#summarize humidity data

hum <- read_csv("./clean_data/humidity_clean.csv")

hum_site <- ggplot(aes(x = site, y = humidity, color = site),
                   data = hum) +
  geom_boxplot(width = 0.8, lwd = 1) +
  geom_jitter(width = 0.25) +
  scale_color_manual(values = c("steelblue3", "indianred3")) +
  labs(x = "Study region", y = "Mean relative humidity (%)") +
  fig1bd_theme
hum_site

#ggsave("./figures/Figure_1d.tiff", plot = hum_site, 
      #width = 5, height = 6, units = "in",
      #dpi = 600)