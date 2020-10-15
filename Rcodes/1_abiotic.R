# 1: Analysis and visualization of temperature and humidity data


## First load packages of interest

packages <- c("tidyverse")

if (!require(packages, character.only = TRUE)) {
      install.packages(packages, dependencies = TRUE)
      library(packages, character.only = TRUE)
    }

# define standard error function for eventual plotting
se <- function(x){
  sd(x)/sqrt(length(x))
}

# First, read in dataframes and calculate the time since/until summer solstice 
# for each measurement in each hemisphere

bc <- read_csv("./raw_data/bc_ibutton.csv") %>% 
  mutate(site = "BP") %>% 
  mutate(date = as.Date(date, format = "%m/%d/%Y")) %>% 
  # here, we calculate dss as the 'days since solstice', in N. America, June 21
  mutate(dss = difftime(date, "2006-06-21", units = "days")) %>% 
  na.omit()

arg <- read_csv("./raw_data/arg_ibutton.csv") %>% 
  mutate(site = "PA") %>% 
  # and in the southern hemisphere, solstice is instead Sept. 21
  mutate(dss = difftime(date, "2009-12-21", units = "days"))

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
bc_hum <-  read_csv("./raw_data/bc_hum.csv") %>% 
  rename(year = "Year", month = "Month", day = "Day", temp = "Temp...C.",
         hum = "Rel.Hum....") %>% 
  group_by(year, month) %>% 
  summarize(humidity = mean(hum, na.rm=TRUE)) %>% 
  na.omit() %>% 
  mutate(site = "BS") %>% 
  # only have data from Jan 2015-May 2020 for Nuevo Gulf, so filter out other BC data
  mutate(remove = ifelse(year == 2020 & month > 5, TRUE, FALSE)) %>% 
  filter(remove == FALSE) %>% 
  select(-remove)

# read in data for Argentina
arg_hum <- read_csv("./raw_data/arg_hum.csv") %>% 
  select(1:3) %>% 
  mutate(site = "NG") %>% 
  filter(year >= 2015)

# join together and write combined clean dataframe
hum <- bc_hum %>% 
  full_join(arg_hum)
#write.csv(hum, "./clean_data/humidity_clean.csv")

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
  summarize(adm = mean(dm), se_adm = se(dm)) %>% 
  rename(Site = site)

temp_time <- ggplot(aes(x = dss, y = adm, colour = Site, fill = Site), data = temp_time_fig) +
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
  ylab("Daily maximum temperature (˚C)") +
  xlab("Time until/since summer solstice (days)")

ggsave("./figures/Figure_1a.tiff", plot = temp_time, 
       width = 11, height = 5, units = "in",
       dpi = 600)

# summarize the average daily maximum temperature average through time
adm_fig <- read_csv("./clean_data/temperature_clean.csv") %>% 
  # calculate maximum temp each day for each ibutton at each site
  group_by(site, dss, number) %>% 
  summarize(dm = max(temp)) %>% 
  # then average over time, leaving ibutton as the plotting replicate
  group_by(site, number) %>% 
  summarize(adm = mean(dm)) %>% 
  rename(Site = site)

# also want to summarize the upper 99th quantile of temperature
temp_site_summary <- read_csv("./clean_data/temperature_clean.csv") %>% 
  # calculate the 99th quantile of temperature for each site with ibutton as replicate
  group_by(site, number) %>% 
  summarize(nnquant = quantile(temp , prob = 0.99)) %>% 
  rename(Site = site) %>% 
  full_join(adm_fig)

# making and saving figures for 1) average daily maximum temperature
adm_summary <- ggplot(aes(x = Site, y = adm, colour = Site), data = temp_site_summary) +
  geom_boxplot(width = 0.8, lwd = 1) +
  theme_classic() +
  geom_jitter(width = 0.25) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("steelblue3", "indianred3")) +
  theme(axis.text.x = element_text(size = 22)) +
  theme(axis.text.y = element_text(size = 22)) +
  theme(axis.title.x = element_text(size = 24)) +
  theme(axis.title.y = element_text(size = 24)) +
  theme(legend.text = element_text(size = 22)) +
  theme(legend.title = element_text(size = 24)) +
  ylab("Mean daily max temp. (˚C)") +
  xlab("Study site")

ggsave("./figures/Figure_1b.tiff", plot = adm_summary, 
       width = 5, height = 6, units = "in",
       dpi = 600)

# and a figure for upper 99th quantile of temperature
nnq_summary <- ggplot(aes(x = Site, y = nnquant, colour = Site), data = temp_site_summary) +
  geom_boxplot(width = 0.8, lwd = 1) +
  theme_classic() +
  geom_jitter(width = 0.25) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("steelblue3", "indianred3")) +
  theme(axis.text.x = element_text(size = 22)) +
  theme(axis.text.y = element_text(size = 22)) +
  theme(axis.title.x = element_text(size = 24)) +
  theme(axis.title.y = element_text(size = 24)) +
  theme(legend.text = element_text(size = 22)) +
  theme(legend.title = element_text(size = 24)) +
  ylab("Upper 99th quantile temp. (˚C)") +
  xlab("Study site")

ggsave("./figures/Figure_1c.tiff", plot = nnq_summary, 
       width = 5, height = 6, units = "in",
       dpi = 600)

#summarize humidity data

hum <- read_csv("./clean_data/humidity_clean.csv") %>% 
  group_by(site, month) %>% 
  summarize(mean_hum = mean(humidity), se_hum = se(humidity)) %>% 
  mutate(month = as.factor(month))

hum_site <- ggplot(aes(x = site, y = mean_hum, color = site),
                   data = hum) +
  theme_classic() +
  geom_boxplot(width = 0.8, lwd = 1) +
  geom_jitter(width = 0.25) +
  ylab("Mean relative humidity (%)") +
  scale_color_manual(values = c("steelblue3", "indianred3")) +
  theme(axis.text.x = element_text(size = 22)) +
  theme(axis.text.y = element_text(size = 22)) +
  theme(axis.title.x = element_text(size = 24)) +
  theme(axis.title.y = element_text(size = 24)) +
  theme(legend.text = element_text(size = 22)) +
  theme(legend.title = element_text(size = 24)) +
  theme(legend.position = "none") +
  xlab("Study region")

ggsave("./figures/Figure_1d.tiff", plot = hum_site, 
       width = 5, height = 6, units = "in",
       dpi = 600)
