## April 2019
## Script to make Argentina & BC data analyzable by R
packages <- c("tidyverse")

for (i in 1:length(packages)){
  if (!require(packages[i], character.only = TRUE)) {
    install.packages(packages[i], dependencies = TRUE)
    library(packages[i], character.only = TRUE)
  }
}

# convert from wide format to long format
bc_gather <- read_csv("./raw_data/bc_bal_lot.csv") %>% 
  # algae need percent cover as measurement
  gather(key = "species", value = "percent_cover", 15:19) %>% 
  # and invertebrates need count data
  gather(key = "species_2", value = "count", 6:10) %>% 
  select(1:5, 24:27)

# need to merge species columns and keep the response columns
bc_tidy <- bc_gather %>% 
  # remove the count data
  select(-6,-7) %>% 
  # rename the column
  rename(species = species_2) %>%
  # omit empty data from timepoint = 0 for recruitment
  na.omit() %>% 
  # join back with original dataframe to get species in one column
  full_join(bc_gather) %>% 
  mutate(location = "BP") %>% 
  select(-species_2) %>% 
  mutate(date = format(as.Date(date, "%d-%b-%y"))) %>% 
  #calculate time between start of experiment and each survey date
  mutate(timediff = difftime(date,min(date), unit = c("weeks"))) %>%
  # and rename all the columns
  mutate(species = str_replace_all(species,
                  c("Bgrec" = "Balanus_recruits", "Cdrec" = "Chthamalus_recruits",
                    "#Lot" = "Lottia", "#Lot+added" = "Lottia_replaced",
                    "%U" = "Ulva", "%P"="Pyropia", "%s"="Urospora",
                    "%ephem"="Ephemerals","%peren"="Perennials")))

## repeat similar process for Argentina data

arg_gather <- read_csv("../data/arg_bal_lot.csv") %>% 
  mutate(location = "PA") %>% 
  gather(key = "species", value = "percent_cover", 10:14) %>% 
  gather(key = "species_2", value = "count", 7:10) %>% 
  select(2:6, 11:15)
  

# then take out the species that have percent cover (just counts)

arg_area_out <- arg_gathered %>% 
  select(-7,-8) %>% 
  rename("species" = "species_2")

arg_area_out[is.na(arg_area_out)] <- 0

arg_counts_out <- arg_gathered %>% 
  select(-9,-10)

arg_counts_out[is.na(arg_counts_out)] <- 0

# then merge them together to get all species in the same column

arg_all <- arg_area_out %>% 
  full_join(arg_counts_out)

arg_all <- unique(arg_all)

# rename those columns

arg_all <- arg_all %>% 
  rename("limpets" = "Siphonaria", "barnacles" = "Balanus",
         "replicate" = "REPLICATE", "date" = "MONTH")

# then figure out the dates and get in easy format for R to read ... yyyy-mm-dd

arg_date <- arg_all

arg_date$date <- gsub("Apr-06", "06-Apr", arg_date$date)
arg_date$date <- gsub("Jun-06", "06-Jun", arg_date$date)
arg_date$date <- gsub("Aug-06", "06-Aug", arg_date$date)
arg_date$date <- gsub("Dec-06", "06-Dec", arg_date$date)
arg_date$date <- gsub("Jan-07", "07-Jan", arg_date$date)

# put in dummy date of the middle of each month

arg_date <- arg_date %>% 
  mutate(day = 15) %>% 
  unite(date, day, date, sep = " ")

arg_date$date <- format(as.Date(arg_date$date, "%d %y-%b"))

# replace dummy dates with actual dates from Balanus limpets Argentina.xls sheet 'sampling dates'
arg_date$date <- gsub("2006-04-15", "2006-04-26", arg_date$date)
arg_date$date <- gsub("2006-05-15", "2006-05-12", arg_date$date)
arg_date$date <- gsub("2006-06-15", "2006-06-12", arg_date$date)
arg_date$date <- gsub("2006-07-15", "2006-07-11", arg_date$date)
arg_date$date <- gsub("2006-08-15", "2006-08-09", arg_date$date)
arg_date$date <- gsub("2006-09-15", "2006-09-09", arg_date$date)
arg_date$date <- gsub("2006-10-15", "2006-10-07", arg_date$date)
arg_date$date <- gsub("2006-11-15", "2006-11-09", arg_date$date)
arg_date$date <- gsub("2006-12-15", "2006-12-04", arg_date$date)
arg_date$date <- gsub("2007-01-15", "2007-01-06", arg_date$date)
arg_date$date <- gsub("2007-02-15", "2007-02-19", arg_date$date)

arg_all <- arg_date

arg_all$species <- gsub("# siphonaria recruits", "Siphonaria_recruits", arg_all$species)
arg_all$species <- gsub("#Balanus recruits", "Balanus_recruits", arg_all$species)
arg_all$species <- gsub("mussels", "Mussels", arg_all$species)
arg_all$species <- gsub("cover Ralfsia sp.", "Ralfsia", arg_all$species)
arg_all$species <- gsub("cover Blindigia", "Blindigia", arg_all$species)
arg_all$species <- gsub("cover Porphyra", "Pyropia", arg_all$species)
arg_all$species <- gsub("\\?", " ", arg_all$species)
arg_all$species <- gsub("polysiphonia", "Polysiphonia", arg_all$species)
arg_all$species <- gsub("Ephemeral Algae", "Ephemerals", arg_all$species)

arg_all$plot <- as.factor(arg_all$plot)

arg_all <- arg_all %>% 
  select(-replicate) %>% 
  mutate(block = "A")

# add in limpet densities to argentina data

siphonaria_density <- read_csv("../data/limpet_density.csv", col_names = FALSE)

siphonaria_density <- siphonaria_density %>% 
  slice(-1) %>% 
  rename("balanus" = "X1", "siphonaria" = "X2", "replicate" = "X3", "plot" = "X4",
         "date" = "X5", "density" = "X6")

siphonaria_1 <- siphonaria_density %>%
  select(1:4,7:8) %>% 
  rename("date" = "X7", "density" = "X8") %>% 
  full_join(siphonaria_density) %>% 
  select(-7,-8)

siphonaria_2 <- siphonaria_1 %>% 
  select(1:4,7:8) %>% 
  rename("date" = "X9", "density" = "X10") %>% 
  full_join(siphonaria_1) %>% 
  select(-7,-8) %>% 
  slice(35:length(balanus))
  
siphonaria_3 <- siphonaria_2 %>% 
  select(1:4,7:8) %>% 
  rename("date" = "X11", "density" = "X12") %>% 
  full_join(siphonaria_2) %>% 
  select(-7,-8) %>% 
  slice(69:length(balanus))

siphonaria_4 <- siphonaria_3 %>% 
  select(1:4,7:8) %>% 
  rename("date" = "X13", "density" = "X14") %>% 
  full_join(siphonaria_3) %>% 
  select(-7,-8) %>% 
  slice(103:length(balanus))

siphonaria_5 <- siphonaria_4 %>% 
  select(1:4,7:8) %>% 
  rename("date" = "X15", "density" = "X16") %>% 
  full_join(siphonaria_4) %>% 
  select(-7,-8) %>% 
  slice(137:length(balanus))

siphonaria_6 <- siphonaria_5 %>% 
  select(1:4,7:8) %>% 
  rename("date" = "X17", "density" = "X18") %>%
  full_join(siphonaria_5) %>% 
  select(-7,-8) %>% 
  slice(171:length(balanus))

siphonaria_7 <- siphonaria_6 %>% 
  select(1:4,7:8) %>% 
  rename("date" = "X19", "density" = "X20") %>% 
  full_join(siphonaria_6) %>% 
  select(-7,-8) %>% 
  slice(205:length(balanus))

siphonaria_8 <- siphonaria_7 %>% 
  select(1:4,7:8) %>% 
  rename("date" = "X21", "density" = "X22") %>%  
  full_join(siphonaria_7) %>% 
  select(-7,-8) %>% 
  slice(239:length(balanus))
  
siphonaria_9 <- siphonaria_8 %>% 
  select(1:4,7:8) %>% 
  rename("date" = "X23", "density" = "X24") %>%  
  full_join(siphonaria_8) %>% 
  select(-7,-8) %>% 
  slice(274:length(balanus))

siphonaria <- siphonaria_9 %>% 
  mutate(day = 1) %>% 
  unite(date, day, date, sep = " ")

siphonaria$date <- format(as.Date(siphonaria$date, "%d %y-%b"))
  
siphonaria <- siphonaria %>% 
  rename("count" = "density", "limpets" = "siphonaria",
         "barnacles" = "balanus") %>% 
  mutate(species = "Siphonaria") %>% 
  mutate(location = "PA") %>% 
  select(-replicate)

# reformat initial timepoint to be joined with rest of data

initial_timept <- read_csv("../data/arg_bal_lot_initial.csv")

initial <- initial_timept %>%
  select(-3, -4, -8, -9) %>% 
  mutate(location = "PA", block = "A") %>% 
  rename(barnacles = 'Balanus', limpets = 'Siphonaria') %>% 
  gather(key = "species", value = "count", 3) %>% 
  gather(key = "species_2", value = "percent_cover", 3:4)

## separate the count and percent cover species and do some renaming

initial_counts <- initial %>%
  select(-8,-9)

initial_cover <- initial %>% 
  select(-6,-7) %>% 
  rename("species" = "species_2")

initial_joint <- initial_counts %>% 
  full_join(initial_cover)

initial_joint$species <- gsub("Dens Siph", "Siphonaria", initial_joint$species)
initial_joint$species <- gsub("Cover Eph", "Ephemerals", initial_joint$species)
initial_joint$species <- gsub("Cover Per", "Perennials", initial_joint$species)

initial_complete <- initial_joint %>% 
  mutate(date = "2005-12-12")

initial_complete$limpets <- gsub("no","out", initial_complete$limpets)
initial_complete$limpets <- gsub("yes","in", initial_complete$limpets)

initial_complete$plot <- as.factor(initial_complete$plot)

# put together limpet densities and argentinian dataset with initial

siphonaria$count <- as.numeric(siphonaria$count)
siphonaria$plot <- as.factor(siphonaria$plot)

arg_II <- arg_all %>% 
  full_join(siphonaria) %>% 
  full_join(initial_complete) %>% 
  mutate(block = "A")

arg_II$plot <- as.factor(arg_II$plot)

arg_II$timediff <- difftime(arg_II$date, "2005-12-12", units = c("weeks"))

write_csv(arg_III, "argentina2.csv")


## unite all the data into one? Might not make sense for all variables
#, but could do for comparable things like ephem/peren/barnacle recuits, etc.

limpets_barnacles <- bc_all %>% 
  full_join(arg_II)

limpets_barnacles$date <- as.Date(limpets_barnacles$date)
limpets_barnacles$block <- as.factor(limpets_barnacles$block)
limpets_barnacles$barnacles <- as.factor(limpets_barnacles$barnacles)
limpets_barnacles$location <- as.factor(limpets_barnacles$location)
limpets_barnacles$plot <- as.factor(limpets_barnacles$plot)
limpets_barnacles$limpets <- as.factor(limpets_barnacles$limpets)

## make df for looking at effect of treatments on ephemerals
# note that in the df, ephemerals for BC is the sum of all ephemeral spp.

ephemerals <- limpets_barnacles %>% 
  filter(species == "Ephemerals") %>% 
  dplyr::select(-species, -count) %>% 
  filter(!c(timediff == 0 & location == "PA"))

write.csv(ephemerals, "../data/ephemerals.csv")

## for perennials

perennials <- limpets_barnacles %>% 
  filter(species == "Ralfsia" | species == "Perennials") %>% 
  dplyr::select(-species, -count)

write.csv(perennials, "../data/perennials.csv")

glimpse(limpets_barnacles)

## herbivores

herbivores <- limpets_barnacles %>% 
  filter(species == "Lottia" | species == "Siphonaria" | species == "Littorina") %>%
  mutate(remove = ifelse(species == "Lottia" & limpets != "con", TRUE, FALSE)) %>% 
  filter(remove == FALSE) %>% 
  mutate(remove = ifelse(species == "Siphonaria" & limpets != "con", TRUE, FALSE)) %>% 
  filter(remove == FALSE) %>% 
  dplyr::select(-percent_cover, -limpets, -remove)

write.csv(herbivores, "../data/herbivore_counts.csv")

# barnacle recruitment

balanus <- limpets_barnacles %>% 
  filter(species == "Balanus_recruits") %>% 
  dplyr::select(-species, -percent_cover, -date)

chthamalus <- limpets_barnacles %>% 
  filter(species == "Chthamalus_recruits") %>% 
  dplyr::select(-species, -percent_cover, -date, -location)

write.csv(balanus, "../data/balanus_recruitment.csv")
write.csv(chthamalus, "../data/chthamalus_recruitment.csv")

## read in fucus data too

fucus <- read_csv("../data/fucus.csv")

fucus$date <- format(as.Date(fucus$date, "%y-%m-%d"))
fucus$timediff <- difftime(fucus$date, "2006-06-12", units = c("weeks"))
fucus$limpets <- gsub('no','out', fucus$limpets)
fucus$balanus <- gsub('out','no', fucus$balanus)

fucus$timediff <- as.numeric(fucus$timediff)
fucus$limpets <- as.factor(fucus$limpets)
fucus$balanus <- as.factor(fucus$balanus)
fucus$block <- as.factor(fucus$block)
fucus$plot <- as.factor(fucus$plot)

fucus <- fucus %>% 
  rename(barnacles = balanus) %>% 
  rename(prop_cover = fucus_cover)
