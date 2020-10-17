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
  # remove the cover data
  select(-6,-7) %>% 
  # rename the column
  rename(species = species_2) %>% 
  # join back with original dataframe to get species in one column
  full_join(bc_gather) %>% 
  mutate(location = "BP") %>% 
  select(-species_2) %>% 
  mutate(date = format(as.Date(date, "%d-%b-%y")),
  #calculate time between start of experiment and each survey date
  timediff = difftime(date,min(date), unit = c("weeks")),
  # and rename all the columns
  species = str_replace_all(species,
                  c("Bgrec" = "Balanus_recruits", "Cdrec" = "Chthamalus_recruits",
                    "#Lot" = "Lottia", "#Lot+added" = "Lottia_replaced",
                    "#Lit" = "Littorina", "Ephemeral Algae" = "Ephemerals",
                    "%U" = "Ulva", "%P"="Pyropia", "%s"="Urospora",
                    "%ephem"="Ephemerals","%peren"="Perennials")),
  # convert all NA values to zeroes except for the cover column for countable species and vice versa
  percent_cover = if_else(is.na(percent_cover) & is.na(count), 0, percent_cover),
  count = if_else(is.na(count) & is.na(percent_cover), 0, count),
  date = as.Date(date, format = "%Y-%m-%d"))

# repeat same process for Argentina data
# NA values at this location are a mix of zeroes and missing data from plots getting dislodged
# remove data rows for lost plots, rest of NA are zeroes
arg_gather <- read_csv("./raw_data/arg_bal_lot.csv") %>%  
  gather(key = "species", value = "percent_cover", 10:14) %>% 
  gather(key = "species_2", value = "count", 7:10) %>% 
  filter(!grepl("lost", notes)) %>% 
  select(2:3, 5:6, 11:14)

arg_tidy <- arg_gather %>% 
  select(-5,-6) %>%
  rename("species" = "species_2") %>% 
  full_join(arg_gather) %>%
  select(-8) %>% 
  rename("limpets" = "Siphonaria", "barnacles" = "Balanus", "date" = "MONTH") %>% 
  mutate(date = as.character(date)) %>% 
  mutate(date = str_replace_all(date, c("06-Apr" = "2006-04-26", 
                         "06-May" = "2006-05-12",
                         "06-Jun" = "2006-06-12",
                         "06-Aug" = "2006-08-09",
                         "06-Sep" = "2006-09-09",
                         "06-Oct" = "2006-10-07",
                         "06-Nov" = "2006-11-09",
                         "06-Dec" = "2006-12-04",
                         "07-Jan" = "2007-01-06",
                         "07-Feb" = "2007-02-19")),
         # rename species to more comparable labels. Note Ralfsia is the only perennial here.
          species = str_replace_all(species, c("cover Ralfsia sp." = "Perennials",
                                    "# siphonaria recruits" = "Siphonaria_recruits",
                                    "cover Blindigia" = "Blindigia",
                                    "cover Porphyra" = "Porphyra",
                                    "polysiphonia?" = "Polysiphonia",
                                    "Ephemeral algae" = "Ephemerals",
                                    "#Balanus recruits" = "Balanus_recruits"
                                    )),
         # replace NA values with zeroes where they remain
  percent_cover = if_else(is.na(percent_cover) & is.na(count), 0, percent_cover),
  count = if_else(is.na(count) & is.na(percent_cover), 0, count),
  timediff = difftime(date, "2005-12-12", units = c("weeks"), tz = "PST"),
  date = as.Date(date, format = "%Y-%m-%d"))

# add in siphonaria densities to argentina data
siphonaria_density <- read_csv("./raw_data/limpet_density.csv", skip = 1,
         # name all density columns the dates that data were collected
          col_names = c("barnacles", "limpets", "replicate", "plot",
                       "x","2006-04-26", "x2", "2006-05-12",
                       "x3", "2006-06-12", "x4","2006-08-09",
                       "x5","2006-09-09", "x6", "2006-10-07",
                       "x7", "2006-11-09", "x8", "2006-12-04",
                       "x9", "2007-01-06", "x10", "2007-02-19")) %>% 
  # take all columns with density information and put into one
  pivot_longer(cols = starts_with("200"),
               names_to = "date",
               values_to = "count") %>% 
  select(1:2, 4, 15, 16) %>% 
  mutate(species = "Siphonaria", date = as.Date(date, format = "%Y-%m-%d")) %>% 
  mutate(timediff = difftime(date, "2005-12-12", units = c("weeks"), tz = "PST"))
# note that unlike previous dataframes, missing values here are lost plots and so are real NAs
  

# tidy initial timepoint data to be joined with rest of data
initial_timept <- read_csv("./raw_data/arg_bal_lot_initial.csv") %>% 
  select(-3, -4, -8, -9) %>% 
  rename(barnacles = 'Balanus', limpets = 'Siphonaria') %>%
  # get cover and count into their own columns w associated species column
  pivot_longer(cols = 4:5, names_to = "species", values_to = "percent_cover") %>% 
  pivot_longer(cols = 3, names_to = "species_2", values_to = "count")

initial_tidy <- initial_timept %>% 
  # remove cover data
  select(-4,-5) %>% 
  # to rename species column
  rename(species = "species_2") %>% 
  # and then rejoin and remove species_2 column
  full_join(initial_timept) %>% 
  select(-7) %>% 
  mutate(species = str_replace_all(species, c("Dens Siph" = "Siphonaria",
                                              "Cover Eph" = "Ephemerals",
                                              "Cover Per" = "Perennials")),
         limpets = str_replace_all(limpets, c("no" = "out", "yes" = "in")),
         date = as.Date("2005-12-12"), 
         timediff = difftime(date, "2005-12-12", units = c("weeks"), tz = "PST"))


# put together all of the argentina datasets
argentina_tidy <- arg_tidy %>% 
  full_join(siphonaria_density) %>% 
  full_join(initial_tidy) %>% 
  mutate(location = "PA", block = "A") %>% 
  unique()

# unite the data from both sites into one large df

bio_responses <- bc_tidy %>% 
  full_join(argentina_tidy)

#write_csv(bio_responses, "./clean_data/bio_responses.csv")

###

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

#write.csv(herbivores, "../data/herbivore_counts.csv")

#barnacle recruitment

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
