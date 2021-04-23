## April 2019
## Script to make Argentina & BC data analyzable by R
pkgs <- c("tidyverse", "lubridate")
lapply(pkgs, install.packages, character.only = TRUE)
lapply(pkgs, library, character.only = TRUE)

# convert from wide format to long format
bc_cover <- read_csv("./raw_data/community_BP_200707.csv") %>% 
  select(1:5, 15:19) %>% 
  # algae need percent cover as measurement
  pivot_longer(cols = 6:10, names_to = "species", values_to = "percent_cover") %>% 
  mutate(percent_cover = if_else(is.na(percent_cover),0, percent_cover))
  
bc_tidy <- read_csv("./raw_data/community_BP_200707.csv") %>% 
  select(1:10) %>% 
  # and invertebrates need count data
  pivot_longer(cols = 6:10, names_to = "species", values_to = "count") %>% 
  mutate(count = if_else(is.na(count),0, count)) %>% 
  # join back with original dataframe to get species in one column
  full_join(bc_cover) %>% 
  mutate(location = "BP", date = ymd(date),
  #calculate time between start of experiment and each survey date
  timediff = difftime(date,min(date), unit = c("weeks")),
  # and rename all the columns
  species = str_replace_all(species,
                  c("Bgrec" = "Balanus_recruits", "Cdrec" = "Chthamalus_recruits",
                    "Lot" = "Lottia", "Lot_plus_added" = "Lottia_replaced",
                    "Lit" = "Littorina",
                    "%U" = "Ulva", "%P"="Pyropia", "%s"="Scunge",
                    "%ephem"="Ephemerals","%peren"="Perennials")))

# repeat same process for Argentina data
# NA values at this location are a mix of zeroes and missing data from plots getting dislodged
# remove data rows for lost plots, rest of NA are zeroes
arg_cover <- read_csv("./raw_data/community_PA_200602.csv") %>% 
  select(1:4, 8:12, 14) %>% 
  pivot_longer(names_to = "species", values_to = "percent_cover", cols = 5:9) %>% 
  filter(!grepl("lost", notes)) %>% 
  mutate(percent_cover = if_else(is.na(percent_cover), 0, percent_cover))

arg_tidy <- read_csv("./raw_data/community_PA_200602.csv") %>% 
  select(1:7,13:14) %>% 
  pivot_longer(names_to = "species", values_to = "count", cols = 5:8) %>% 
  filter(!grepl("lost", notes)) %>% 
  mutate(count = if_else(is.na(count),0, count)) %>% 
  full_join(arg_cover) %>%
  mutate(date = ymd(date), 
         # rename species to more comparable labels. Note Ralfsia is the only perennial here.
          species = str_replace_all(species, c("cover Ralfsia sp." = "Perennials",
                                    "Dens Siph" = "Siphonaria_recruits",
                                    "cover Blindigia" = "Blindigia",
                                    "cover Porphyra" = "Porphyra",
                                    "polysiphonia" = "Polysiphonia",
                                    "%ephem" = "Ephemerals",
                                    "Bgrec" = "Balanus_recruits",
                                    "mussels" = "Mussels"
                                    )),
         # replace NA values with zeroes where they remain
  timediff = difftime(date, "2005-12-12", units = c("weeks"), tz = "PST"),
  date = as.Date(date, format = "%Y-%m-%d"))


# add in siphonaria densities to argentina data
siphonaria_density <- read_csv("./raw_data/siphonaria_PA_200702.csv", skip = 1,
         # name all density columns the dates that data were collected
          col_names = c("barnacles", "herbivores", "plot",
                       "x","2006-04-26", "x2", "2006-05-12",
                       "x3", "2006-06-12", "x4","2006-08-09",
                       "x5","2006-09-09", "x6", "2006-10-07",
                       "x7", "2006-11-09", "x8", "2006-12-04",
                       "x9", "2007-01-06", "x10", "2007-02-19")) %>% 
  # take all columns with density information and put into one
  pivot_longer(cols = starts_with("200"),
               names_to = "date",
               values_to = "count") %>% 
  select(1:3, 14:15) %>% 
  mutate(species = "Siphonaria", date = ymd(date)) %>% 
  mutate(timediff = difftime(date, "2005-12-12", units = c("weeks"), tz = "PST"))
# note that unlike previous dataframes, missing values here are lost plots and so are real NAs

# tidy initial timepoint data to be joined with rest of data
initial_timept <- read_csv("./raw_data/community_PA_200512.csv") %>% 
  # get cover and count into their own columns w associated species column
  pivot_longer(cols = 4:5, names_to = "species", values_to = "percent_cover") %>% 
  pivot_longer(cols = 3, names_to = "species_2", values_to = "count") %>% 
  mutate(percent_cover = if_else(is.na(percent_cover), 0, percent_cover),
         count = if_else(is.na(count),0, count))

initial_cover <- initial_timept %>% 
  select(1:5)

initial_counts <- initial_timept %>% 
  # remove cover data
  select(-4,-5) %>% 
  # to rename species column
  rename(species = "species_2") %>% 
  # and then rejoin and remove species_2 column
  full_join(initial_cover) %>% 
  mutate(species = str_replace_all(species, c("Dens Siph" = "Siphonaria",
                                              "%ephem" = "Ephemerals",
                                              "%peren" = "Perennials")),
         date = as.Date("2005-12-12"), 
         timediff = difftime(date, "2005-12-12", units = c("weeks"), tz = "PST")) %>% 
  unique()

# put together all of the argentina datasets
argentina_tidy <- arg_tidy %>% 
  full_join(siphonaria_density) %>% 
  full_join(initial_tidy) %>% 
  mutate(location = "PA", block = "A") %>% 
  select(-notes)

# unite the data from both sites into one large df
bio_responses <- bc_tidy %>% 
  full_join(argentina_tidy) %>% 
  unique()

#write_csv(bio_responses, "./clean_data/bio_responses.csv")

## fucus cover data

fucus <- read_csv("./raw_data/fucus_BP_200707.csv") %>% 
  mutate(date = ymd(date),
         timediff = difftime(date, min(date), units = c("weeks")),
         percent_cover = fucus_cover*100)

#write_csv(fucus, "./clean_data/fucus_clean.csv")
