rm(list = ls())
# install.packages("librarian")
librarian::shelf(plyr, dplyr, tidyverse, reshape2, here)


#### READ AND FORMAT DATA ####
demo <- read.csv(here::here("data/raw", "demo.csv")  )

# replace NA values with 0s
demo <- demo %>% mutate(count = replace_na(count, 0))

# remove total colony counts from dataset - calculate those manually later
demo <- demo %>% subset(species != "tot colonies")


#### CLEAN DATES AND SITE NAMES ####
demo <- demo %>%
  mutate(site = str_to_lower(str_squish(site)),        # removes extra whitespace
         site = str_replace_all(site, " ", ""),        # removes spaces
         site = str_remove_all(site, "[^a-z0-9_]"))  #removes extra special characters

demo$site <- dplyr::recode(demo$site, "cbclagoon" = "lagoon",
                    "cbclagoonreef" = "lagoon",
                    "cbclagoonpatch" = "lagoon",
                    "southwatersw" = "southwater",
                    "sw" = "southwater",
                    "srdisease" = "sr_disease")


levels(as.factor(demo$site))

levels(as.factor(demo$species))

sums <- demo %>% group_by(species) %>%
  summarize(total = sum(count)) %>% arrange(total)

demo$species <- dplyr::recode(demo$species, "Ag humilis" = "Ag. agaricites",
                       "O. annularis" = "Orbicella spp.",
                       "Orbicella" = "Orbicella spp.",
                       "O. franksi" = "Orbicella spp.",
                       "O. faveolata" = "Orbicella spp.",
                       "P. furcata" = "P. porites",
                       "Ag. tenufolia" = "Ag. tenuifolia")

demo <- demo %>% mutate("SCTLD_suscep" = species)

demo$SCTLD_suscep <- dplyr::recode(demo$SCTLD_suscep, 
                            "Ag fragilis" = "low",
                            "Ag. agaricites" = "low",
                            "Ag. tenuifolia" = "low",
                            "Colpo. natans" = "high",
                            "Dichocoenia " = "very high",
                            "Dipl. clivosa" = "high",
                            "Dipl. labyrinth." = "high",
                            "Dipl. strigosa" = "high",
                            "Eusmilia fast." = "very high",
                            "Favia fragum" = "low",
                            "Helio cuc." = "low",
                            "Isophyllia" = "low",
                            "M. cavernosa" = "med",
                            "Madracis decac." = "low",
                            "Madracis mir." = "low",
                            "Meandrina mea." = "very high",
                            "Mycetophyllia" = "low",
                            "Orbicella spp." = "med",
                            "P. astreoides" = "low",
                            "P. porites" = "low",
                            "Sid. radians" = "low",
                            "Sid. siderea" = "med",
                            "Stephanocoenia" = "med")

demo <- demo %>% group_by(site, year, transect, species, SCTLD_suscep) %>%
  summarize(count = sum(count))

# write.csv(demo, here::here("data", "coral_demographics.csv"), row.names = F  )