rm(list = ls())
# install.packages("librarian")
librarian::shelf(here, stringr, tidyverse, ggplot2, ggpubr, vegan, glmmTMB,
                 dplyr, janitor, DHARMa, car, emmeans, performance)

#### FUNCTIONS AND FORMATTING ####
# smithson verkuilen transformation for 0-1 continuous data
sv_trans <- function(prop, s = 0.000005){
  (prop*(length(prop) - 1) + s)/length(prop)
  # where prop is the vector of the proportional value you're transforming, 
  # N is the sample size, which is specified by taking the number of rows/observations from a given dataframe,
  # and s is a small offset 
}


#### READ AND FORMAT DATA ####
benthic <- read.csv(here("data", "benthic_cover.csv"), stringsAsFactors = T)

# convert benthic data to proportional data instead of the number of points
benthic_prop <- benthic %>%
  # filter empty rows
  filter(substrate != "") %>%
  group_by(site, year, transect) %>%
  mutate(
        # for model
        year = factor(year, levels = c(2024, 2022, 2019)),
        
        # to model as cover instead of points
        prop_cover = n_points / sum(n_points, na.rm = TRUE),
         
         # transform for modeling
         cover_trans = sv_trans(prop_cover),
         
         # for random effect
         transect_id = paste0(site, "-", transect)) %>%
  ungroup()

coral <- subset(benthic_prop, substrate == "hard coral")
soft_coral <- subset(benthic_prop, substrate == "soft coral")
macroalgae <- subset(benthic_prop, substrate == "macroalgae")
sponge <- subset(benthic_prop, substrate == "sponge")
other <- subset(benthic_prop, substrate == "other")

#### UNIVARIATE MODELS ####
# coral
coral_mod <- glmmTMB(cover_trans ~ year + (1|site/transect),
                     dispformula = ~year, # model residuals on year to handel overdispersion
                     family = beta_family(), data = coral)
plot(simulateResiduals(coral_mod)) # no more issues
summary(coral_mod)
Anova(coral_mod) # X2 = 0.8, P = 0.66

# extract estimated marginal means
coral_emm <- as.data.frame( emmeans(coral_mod, ~ year, type = "response")) %>%
  mutate(substrate = rep("Hard coral", n()))

# soft coral
soft_coral_mod <- glmmTMB(cover_trans ~ year + (1|site/transect),
                     dispformula = ~year,
                     family = beta_family(), data = soft_coral)
plot(simulateResiduals(soft_coral_mod)) # positive skew in residuals
summary(soft_coral_mod)
Anova(soft_coral_mod) # X2 = 3.7, P = 0.15

# extract estimated marginal means
soft_coral_emm <- as.data.frame( emmeans(soft_coral_mod, ~ year, type = "response")) %>%
  mutate(substrate = rep("Soft coral", n()))

# macroalgae
macroalgae_mod <- glmmTMB(cover_trans ~ year + (1|site/transect),
                     family = beta_family(), data = macroalgae)
plot(simulateResiduals(macroalgae_mod)) # no  issues
summary(macroalgae_mod)
Anova(macroalgae_mod) # X2 = 24.8, P < 0.001

macroalgae_emm <- as.data.frame( emmeans(macroalgae_mod, ~ year, type = "response")) %>%
  mutate(substrate = rep("Macroalgae", n())
  )

# sponge
sponge_mod <- glmmTMB(cover_trans ~ year + (1|site/transect),
                          family = beta_family(), data = sponge)
plot(simulateResiduals(sponge_mod)) # no  issues
summary(sponge_mod)
Anova(sponge_mod) # X2 = 10.8, P = 0.004

sponge_emm <- as.data.frame( emmeans(sponge_mod, ~ year, type = "response")) %>%
  mutate(substrate = rep("Sponges", n())
  )

##### PLOT #####
combined_emm <- rbind(coral_emm, soft_coral_emm, macroalgae_emm, sponge_emm) 

(univariate_plot <-
  ggplot() +
    geom_point(data = combined_emm, aes(x = response, y = substrate, colour = year),
               position = position_dodge(0.8), size = 4) +
    geom_errorbar(data = combined_emm, aes(xmin = response - asymp.LCL, xmax = response + asymp.UCL,
                                           y = substrate,
                                           colour = year),
                  position = position_dodge(0.8), width = 0.1) +
    scale_y_discrete(limits = c("Hard coral", "Soft coral", "Sponges", "Macroalgae")) +
    labs(x = "Estimated marginal mean", y = "") +
    scale_colour_manual(name = "Year",
                        limits = c("2019", "2022", "2024"),
                        values = c("2024" = "#D55E00",
                                   "2022" = "purple3",
                                   "2019" = "yellow4")) +
    theme_classic() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.title = element_text(colour = "black", size = 24),
          legend.text = element_text(colour = "black", size = 20),
          legend.position = "bottom",
          axis.line = element_line(colour = "black", size = 1.2),
          plot.title = element_text(color = "black", size = 25, hjust = 0, vjust = 0, face = "plain"),
          axis.text.x = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain"),
          axis.title.x = element_text(color = "black", size = 25, hjust = .5, vjust = 0, face = "plain"),
          axis.text.y = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain"),
          axis.title.y = element_text(color = "black", size = 25, hjust = .5, vjust = 0, face = "plain", margin = margin(r = 15))
    )
)

# ggsave(filename = here::here("output", "benthic_cover_effects.png"), univariate_plot, width = 12, height = 8,
#        dpi = "retina")
