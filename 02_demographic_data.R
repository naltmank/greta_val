rm(list = ls())
# install.packages("librarian")
librarian::shelf(here, stringr, tidyverse, ggplot2, ggpubr, vegan, glmmTMB,
                 dplyr, janitor, DHARMa, car, emmeans, performance)


#### READ AND FORMAT DATA ####
demo_raw <- read.csv(here("data", "coral_demographics.csv"), stringsAsFactors = T)
naming_key <- read.csv(here("data", "naming_key.csv"), stringsAsFactors = T)

demo <- demo_raw %>%
  mutate(
    # for model
    year = factor(year, levels = c(2024, 2022, 2019)),

    # for random effect
    transect_id = paste0(site, "-", transect)) %>%
  ungroup()

demo_newname <- demo %>%
  left_join(naming_key, by = "species") %>%
  select(-species) %>%
  rename(species = corrected_name)

# for total density
total_density <- demo_newname %>%
  group_by(site, year, transect, transect_id) %>%
  summarise(total_count = sum(count), .groups = "drop")

species_obs <- demo_newname %>%
  group_by(species) %>%
  summarise(nonzero = sum(count > 0),
            total = sum(count))

#### UNIVARIATE MODELS ####
##### TOTAL DENSITY #####
total_density_mod <- glmmTMB(total_count ~ year + (1|site/transect_id),
                 family = poisson(),
                 data = total_density)
plot(simulateResiduals(total_density_mod)) # no issues
summary(total_density_mod)
Anova(total_density_mod) # X2 = 1815.2, P < 0.0001

emmeans(total_density_mod, pairwise ~ year)
total_emm <- as.data.frame( emmeans(total_density_mod, ~ year, type = "response")) %>%
  mutate(response = rep("Total density", n()))

##### PLOT #####
(total_density_plot <- ggplot() +
   geom_boxplot(data = total_density, aes(x = year, y = total_count, colour = year), outlier.shape = NA) +
   geom_jitter(data = total_density, aes(x = year, y = total_count, colour = year)) +
   scale_colour_manual(name = "Year",
                       limits = c("2019", "2022", "2024"),
                       values = c("2024" = "#D55E00",
                                  "2022" = "purple3",
                                  "2019" = "yellow4")) +
   scale_x_discrete(limits = c("2019", "2022", "2024")) + 
   scale_x_discrete(limits = c("2019", "2022", "2024")) + 
   labs(y = "Number of colonies", x = "") +
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

# ggsave(filename = here::here("output", "total_density.png"), total_density_plot, width = 10, height = 8,
#        dpi = "retina")

#### SPECIES INTERACTION MODEL ####
interaction_mod <- glmmTMB(count ~ year * species + (1|site) + (1 | transect_id),
                             ziformula = ~ 1, 
                             family = nbinom2(),
                             data = demo_newname)
plot(simulateResiduals(interaction_mod)) # dispersion test sig but also too much power - looks fine
summary(interaction_mod)
Anova(interaction_mod)

emmeans(interaction_mod, pairwise ~ year | species, type = "response")
interaction_emm <- as.data.frame( emmeans(interaction_mod, ~ year | species, type = "response"))  %>%
  # low values makes some CIs not estimate properly, resulting in Infinite values. Replace with NAs
  mutate(
    asymp.LCL = ifelse(is.infinite(asymp.LCL), NA, asymp.LCL),
    asymp.UCL = ifelse(is.infinite(asymp.UCL), NA, asymp.UCL)
  )

##### PLOT #####
# subset disease susceptibility info
susceptibility <- demo_newname %>%
  select(species, SCTLD_suscep) %>%
  group_by(species, SCTLD_suscep) %>%
  summarise(species = unique(species),
            SCTLD_suscep = unique(SCTLD_suscep))

# merge with emmeans object
interaction_plot <- interaction_emm %>%
  left_join(susceptibility, by = "species")

# subset by disease susceptibility
l_suscep_plot_df <- subset(interaction_plot, SCTLD_suscep == "low") 
m_suscep_plot_df <- subset(interaction_plot, SCTLD_suscep == "med")
h_suscep_plot_df <- subset(interaction_plot, SCTLD_suscep == "high") 
v_suscep_plot_df <- subset(interaction_plot, SCTLD_suscep == "very high") 

# create plots for each

(l_suscep_plot <-
    ggplot() +
    geom_point(data = l_suscep_plot_df, aes(x = response, y = species, colour = year),
               position = position_dodge(0.8), size = 4) +
    geom_errorbar(data = l_suscep_plot_df, aes(xmin = response - asymp.LCL, xmax = response + asymp.UCL,
                                           y = species,
                                           colour = year),
                  position = position_dodge(0.8), width = 0.5, linewidth = 1.2) +
    labs(x = "Estimated marginal mean", y = "", title = "a. Low susceptibility") +
    scale_colour_manual(name = "Year",
                        limits = c("2019", "2022", "2024"),
                        values = c("2024" = "#D55E00",
                                   "2022" = "purple3",
                                   "2019" = "yellow4")) +
#    scale_x_discrete(limits = c("2019", "2022", "2024")) + 
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
          axis.text.y = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "italic"),
          axis.title.y = element_text(color = "black", size = 25, hjust = .5, vjust = 0, face = "plain", margin = margin(r = 15))
    )
)

(m_suscep_plot <-
    ggplot() +
    geom_point(data = m_suscep_plot_df, aes(x = response, y = species, colour = year),
               position = position_dodge(0.8), size = 4) +
    geom_errorbar(data = m_suscep_plot_df, aes(xmin = response - asymp.LCL, xmax = response + asymp.UCL,
                                               y = species,
                                               colour = year),
                  position = position_dodge(0.8), width = 0.5, linewidth = 1.2) +
    labs(x = "Estimated marginal mean", y = "", title = "b. Medium susceptibility") +
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
          axis.text.y = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "italic"),
          axis.title.y = element_text(color = "black", size = 25, hjust = .5, vjust = 0, face = "plain", margin = margin(r = 15))
    )
)


(h_suscep_plot <-
    ggplot() +
    geom_point(data = h_suscep_plot_df, aes(x = response, y = species, colour = year),
               position = position_dodge(0.8), size = 4) +
    geom_errorbar(data = h_suscep_plot_df, aes(xmin = response - asymp.LCL, xmax = response + asymp.UCL,
                                               y = species,
                                               colour = year),
                  position = position_dodge(0.8), width = 0.5, linewidth = 1.2) +
    labs(x = "Estimated marginal mean", y = "", title = "c. High susceptibility") +
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
          axis.text.y = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "italic"),
          axis.title.y = element_text(color = "black", size = 25, hjust = .5, vjust = 0, face = "plain", margin = margin(r = 15))
    )
)


(v_suscep_plot <-
    ggplot() +
    geom_point(data = v_suscep_plot_df, aes(x = response, y = species, colour = year),
               position = position_dodge(0.8), size = 4) +
    geom_errorbar(data = v_suscep_plot_df, aes(xmin = response - asymp.LCL, xmax = response + asymp.UCL,
                                               y = species,
                                               colour = year),
                  position = position_dodge(0.8), width = 0.5, linewidth = 1.2) +
    labs(x = "Estimated marginal mean", y = "", title = "d. Very high susceptibility") +
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
          axis.text.y = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "italic"),
          axis.title.y = element_text(color = "black", size = 25, hjust = .5, vjust = 0, face = "plain", margin = margin(r = 15))
    )
)

(demo_panel <- ggarrange(l_suscep_plot, m_suscep_plot, h_suscep_plot, v_suscep_plot,
          ncol = 2, nrow = 2, common.legend = T, legend = "bottom"))

# ggsave(filename = here::here("output", "demo_effects_panel.png"), demo_panel, width = 22, height = 16,
#        dpi = "retina")

#### BETA DISPERSION ####
# create unique sample ID
demo2 <- demo %>%
  mutate(sample_id = paste0(site, "-", transect_id, "-", year))

# create species matrix
demo_wide <- demo2 %>%
  select(sample_id, year, site, transect_id, species, count) %>%
  pivot_wider(
    names_from = species,
    values_from = count,
    values_fill = 0
  )

# save metadata
meta <- demo_wide %>%
  select(sample_id, year, site, transect_id,)



# abundance matrix only
comm_df <- demo_wide[,-which(names(meta) %in% names(meta))]

# square root transform
comm_mat <- sqrt(as.matrix(comm_df))

# Bray-Curtis dissimilarity
demo_bray <- vegdist(comm_mat, method = "bray")

# beta dispersion
demo_disp <- betadisper(demo_bray, group = meta$year)

# merge distances back to metadata
disp_df <- meta %>%
  mutate(distance = demo_disp$distances)

# model
disp_mod <- glmmTMB(distance ~ year + (1|site/transect_id),
                    family = Gamma("log"),
                    data = disp_df)
plot(simulateResiduals(disp_mod)) # no issues
summary(disp_mod)
Anova(disp_mod) # X2 = 23.003, P = 1.011e-05
emmeans(disp_mod, pairwise ~ year)

(dispersion_plot <- ggplot() +
    geom_boxplot(data = disp_df, aes(x = year, y = distance, colour = year), outlier.shape = NA) +
    geom_jitter(data = disp_df, aes(x = year, y = distance, colour = year)) +
    scale_colour_manual(name = "Year",
                        limits = c("2019", "2022", "2024"),
                        values = c("2024" = "#D55E00",
                                   "2022" = "purple3",
                                   "2019" = "yellow4")) +
    scale_x_discrete(limits = c("2019", "2022", "2024")) + 
    labs(y = "Distance to centroid", x = "") +
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

# ggsave(filename = here::here("output", "beta_dispersion.png"), dispersion_plot, width = 10, height = 8,
#        dpi = "retina")
