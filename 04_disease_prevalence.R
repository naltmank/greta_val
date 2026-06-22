rm(list = ls())
# install.packages("librarian")
librarian::shelf(here, ggplot2, ggpubr, glmmTMB, DHARMa, car,
                 mgcv, MuMIn, mgcViz, sf, spdep)

#### READ AND FORMAT DATA ####
disease <- read.csv(here::here("data", "disease_prevalence.csv"))
disease$healthy <- disease$total - disease$diseased
disease$prevalence <- disease$diseased / disease$total

#### SIMPLE LOGISTIC REGRESSION ####
disease_mod <- glmmTMB(
  cbind(diseased, healthy) ~ latitude,
  family = binomial,
  data = disease
)

summary(disease_mod)
plot(simulateResiduals(disease_mod)) # massive residual issues


Anova(disease_mod) # X2 = 20.8, P = 6.2E-07

#### PLOT ####
# extract model response
emm <- emmeans(
  disease_mod, ~ latitude,
  at = list(latitude = seq(16.7,17.1, 0.01)), # treating position as true continuouss
  type = "response"
)

emm_df <- as.data.frame(emm)

(disease_plot <-
  ggplot() +
    geom_point(data = disease, aes(x = latitude, y = prevalence, size = total) ) +
    geom_ribbon(data = emm_df,
                aes( x = latitude, ymin = asymp.LCL, ymax = asymp.UCL),
                alpha = 0.2) +
    geom_line(data = emm_df, aes(x = latitude, y = prob )) +
  labs(y = "Proportion diseased", x = "Latitude",
       size = "No. individuals surveyed") +
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

# ggsave(filename = here::here("output", "disease_plot_glm.png"), disease_plot, width = 10, height = 8,
#        dpi = "retina")

#### GAM ####

gam_lat <- gam(
  cbind(diseased, healthy) ~ s(latitude, k = 5),
  family = binomial,
  data = disease,
  method = "REML"
)

summary(gam_lat) 
gam.check(gam_lat) # k' = 4, edf = 3.95, k-index = 1.63, p = 1
plot(simulateResiduals(gam_lat)) 
plot(gam_lat)
# non heterogeneous residuals, but could be due to lack of power

gam_spat <- gam(
  cbind(diseased, healthy) ~ s(longitude, latitude, k = 5),
  family = binomial,
  data = disease,
  method = "REML"
)

summary(gam_spat)
gam.check(gam_spat)

AICc(gam_lat, gam_spat, disease_mod) # latitude model best

# try to determine if residual problems are due to spatial autocorrelation
disease$resid <- residuals(
  gam_lat,
  type = "deviance"
)

# quick and dirty plot
ggplot(data = disease,
       # plot lat and lon
       aes(x = longitude, y = latitude,
           # color based on deviation from fittend trendline
           color = resid)) +
  geom_point(size = 4) +
  coord_equal() +
  scale_color_gradient2() +
  theme_bw() # no clear trends to naked eye

# try Moran's I for more quantitative
coords <- cbind(
  disease$longitude,
  disease$latitude
)

nb <- knearneigh(coords, k = 3)
nb <- knn2nb(nb)
lw <- nb2listw(nb)
moran.test(disease$resid, lw) # P = 0.89 - no real autocorrelation

#### GAM PLOT PREP ####
# create df of possible latitude with length = 200 for smooth appearance
lat_df <- data.frame(
  latitude = seq(
    min(disease$latitude),
    max(disease$latitude),
    length.out = 200
  )
)

# create prediction df
pred <- predict(
  gam_lat,
  newdata = lat_df,
  type = "link",
  se.fit = TRUE
)

# add model fit to lat_df and 95% CI based on predicted SEs
lat_df_final <- lat_df %>%
  mutate(
    fit = plogis(pred$fit),
    lwr = plogis(pred$fit - 1.96 * pred$se.fit),
    upr = plogis(pred$fit + 1.96 * pred$se.fit)
  )

#### GAM PLOT ####
(gam_plot <- 
   ggplot() +
      # add raw data
      geom_point(
        data = disease,
        aes(x = latitude, y = prevalence, size = total)) +
      # add trendline
      geom_line(data = lat_df_final, aes(x = latitude, y = fit)) +
      # add CI
      geom_ribbon(
        data = lat_df_final,
        aes(x = latitude, ymin = lwr, ymax = upr), alpha = 0.25) +
      labs(
        x = "Latitude",
        y = "Proportion diseased",
        size = "No. individuals surveyed"
      ) +
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

# ggsave(filename = here::here("output", "disease_plot_gam.png"), gam_plot, width = 10, height = 8,
#        dpi = "retina")
