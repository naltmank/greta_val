rm(list = ls())
# install.packages("librarian")
librarian::shelf(here, ggplot2, ggpubr, glmmTMB, DHARMa, car)

#### READ AND FORMAT DATA ####
disease <- read.csv(here::here("data", "disease_prevalence.csv"))
disease$healthy <- disease$total - disease $diseased

#### SIMPLE LOGISTIC REGRESSION ####
disease_mod <- glmmTMB(
  cbind(diseased, healthy) ~ position,
  family = binomial,
  data = disease
)

summary(disease_mod)
plot(simulateResiduals(disease_mod)) # massive residual issues
# treating categorical data as continuous
# likely spatial clustering in sites based on the bimodal residuals

Anova(disease_mod) # X2 = 24.8, P = 6.2E-07

#### PLOT ####
# extract model response
emm <- emmeans(
  disease_mod, ~ position,
  at = list(position = seq(0,17, 0.1)), # treating position as true continuouss
  type = "response"
)

emm_df <- as.data.frame(emm)

(disease_plot <-
  ggplot() +
    geom_point(data = disease, aes(x = position, y = diseased/total) ) +
    geom_ribbon(data = emm_df,
                aes( x = position, ymin = asymp.LCL, ymax = asymp.UCL),
                alpha = 0.2) +
    geom_line(data = emm_df, aes(x = position, y = prob )) +
  labs(y = "Proportion diseased", x = "Relative position") +
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

 ggsave(filename = here::here("output", "disease_plot.png"), disease_plot, width = 10, height = 8,
        dpi = "retina")
