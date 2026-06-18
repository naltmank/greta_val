rm(list = ls())
# install.packages("librarian")
librarian::shelf(nlme, plyr, car, tidyverse, reshape2, vegan, ggrepel, viridis, kableExtra,
                 patchwork, ggnewscale, here)

#### SET UP - AES AND FUNCTIONS ####
pd <- position_dodge(width = 0.93)
se<-function(x)sqrt(var(x)/length(x))

#### READ AND PREP DATA ####
demo <- read.csv(here::here("data", "coral_demographics.csv"), stringsAsFactors = T)
naming_key <- read.csv(here("data", "naming_key.csv"), stringsAsFactors = T)

demo <- demo %>%
  unite("event", c("site", "year", "transect"), sep = "_", remove = FALSE) 

demo_newname <- demo %>%
  left_join(naming_key, by = "species") %>%
  select(-species) %>%
  rename(species = corrected_name)

cast <- demo_newname %>% dplyr::select(-SCTLD_suscep) %>%
  pivot_wider(id_cols = c(site, year, transect, event), names_from = species, 
              values_from = count)


cast$year <- as.factor(cast$year)

#### PREP MULTIVARIATE STATS ####

# create bray-curtis dissimilarity matrix
mat = cast[,4:ncol(cast)]
mat <- mat %>% remove_rownames %>% column_to_rownames(var="event")
mat[is.na(mat)] <- 0
#mat[mat > 0] <- 1
mat <- as.matrix(mat)
mat <- sqrt(mat)
set.seed(123456)


NMDS1 <-
  metaMDS(mat,
          distance = "bray",
          k = 2,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)

goodness(NMDS1)
stressplot(NMDS1)
plot(NMDS1, type = "t")

#set up grouping variables
data.scores1 = as.data.frame(scores(NMDS1)$sites)

data.scores1$site = cast$site
data.scores1$year = cast$year
data.scores1$event = cast$event

# data.scores1 <- data.scores1 %>% 
#   mutate(survey = fct_relevel(survey,
#                               "January20", "May22", "December22")) 



species.scores1 <- as.data.frame(scores(NMDS1, "species")) 
species.scores1$species <- rownames(species.scores1)  # create a column of species, from the rownames of species.scores

yr1 <- data.scores1[data.scores1$year == "2019", ][chull(data.scores1[data.scores1$year == 
                                                                               "2019", c("NMDS1", "NMDS2")]), ]
yr2 <- data.scores1[data.scores1$year == "2022", ][chull(data.scores1[data.scores1$year == 
                                                                           "2022", c("NMDS1", "NMDS2")]), ]
yr3 <- data.scores1[data.scores1$year == "2024", ][chull(data.scores1[data.scores1$year == 
                                                                                 "2024", c("NMDS1", "NMDS2")]), ]

hull.data1 <- rbind(yr1, yr2, yr3) #%>% 
#   mutate(survey = fct_relevel(survey,
#                               "January20", "May22", "December22")) 


hull.data1


png("output/NMDS_bray.png",width = 7, height = 7, units = "in", res = 400)
ggplot() + 
  geom_polygon(data=hull.data1,aes(x=NMDS1,y=NMDS2,fill=as.factor(year),group=as.factor(year)),alpha=0.30, color = "black") + 
  geom_point(data=data.scores1,aes(x=NMDS1,y=NMDS2,colour = as.factor(year), fill = as.factor(year), shape = site),size=4, alpha = 0.8) +
  geom_point(data=species.scores1,aes(x=NMDS1,y=NMDS2), pch = 20, color = "black", size=2) +
  geom_text_repel(data=species.scores1,aes(x=NMDS1,y=NMDS2,label=species),size = 2, min.segment.length = 0.5, max.overlaps = 18) + 
  #geom_text_repel(data=data.scores1,aes(x=NMDS1,y=NMDS2,label=location_name),size = 4, colour = "black", min.segment.length = 0.5, max.overlaps = 18) + 
  scale_fill_viridis("Year", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  scale_color_viridis("Year", begin = 0.25, end = 0.9, option = "magma", discrete = TRUE) +
  scale_shape_manual("Site",values=c("cbc30c" = 21, "cbc30n" = 22, "cbc30s" = 23, "lagoon" = 24,
                                     "southwater" = 25, "sr_disease" = 3, "sr30c" = 9, "sr30n" = 7, "sr30s" = 10)) +
  #xlim(-1,2) +  
  ggtitle("Bray-Curtis Dissimilarity") +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11), 
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        plot.title = element_text(hjust = 0.5))
dev.off()


cast$year <- as.factor(cast$year)

braydf <- adonis2(mat ~ year, 
                  data = cast,
                  permutation = how(nperm = 9999, blocks = cast$site))
rownames(braydf) <- c("Year", "Residual", "Total")

braytable <- braydf %>%
  kbl(caption = "Table 2: Bray-Curtis dissimilarity, PERMANOVA results") %>%
  kable_styling()

braytable

write.csv(braydf, "output/PERMANOVAresults.csv")

# run simper
simp <- simper(mat, 
               group = cast$year,
               permutations = how(nperm = 9999, blocks = cast$site))

# summary gives you contributing species per pairwise comparison
simp_df <- as.data.frame(summary(simp)[["2019_2022"]])
simp_df1 <- simp_df %>% dplyr::select(p, cumsum, average, sd, ratio, ava, avb) %>%
  #subset(cumsum < 0.71) %>%
  arrange(cumsum)


simptable1 <- simp_df1 %>%
  kbl(caption = "Table 3: SIMPER contributions to Bray-Curtis dissimilarity between 2019 and 2022") %>%
  kable_styling()

simptable1

write.csv(simp_df1, "output/SIMPER_2019to2022.csv")

simp_df2 <- as.data.frame(summary(simp)[["2022_2024"]])
simp_df2 <- simp_df2 %>% dplyr::select(p, cumsum, average, sd, ratio, ava, avb) %>%
  #subset(cumsum < 0.71) %>%
  arrange(cumsum)


simptable2 <- simp_df2 %>%
  kbl(caption = "Table 4: SIMPER contributions to Bray-Curtis dissimilarity between 2022 and 2024
      (up to 70% cumulative sum of contributions)") %>%
  kable_styling()

simptable2

write.csv(simp_df2, "output/SIMPER_2022to2024.csv")

simp_df3 <- as.data.frame(summary(simp)[["2019_2024"]])
simp_df3<- simp_df3 %>% dplyr::select(p, cumsum, average, sd, ratio, ava, avb) %>%
  #subset(cumsum < 0.71) %>%
  arrange(cumsum)


simptable3 <- simp_df3 %>%
  kbl(caption = "Table 5: SIMPER contributions to Bray-Curtis dissimilarity between 2019 and 2024
      (up to 70% cumulative sum of contributions)") %>%
  kable_styling()

simptable3

write.csv(simp_df3, "output/SIMPER_2019to2024.csv")

# Add vertical space between tables using LaTeX \bigskip or a fixed space
spacer <- "<br><br><br>" # use more \\bigskip for even more space

combined_tables <- paste0(braytable, spacer, simptable1, spacer, simptable2, spacer, simptable3)
save_kable(combined_tables, file = "output/combined_tables.pdf")



# # see what comparisons are available
# names(summary(simp))
# 
# # loop through all comparisons into a list of dfs
# simp_results <- lapply(names(summary(simp)), function(comp) {
#   df <- as.data.frame(summary(simp)[[comp]])
#   df$species <- rownames(df)
#   df$comparison <- comp
#   df
# }) %>% bind_rows()
# 
# # filter to species explaining cumulative 70% of dissimilarity
# simp_70 <- simp_results %>% 
#   filter(cumsum <= 0.70)
# 

#remake Bray with only significant simper species, colored by suscep

simp_df4 <- as.data.frame(summary(simp)[["2019_2024"]])
simp_df4<- simp_df4 %>% dplyr::select(p, cumsum, average, sd, ratio, ava, avb) %>%
  subset(p < 0.05) %>%
  rownames_to_column(var = "species")

species.scores1 <- as.data.frame(scores(NMDS1, "species"))
species.scores1$species <- rownames(species.scores1)

species.scores1 <- species.scores1 %>% inner_join(simp_df4, by = "species") 

sus <- demo_newname %>% ungroup() %>% dplyr::select(species, SCTLD_suscep) %>% distinct()

species.scores1 <- species.scores1 %>% inner_join(sus, by = "species")

species.scores1 <- species.scores1 %>%
   mutate(SCTLD_suscep = fct_relevel(SCTLD_suscep,
                              "very high", "high", "med", "low")) 


site <- ggplot() + 
  geom_polygon(data=hull.data1,aes(x=NMDS1,y=NMDS2,fill=as.factor(year),group=as.factor(year)),alpha=0.30, color = "black") + 
  geom_point(data=data.scores1,aes(x=NMDS1,y=NMDS2,colour = as.factor(year), 
                                   fill = as.factor(year), shape = site),size=4, alpha = 0.8) +
  geom_point(data=species.scores1,aes(x=NMDS1,y=NMDS2), pch = 20, color = "black", size=2, alpha = 0) +
  #geom_text_repel(data=species.scores1,aes(x=NMDS1,y=NMDS2,label=species),size = 2, min.segment.length = 0.5, max.overlaps = 18) + 
  #geom_text_repel(data=data.scores1,aes(x=NMDS1,y=NMDS2,label=location_name),size = 4, colour = "black", min.segment.length = 0.5, max.overlaps = 18) + 
  scale_fill_manual(name = "Year",
                      limits = c("2019", "2022", "2024"),
                      values = c("2024" = "#D55E00",
                                 "2022" = "purple3",
                                 "2019" = "yellow4")) +
  scale_colour_manual(name = "Year",
                      limits = c("2019", "2022", "2024"),
                      values = c("2024" = "#D55E00",
                                 "2022" = "purple3",
                                 "2019" = "yellow4"), guide = "none") +
  scale_shape_manual("Site",values=c("cbc30c" = 21, "cbc30n" = 22, "cbc30s" = 23, "lagoon" = 24,
                                     "southwater" = 25, "sr_disease" = 3, "sr30c" = 9, "sr30n" = 7, "sr30s" = 10)) +
  #xlim(-1,2) +  
  #ggtitle("Bray-Curtis Dissimilarity") +
  new_scale_color() +
  coord_equal() +
  theme_bw() +
  theme(  # remove axis ticks
        axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11), 
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.box.just = "left",
        plot.title = element_text(hjust = 0.5))

spec <- ggplot() + 
  
  # Site points by year
  geom_point(data = data.scores1,
             aes(x = NMDS1, y = NMDS2),
             size = 2, alpha = 0) +
  
  # Species points
  geom_point(data = species.scores1,
               aes(x = NMDS1, y = NMDS2, color = SCTLD_suscep), alpha = 0.7) +
  geom_text_repel(data = species.scores1,
                  aes(x = NMDS1, y = NMDS2, label = species, color = SCTLD_suscep),
                  size = 3, min.segment.length = 0.5, max.overlaps = 18, show.legend = FALSE) + 
    scale_color_manual("SCTLD\nSusceptibility",
                       values = c("very high" = "red4",
                                  "high" = "darkorange3",
                                  "med" = "#2CA02C",
                                  "low" = "cadetblue")) +

 # ggtitle("Bray-Curtis Dissimilarity") +
  coord_equal() +
  theme_bw() +
  theme(
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11), 
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))



NMDS_combo <- site + spec #+ 
  # theme(legend.position = "bottom",
  #       legend.box = "vertical",        # stack legend groups vertically
  #       legend.box.just = "left") &     # align them to the left
  # guides(
  #   shape = guide_legend("Site", nrow = 2, order = 3),
  #   fill  = guide_legend("Year", nrow = 1, order = 1),
  #   color = guide_legend("SCTLD\nSusceptibility", nrow = 1, order = 2)
  #)


png("output/NMDS_bray2.png",width = 11, height = 8, units = "in", res = 400)
NMDS_combo
dev.off()


# # Build ellipse polygons using vegan's veganCovEllipse
# library(vegan)
# library(ggrepel)
# library(viridis)
# 
# veganCovEllipse <- function(cov, center = c(0,0), scale = 1, npoints = 100) {
#   theta <- seq(0, 2 * pi, length = npoints)
#   Circle <- cbind(cos(theta), sin(theta))
#   t(center + scale * t(Circle %*% chol(cov)))
# }
# 
# ellipse.data <- data.frame()
# for (g in unique(data.scores1$year)) {
#   group_data <- data.scores1[data.scores1$year == g, c("NMDS1", "NMDS2")]
#   if (nrow(group_data) >= 3) {
#     cov_mat <- cov(group_data)
#     center  <- colMeans(group_data)
#     # scale factor for 95% CI ellipse: sqrt of 95th percentile of chi-squared with 2 df
#     scale   <- sqrt(qchisq(0.95, df = 2))
#     ell     <- veganCovEllipse(cov_mat, center = center, scale = scale)
#     ell_df  <- as.data.frame(ell)
#     colnames(ell_df) <- c("NMDS1", "NMDS2")
#     ell_df$year <- g
#     ellipse.data <- rbind(ellipse.data, ell_df)
#   }
# }
# 
# library(ggnewscale)
# 
# png("Figures/NMDS_bray2.png", width = 7, height = 7, units = "in", res = 400)
# ggplot() + 
#   geom_polygon(data = ellipse.data,
#                aes(x = NMDS1, y = NMDS2, fill = as.factor(year), group = as.factor(year)),
#                alpha = 0.30, color = "black") + 
#   scale_fill_manual("Year", values = c("2019" = "gray10", "2022" = "gray45", "2024" = "gray80")) +
#   
#   # First colour scale: site points by year
#   geom_point(data = data.scores1,
#              aes(x = NMDS1, y = NMDS2, colour = as.factor(year),
#                  fill = as.factor(year), shape = site),
#              size = 2, alpha = 0.5) +
#   scale_color_manual("Year", values = c("2019" = "gray10", "2022" = "gray50", "2024" = "gray80")) +
#   scale_shape_manual("Site", values = c("cbc30c" = 21, "cbc30n" = 22, "cbc30s" = 23,
#                                         "lagoon" = 24, "southwater" = 25, "sr_disease" = 3,
#                                         "sr30c" = 9, "sr30n" = 7, "sr30s" = 10)) +
#   
#   # Reset colour scale for species points
#   new_scale_color() +
#   
#   # Second colour scale: species points by SCTLD_suscep
#   geom_point(data = species.scores1,
#              aes(x = NMDS1, y = NMDS2, color = SCTLD_suscep),
#              pch = 20, size = 4) +
#   scale_color_manual("SCTLD\nSusceptibility",
#                      values = c("very high" = "#D62728",
#                                 "high" = "#FF7F0E", 
#                                 "med" = "yellow",
#                                 "low" = "#2CA02C")) +
#   
#   geom_text_repel(data = species.scores1,
#                   aes(x = NMDS1, y = NMDS2, label = species),
#                   size = 2, min.segment.length = 0.5, max.overlaps = 18) + 
#   ggtitle("Bray-Curtis Dissimilarity") +
#   coord_equal() +
#   theme_bw() +
#   theme(axis.text.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title.x = element_text(size = 11),
#         axis.title.y = element_text(size = 11), 
#         legend.text = element_text(size = 11),
#         legend.title = element_text(size = 12),
#         panel.background = element_blank(), 
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         plot.background = element_blank(),
#         plot.title = element_text(hjust = 0.5))
# dev.off()

