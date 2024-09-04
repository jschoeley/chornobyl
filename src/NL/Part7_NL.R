
## Part7. Measuring spatial autocorrelation, Moran I, and Lisa plotting, unsmoothed SIRs

install.packages("sf")
install.packages("spdep")
install.packages("spatialreg")
install.packages("rgeoda")
install.packages("tidyverse")
install.packages("ggpubr")
install.packages("ggrepel")

library(sf)
library(spdep)
library(spatialreg)
library(tidyverse)
library(rgeoda)
library(ggpubr)
library(ggrepel)

setwd("U:/Documents/Chornobyl/Codes/2023/May_2023/Short/Part7_NL")

## Measuring Spatial Autocorrelation
load("thy_sir15_m.RData") ## has only id_inc, SIR15+ and geometry, males 15+
load("thy_sir15_f.RData") ## has only id_inc, SIR15+ and geometry, females 15+
load("lisa.colors")

##males
data_M1 <- data_M %>%
  filter(!is.na(thy_sir15)) 
queen.nb <- poly2nb(data_M1, queen = TRUE)
listw.counties <- nb2listw(queen.nb, zero.policy = TRUE)
str(listw.counties, list.len = 5)

# Compute the global Moran I test statistic and plot a Moran I scatter plot
moran.test(x = data_M1$thy_sir15, listw = listw.counties, zero.policy = TRUE, alternative = "two.sided")
moran.plot(x = data_M1$thy_sir15, listw = listw.counties, zero.policy = TRUE, 
           plot = TRUE, xlab = "SIR15+", ylab = "Spatially lagged SIR15+", 
           labels = FALSE, col = "orangered", pch = 20, cex = 0.1)


##females
data_F1 <- data_F %>%
  filter(!is.na(thy_sir15)) 

queen.nb <- poly2nb(data_F1, queen = TRUE)
listw.counties <- nb2listw(queen.nb, zero.policy = TRUE)
str(listw.counties, list.len = 5)

# Compute the global Moran I test statistic and plot a Moran I scatterplot
moran.test(x = data_F1$thy_sir15, listw = listw.counties, zero.policy = TRUE, alternative = "two.sided")
moran.plot(x = data_F1$thy_sir15, listw = listw.counties, zero.policy = TRUE, 
           plot = TRUE, xlab = "SIR15+", ylab = "Spatially lagged SIR15+", 
           labels = FALSE, col = "orangered", pch = 20, cex = 0.1)


length(unique(dose_incidence9$thy_sir15[dose_incidence9$sex=="male"], basic = T))
length(unique(data_M$thy_sir15))
summary(dose_incidence9$thy_sir15[dose_incidence9$sex=="male"])
summary(dose_incidence9$thy_sir[dose_incidence9$sex=="male"])
summary(data_M$thy_sir15)

summary(data_F$thy_sir15)
summary(dose_incidence9$thy_sir15[dose_incidence9$sex=="female"])
summary(dose_incidence9$thy_sir[dose_incidence9$sex=="female"])


# Compute LISA statistics (aka Local Moran) and determine the clustering
##males

counties.w <- queen_weights(data_M1, order = 1)  # Queen weights (different structure required by rgeoda)

# Computing actual Moran scores and related statistics
lisa <- local_moran(df = data_M1["thy_sir15"], w = counties.w)

# Extracting labels and assigning the colors (all contained in the object calculated above, just need extraction)
lisa.labels <- lisa_labels(lisa)
lisa.colors <- setNames(lisa_colors(lisa), lisa.labels)

# Generate cluster type designation for each county in the data
lisa.clusters.M <- data_M1 %>% 
  mutate(cluster.num = lisa_clusters(lisa) + 1,  ## LISA marks clusters starting at 0, so adding 1 here for clarity
         cluster = factor(lisa.labels[cluster.num], levels = lisa.labels))

# Map the clusters on the map
ggplot(lisa.clusters.M, aes(fill = cluster)) +
  geom_sf(color = "white", size = 0) +
  scale_fill_manual(values = lisa.colors, na.value = "black") +
  guides(fill = guide_legend(title = "SIRs 15+, males")) +
  theme_void()

save(lisa.clusters.M, file="lisa.clusters_M.RData")

load("chornobyl.RData")
##this plot does not work
ggplot(lisa.clusters.M, aes(fill = cluster)) +
  geom_sf(color = "white", size = 0) +
  scale_fill_manual(values = lisa.colors, na.value = "black") +
  guides(fill = guide_legend(title = "SIRs 15+, males")) +
  geom_point(data = Chornobyl, mapping = aes(x = X, y = Y),
             shape = 15,
             color = "red",
             fill = "red",
             size = 3,
             stroke = 1, # width of the points border
  ) +
  geom_text_repel(data = Chornobyl, aes(x = X, y = Y, label = UnitName),
                  color = "red", size = 7, fontface = "bold",
                  bg.color = "yellow", bg.r = 0.15) +
  theme_void()


##females

counties.w <- queen_weights(data_F1, order = 1)  # Queen weights (different structure required by rgeoda)

# Computing actual Moran scores and related statistics
lisa <- local_moran(df = data_F1["thy_sir15"], w = counties.w)

# Extracting labels and assigning the colors (all contained in the object calculated above, just need extraction)
lisa.labels <- lisa_labels(lisa)
lisa.colors <- setNames(lisa_colors(lisa), lisa.labels)

# Generate cluster type designation for each county in the data
lisa.clusters.F <- data_F1 %>% 
  mutate(cluster.num = lisa_clusters(lisa) + 1,  ## LISA marks clusters starting at 0, so adding 1 here for clarity
         cluster = factor(lisa.labels[cluster.num], levels = lisa.labels))


save(lisa.clusters.F, file="lisa.clusters_F.RData")

# Map the clusters on the map
ggplot(lisa.clusters.F, aes(fill = cluster)) +
  geom_sf(color = "white", size = 0) +
  scale_fill_manual(values = lisa.colors, na.value = "black") +
  guides(fill = guide_legend(title = "SIRs 15+, females")) +
  theme_void()

