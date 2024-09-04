
## Part10. Plotting step bar plot for dosage data (1986 dosage15+ by six population percentiles)

install.packages("dplyr")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("sf")

library(dplyr)
library(tidyverse)
library(ggplot2)
library(sf)


setwd("U:/Documents/Chornobyl/Codes/2023/May_2023/Short/Part10_NL")
load("dose_final.RData")


# Create dataframe based on dose_final.RData

Dose15 <- data.frame(
  sex = rep(c("female", "male"), each = 6),
  quan = rep(c(10, 25, 50, 75, 90, 100), times = 2),
  dose = c(16.40936, 20.72791, 28.43919, 41.17954, 57.38597, 110.67292, 17.38992, 21.99927, 30.18327, 44.34208, 62.50075, 117.45269),
  perc = rep(c(0, 10, 25, 50, 75, 90), times = 2)) 

# Add dummy variable and extra row for geom_step



ggplot() +
  scale_x_continuous(breaks = c(0, 10, 25, 50, 75, 90, 100)) + 
  geom_step(data = Dose15 %>% filter(sex == "female"), 
            aes(x = perc, y = dose, color="females"), linewidth = 1.5) +
  geom_step(data = Dose15%>% filter(sex == "male"),
            aes(x = perc, y = dose, color="males"), linewidth = 2, alpha = .5) +
  scale_color_manual(values=c("females"="#ff5722", "males"="blue"))+
  labs(x = "Population percentage", y = "Dosage") +
  theme_light() +
  theme(legend.title=element_blank(),
        legend.position="right", 
        axis.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

