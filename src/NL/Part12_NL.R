
## Part12. Plotting step bar plot for unsmoothed SIRs at the age 15+ in 2001 by six population percentiles

install.packages("dplyr")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("sf")

library(dplyr)
library(tidyverse)
library(ggplot2)
library(sf)


setwd("U:/Documents/Chornobyl/Codes/2023/May_2023/Short/Part12_NL")
load("incidence_final.RData")


# Dataframe, works
SIR15 <- data.frame(
  sex = rep(c("female", "male"), each = 6),
  quan = rep(c(10, 25, 50, 75, 90, 100), times = 2),
  sir = c(0, 0.21011, 0.63866, 0.98992, 1.68309, 3.011246, 0, 0, 0.04192, 1.07594, 2.159034, 3.4407),
  perc = rep(c(0, 10, 25, 50, 75, 90), times = 2)) 

# Add dummy variable and extra row for geom_step

SIR15 <- SIR15 %>%
  bind_rows(SIR15%>% filter(quan == 100) %>% mutate(perc = perc + 10))


ggplot() +
  scale_x_continuous(breaks = c(0, 10, 25, 50, 75, 90, 100)) + 
  geom_step(data = SIR15 %>% filter(sex == "female"), 
            aes(x = perc, y = sir, color="females"), linewidth = 1.5) +
  geom_step(data = SIR15%>% filter(sex == "male"),
            aes(x = perc, y = sir, color="males"), linewidth = 2, alpha = .5) +
  scale_color_manual(values=c("females"="#f50057", "males"="#304ffe"))+
  labs(x = "Population proportion", y = "SIR") +
  theme_light() +
  theme(legend.title=element_blank(),
        legend.position="right", 
        axis.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())
