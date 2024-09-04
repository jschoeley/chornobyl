

### Part9. Calculating district average thyroid doses in 1986 by six population percentiles 
## (0, .10, .25, .50, .75, .90, 1) 


install.packages("dplyr")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("sf")

library(dplyr)
library(tidyverse)
library(ggplot2)
library(sf)


setwd("U:/Documents/Chornobyl/Codes/2023/May_2023/Short/Part9_NL")

load("thyroid_15.RData") ## this is a file with final interpolated doses by age 


## calculate population and age standardized dosage at age 15+ combined, by sex 

dose_incidence1<-dose_incidence%>% 
  group_by(id_inc, sex)%>% 
  mutate(
    pop15=sum(population_2001), 
    dose_mean15 = sum(dose*population_2001)/pop15) %>%
  ungroup()

##remove some columns
dose_incidence2 <- subset(dose_incidence1, select=-c(type_c, si_product, si_sum, thy_sum, si_product1, si_sum15_29, thy15_29, thy_sum15_29, thy_sir15_29, thy_sr, merge_age))

## sort/arrange dose_mean15 in ascending order along with other variables
dose_incidence3 <- dose_incidence2 %>% 
  group_by(id_inc, sex) %>% 
  arrange(dose_mean15) %>% 
  ungroup()

##calculate aggregated sums of pop15, by sex

dose_incidence4<-dose_incidence3%>% 
  mutate(pop15_0=ifelse(age5==15, pop15, 0))

dose_incidence5<- dose_incidence4 %>%
  group_by(sex) %>% 
  mutate(cum_pop15=cumsum(pop15_0)) %>% 
  ungroup()

rm(dose_incidence2, dose_incidence3, dose_incidence4)

##calculating cumulative proportions of cum_pop15: 

dose_incidence6<-dose_incidence5%>% 
  mutate(cum_pop15_0=ifelse(age5==15, cum_pop15, 0))

dose_incidence7<- dose_incidence6 %>%
  group_by(sex) %>% 
  mutate(tot_pop = sum(pop15_0),
         cum_pop15_pr= (cum_pop15_0)/tot_pop) %>% 
  ungroup()

## cut with breaks: 10, 25, 50, 75, 90%
#here, we first create a variable that says if the age5==15 or not, then we create quantiles for each of these (age5==15) or not groups
#(we are only interested in the age5==15 group, so the other one really doesn't matter)
dose_incidence8<- dose_incidence7 %>%
  mutate(age_15_or_not = ifelse(age5==15,1,0)) %>% #create the binary variable
  group_by(sex, age_15_or_not) %>%  #group by if the age==15 (our variable of interest) or not and sex 
  #this mutate command creates the quantiles for each group of population
  mutate(quant=cut(cum_pop15_pr, breaks=c(0, .10, .25, .50, .75, .90, 1), right = FALSE, include.lowest = TRUE))%>% 
  ungroup()

#create a lookup table of the quants to join back the quantile to all ages of dose_incidence file for the district 
quant_lookup = dose_incidence8 %>% 
  filter(age_15_or_not==1) %>% 
  dplyr::select(id_inc, quant, sex) %>% 
  distinct()

dose_incidence9 <- dose_incidence8 %>% dplyr::select(-quant) %>% left_join(quant_lookup)

summary(dose_incidence8$dose_mean15[dose_incidence8$sex=="male"], basic = T)
summary(dose_incidence8$dose_mean15[dose_incidence8$sex=="female"], basic = T)

##calculate dose15 by quant and sex

dose_incidence10<- dose_incidence9 %>%
  group_by(id_inc, sex)  %>%
  mutate(dose_mean15_age = dose*population_2001) %>% 
  ungroup()


dose_incidence11<- dose_incidence10 %>%
  group_by(sex, quant) %>% 
  mutate(
    pop_total = sum(population_2001),
    dose_mean15_cum=sum(dose_mean15_age),
    dose_mean15_q = dose_mean15_cum/pop_total) %>%
  ungroup()


summary(dose_incidence11$dose_mean15_q[dose_incidence11$sex=="male"], basic = T)
summary(dose_incidence11$dose_mean15_q[dose_incidence11$sex=="female"], basic = T)

##Summarize doses I use for plotting

dose_sum<- dose_incidence10 %>%
  group_by(sex, quant) %>% 
  summarise(pop_total = sum(population_2001),
            dose_mean15_cum=sum(dose_mean15_age),
            dose_mean15_q = dose_mean15_cum/pop_total) %>%
  ungroup()

save(dose_sum, file="dose_final.RData")
