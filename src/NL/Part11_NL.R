
## Part11. Calculate unsmoothed district SIRs of thyroid cancer 15+ in 2001 by six population percentiles

install.packages("dplyr")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("sf")

library(dplyr)
library(tidyverse)
library(ggplot2)
library(sf)

setwd ("U:/Documents/Chornobyl/Codes/2023/May_2023/Short/Part11_NL")
load("thyroid_15.RData") ## dose_incidence

## There are two districts with no data on incidence 
## Lymanskyi raion (id_inc=14230), Dovzhanskyi raion (id_inc=44242), 

## standardized thyroid incidence ratio, 15+, by sex (per 100,000)

dose_incidence1<- dose_incidence%>%
  group_by(id_inc, sex) %>% 
  mutate(si_sum15=sum(si_product)/100000,
         thy_sum15=sum(incidence_2001),
         thy_sir15= thy_sum15/si_sum15,
         pop15=sum(population_2001)) %>% 
  ungroup()


##remove some columns
dose_incidence2 <- subset(dose_incidence1, select=-c(type_c, si_product, si_sum, thy_sum, si_product1, si_sum15_29, thy15_29, thy_sum15_29, thy_sir15_29, thy_sr, merge_age, si_sum15, thy_sum15))

## sort/arrange thy_sir15 in ascending order along with other variables
dose_incidence3 <- dose_incidence2 %>% 
  group_by(id_inc, sex) %>% 
  arrange(thy_sir15) %>% 
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
  #this next mutate command creates the quantiles for each group of population
  mutate(quant=cut(cum_pop15_pr, breaks=c(0, .10, .25, .50, .75, .90, 1), right = FALSE, include.lowest = TRUE))%>% 
  ungroup()

#create a lookup table of the quants to join back the quantile to all ages of dose_incidence file for the district 
quant_lookup = dose_incidence8 %>% 
  filter(age_15_or_not==1) %>% 
  dplyr:: select(id_inc,quant,sex) %>% 
   distinct()

dose_incidence9 <- dose_incidence8 %>% dplyr::select(-quant) %>% left_join(quant_lookup)

summary(dose_incidence1$thy_sir15[dose_incidence1$sex=="male"], basic = T)
summary(dose_incidence1$thy_sir15[dose_incidence1$sex=="female"], basic = T)

summary(dose_incidence8$quant[dose_incidence8$sex=="male"], basic = T)
summary(dose_incidence8$quant[dose_incidence8$sex=="female"], basic = T)


##calculating SIRs by quant 

## standard thyroid incidence rate, Ukraine, 2001, by sex (per 100,000)

Ukr_rate_m15<-c(1.19906, 0.90070, 1.11316, 1.25663, 1.73949, 2.22606, 2.20415, 3.88847, 3.22142, 4.04786, 3.08679, 2.83050, 3.80332, 2.00903, 0.0)
Ukr_rate_f15<-c(2.30235, 4.27059, 6.20360, 6.16168, 8.69239, 9.85421, 12.94928, 15.55536, 10.98955, 12.57365, 7.95880, 9.04459, 5.88348, 1.90925, 1.40940)

#the data.frame to join the sir of Ukraine by sex and age to the dose incidence file.
std_inc_rate = data.frame(sex=rep(c("male","female"),each=15),age5=rep(unique(dose_incidence$age5),times=2),std_rate = c(Ukr_rate_m15,Ukr_rate_f15))

#join standardized incidence rate and calculated expected cases
dose_incidence11<- dose_incidence9 %>% left_join(std_inc_rate) %>% #here is where i join the sir to the file
  group_by(id_inc, sex)  %>%
  mutate(si_product_q=ifelse(sex=="male", population_2001*std_rate, population_2001*std_rate)) %>% ungroup()

#calculate sir by quantile and sex
dose_incidence12<- dose_incidence11 %>%
  group_by(sex, quant) %>% 
  mutate(si_sum_q=sum(si_product_q)/100000,
         thy_sum_q=sum(incidence_2001, na.rm=T),
         thy_sir_q=(thy_sum_q/si_sum_q))%>% 
  ungroup()

summary(dose_incidence12$thy_sir15[dose_incidence12$sex=="male"], basic = T)
summary(dose_incidence12$thy_sir15[dose_incidence12$sex=="female"], basic = T)


# this summarizes the dataframe by group so that it just produces the information we need for the plot

dose_incidence_sum<- dose_incidence11 %>%
  group_by(sex, quant) %>% 
  summarise(si_sum_q=sum(si_product_q)/100000,
         thy_sum_q=sum(incidence_2001, na.rm=T),
         thy_sir_q=(thy_sum_q/si_sum_q))%>% 
  ungroup()

save(dose_incidence_sum, file="incidence_final.Rdata")
