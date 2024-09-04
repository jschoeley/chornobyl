
##Part14. Summarizing Ukraine's population by sex and 5-year age groups
## And calculating standardized incidence rates at the 15+ by sex using Ukraine's population weights, by two sexes combined

packages= c("dplyr","tidyr","ggplot2","viridis","spatialreg","stars", "scales",
            "terra","spdep","rgdal","patchwork","RColorBrewer","sf","maptools",
            "raster","fields","gstat","readr","stringr","readxl","data.table")
new.packages = packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
unlist(lapply(packages, require, character.only=TRUE))


setwd("D:/Institute/Zavdannia_2022/Redim_Chornobyl/Data/Codes/2023/May_2023/Short/Part14_NL")

load("dose_incidence_shapefile.RData")


    dose_incidence<-dose_incidence_shapefile  %>% 
    st_drop_geometry()
 
    ##calculating weights for both sexes combined   

population_Ukr1<-dose_incidence%>% 
  group_by(age5) %>% 
  summarize(ukr_pop=sum(population_2001))%>% 
  mutate(sumpop=sum(ukr_pop),  
         weight=ukr_pop/sumpop,
         weight_r=round(weight, digits=3),
         total_w=sum(weight_r))%>% 
  ungroup() %>% 
  dplyr::select(age5, weight_r) %>% 
  distinct()


save(population_Ukr1, file="Pop2001_weights.RData")

##sum(weight at the age15+)=0.835



#join population_Ukr and dose_incidence_shapefile
dose_incidence_w<- left_join(dose_incidence_shapefile, population_Ukr1, by=c("age5"="age5")) %>% 
  st_as_sf()
                             
class(dose_incidence_w)


##calculating SIrates at the age 15+ by sex
## using 2001 Ukraine's population as a standard, by two sexes combined


dose_incidence_sdr<- dose_incidence_w %>%
  filter(age5>=15) %>% 
  group_by(id_inc, sex) %>% 
  mutate(kx=(incidence_2001/population_2001)*100000) %>% 
  ungroup() %>% 
  group_by(id_inc, sex, age5) %>% 
  mutate(kx_product=kx*weight_r) %>%
  ungroup() %>% 
  group_by(id_inc, sex) %>% 
  mutate(kx_sum = sum(kx_product),
  sdr15 = kx_sum/0.835)%>%
  ungroup()

save(dose_incidence_sdr, file="incidence_rates15.RData")


