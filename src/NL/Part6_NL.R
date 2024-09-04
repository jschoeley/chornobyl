
##Part6. Spatial interpolation of settlement dosage data to get district mean doses and
## interpolation of population weighted doses at the age 15-19


# Load Libraries

packages= c("dplyr","tidyr","ggplot2","viridis","spatialreg","stars", "scales",
            "terra","spdep","rgdal","patchwork","RColorBrewer","sf","maptools",
            "raster","fields","gstat","readr","stringr","readxl","data.table")
new.packages = packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
unlist(lapply(packages, require, character.only=TRUE))

setwd("U:/Documents/Chornobyl/Codes/2023/May_2023/Short/Part6_NL")


# Read dosage data
load("valid_st_data_for_interp.RData") ## "doses" and "new_dist" 
##has two files: "doses" and "new dist", both already have geometry included.

##"doses" has doses by original Age, and Ua_Old_id, Area from the original shapefile
## "new_dist" has population and incidence by original Age5, and merge_id,  type_c


# Read new disrict data file
load("new_dist.RData") ## "new_dist_shape"
##has all variables and thyroid SIR for all ages combined, id_inc (this is our final ID variable) and geometry

# subset 1 instance of each id

shapefile_dat = data.frame(id_inc = unique(new_dist_shape$id_inc)) %>%
  left_join(new_dist_shape[,1],by=c("id_inc"="id_inc"))

# Check NA values for dose
table(is.na(doses$dose)) # 12 NA values ##refer to the city of Slavutysh built after Chornobyl accident

# Remove these NA values so that interpolation regions are all filled.
doses_noNA = doses %>% filter(!is.na(dose))

# REAL INTERPOLATION

interp_age_female=list()
for(k in unique(doses$age)){
  dose_subset = doses_noNA %>% filter(sex=="f"&age==k)
  dist_subset = new_dist_shape %>% filter(age5==k & sex == 2)
  interp_age_female[[paste0("female_",k)]]=st_interpolate_aw(dose_subset["dose"],dist_subset,extensive = F,keep_NA=T)
  interp_age_female[[paste0("female_",k)]]$sex="female"
  interp_age_female[[paste0("female_",k)]]$age=as.numeric(k)
  interp_age_female[[paste0("female_",k)]]$age2001=as.numeric(k)+15
  interp_age_female[[paste0("female_",k)]]$merge_id=dist_subset$id_inc
  print(paste("interp female",k, "completed"))
  
  #save here so that any results pre-r crashing are ok  
  save(interp_age_female, file="dose_interp_females.RData")
}

interp_age_female = rbindlist(interp_age_female)
save(interp_age_female, file="dose_interp_females.RData")


interp_age_male=list()
for(k in unique(doses$age)){
  dose_subset = doses_noNA %>% filter(sex=="m"&age==k)
  dist_subset = new_dist_shape %>% filter(age5==k & sex == 1)
  interp_age_male[[paste0("male_",k)]]=st_interpolate_aw(dose_subset["dose"],dist_subset,extensive = F,keep_NA = T)
  
  interp_age_male[[paste0("male_",k)]]$sex="male"
  interp_age_male[[paste0("male_",k)]]$age=as.numeric(k)
  interp_age_male[[paste0("male_",k)]]$age2001=as.numeric(k)+15
  interp_age_male[[paste0("male_",k)]]$merge_id=dist_subset$id_inc
  print(paste("interp male",k, "completed"))
  
  #save here so that any results pre-r crashing are ok  
  save(interp_age_male, file="dose_interp_males.RData")
}
interp_age_male = rbindlist(interp_age_male)
save(interp_age_male, file="dose_interp_males.RData")


load("dose_interp_females.RData")

#combine dose data into one
dose_data = rbind(interp_age_female,interp_age_male)


###########################################################################################################
################interpolation of population weighted doses at the age 15-19

#load the weightings for ages 15-19
load("weights_15_19.RData")

#separate 15-19 from dosage
dose_15_19 =  dose_data %>% 
  filter(age2001%in%15:16)

dose_data_no15_19 = dose_data %>% 
  filter(!age2001%in%15:16) %>% 
  dplyr::select(-age)

#subset unique geometry to merge back
unique_geom = dose_15_19 %>% 
  dplyr::select(merge_id,geometry) %>% 
  distinct()

weights_15_19 = weights_15_19 %>% 
  mutate(sex=ifelse(sex==1,"male","female"))

#combine dose data with weights, calculated weighted dosage for 15-9 age group for sex and merge/inc_id
dose_15_19 = left_join(dose_15_19,weights_15_19,by=c("merge_id"="id_inc","sex"="sex","age2001"="age_group")) %>% 
  mutate(dose_pct=dose*pct_15_19) %>% group_by(sex,merge_id) %>%
  summarise(dose = sum(dose_pct)) %>% ungroup() %>% mutate(age2001=15)
dose_15_19= left_join(dose_15_19,unique_geom)

#merge the dosage data to get new overall dosage data

dose_dat_combined = rbind(dose_15_19,dose_data_no15_19)


new_dist_shape = new_dist_shape %>% mutate(sex = ifelse(sex==1,"male","female"))

#do some checks: correct unique ids: 627
length(unique(dose_dat_combined$merge_id))

#627 obs (one per dist) for each sex/age
table(dose_dat_combined$age2001,dose_dat_combined$sex)

#combine dosage data with the incidence data
names(new_dist_shape)
names(dose_dat_combined)
class(new_dist_shape$age5)
class(dose_dat_combined$age2001)

#create variable by which to merge age if the population 2001 age is higher than 35 (so they should get dosage of age 20)
new_dist_shape = new_dist_shape %>% 
  mutate(merge_age=ifelse(age5>35,35,age5))

#merge shapefile/incidence with areal weighted dosage be age5 (incidence in 2001 by population) and age 2001 (age in dose + 15 years to get age in 2001)
dose_incidence_shapefile = left_join(new_dist_shape,dose_dat_combined %>% 
        dplyr::select(-geometry),by=c("merge_age"="age2001","sex"="sex","id_inc"="merge_id"))


#save new data
save(dose_incidence_shapefile, file="dose_incidence_shapefile.RData")

#check to make sure the NAs are there/correct
dose_incidence_shapefile %>% filter(id_inc=="44242")

#check to make sure the correct ages have a corresponding dosage amount
table(is.na(dose_incidence_shapefile$dose),dose_incidence_shapefile$age5)

#
dose_incidence_shapefile %>% filter(id_inc=="44242")
