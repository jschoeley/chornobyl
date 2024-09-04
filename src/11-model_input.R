# Prepare Ukrainian data on regional thyroid cancer incidence
# and radiation exposure for modeling

# Init ------------------------------------------------------------

library(yaml)
library(tidyverse)
library(sf)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  glob = 'src/00-global_functions.R',
  config = 'cfg/config.yaml',
  maptemplates = 'out/10-maptemplates.rds',
  # thyroid cancer incidence 2001 and radiation exposure 1986
  thyroid = 'dat/thyroid.rds'
)
paths$output <- list(
  modelinput = 'out/11-modelinput.rds'
)

# global configuration
config <- read_yaml(paths$input$config)

# global objects
source(paths$input$glob)

# list containers for analysis artifacts
dat <- list()
modelinput <- list()

# Input data ------------------------------------------------------

dat$maptemplates <- readRDS(paths$input$maptemplates)

dat$thyroid <-
  readRDS(paths$input$thyroid) |>
  # we exclude the regions where incidence is NA
  filter(!id_inc %in% c(14230, 44242)) |>
  select(region_id = id_inc, age = age5, sex,
         incidence_2001, population_2001, dose)

# Add total sex category ------------------------------------------

dat$thyroid_total_sex <-
  dat$thyroid |>
  group_by(region_id, age) |>
  summarise(
    sex = 'total',
    dose = sum(dose*population_2001)/sum(population_2001),
    incidence_2001 = sum(incidence_2001),
    population_2001 = sum(population_2001)
  ) |>
  ungroup() |>
  select(region_id, age, sex, incidence_2001, population_2001, dose)

dat$thyroid <-
  bind_rows(dat$thyroid, dat$thyroid_total_sex) |>
  arrange(region_id, age, sex)

# Expand data for modeling purposes -------------------------------

# add expected and observed incidences as basis for crude and
# smoothed SIR calculations

model_input <- list()

# calculate district centroids
model_input$centroids <-
  bind_cols(region_id = dat$maptemplates$ukrgeo$region_id,
            st_coordinates(st_centroid(dat$maptemplates$ukrgeo)))

# calculate ukrainian average
# age-sex-specific thyroid cancer incidence
model_input$national_incidence <-
  dat$thyroid |> 
  group_by(age, sex) |>
  summarise(
    national_incidence_2001 = sum(incidence_2001),
    national_population_2001 = sum(population_2001),
    national_rate = national_incidence_2001/national_population_2001
  ) |>
  ungroup()

# calculate total 2001 population by sex and region
model_input$population_2001_by_region_sex <-
  dat$thyroid |>
  group_by(region_id, sex) |>
  summarise(population_2001 = sum(population_2001)) |>
  ungroup()

# calculate average dosage received by region and sex
# weighted by age-distribution
model_input$avg_dosage_by_region_sex <-
  dat$thyroid |>
  group_by(region_id, sex) |>
  summarise(
    # age-structure weighted average dosage
    average_dose = sum(dose*population_2001)/sum(population_2001)
  ) |>
  ungroup()

model_input$national_avg_dosage <-
  dat$thyroid |>
  summarise(
    average_dose = sum(dose*population_2001)/sum(population_2001)
  )

# expected incidences and crude SIRs by region, sex, and age
model_input$expected_observed_by_region_age_sex <-
  dat$thyroid |>
  left_join(model_input$national_incidence, by = c('age', 'sex')) |>
  mutate(
    # calculate expected incidence by region, sex and age,
    # given that ukrainian avg. thyroid cancer rates apply
    incidence_expected = population_2001*national_rate,
    incidence_observed = incidence_2001
  )

# expected incidences and raw SIRs by region, and sex
model_input$expected_observed_by_region_sex <-
  model_input$expected_observed_by_region_age_sex |>
  group_by(region_id, sex) |>
  # sum over age
  summarise(
    incidence_expected = sum(incidence_expected),
    incidence_observed = sum(incidence_observed)
  ) |>
  ungroup()

# Merge model input data ------------------------------------------

model_input$region_sex <-
  model_input$expected_observed_by_region_sex |>
  # add average dosage
  left_join(model_input$avg_dosage_by_region_sex) |>
  # add population
  left_join(model_input$population_2001_by_region_sex) |>
  # add region centroids
  left_join(model_input$centroids)

model_input$region_sex_age <-
  dat$thyroid |>
  # add region centroids
  left_join(model_input$centroids)

# Export ----------------------------------------------------------

saveRDS(model_input, paths$output$modelinput)
