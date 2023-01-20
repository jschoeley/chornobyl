# Estimate SIR of Thyroid cancer incidence in Ukraine 2001
# and correlate with radiation exposure received in 1986

# Init ------------------------------------------------------------

library(yaml)
library(tidyverse)
library(sf)
library(mgcv)
# detect spatial autocorrelation via residual dignostics
library(DHARMa)
library(mgcViz)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  tmpdir = 'tmp',
  glob = 'src/00-global_functions.R',
  config = 'cfg/config.yaml',
  # district level geodata (sf) of Ukraine
  ukrgeo = 'dat/ukrgeo.rds',
  # thyroid cancer incidence 2001 and radiation exposure 1986
  thyroid = 'dat/thyroid.rds'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  out = 'out/'
)

# global configuration
config <- read_yaml(paths$input$config)

# global objects
source(paths$input$glob)

# constants specific to this analysis
cnst <- within(list(), {
  sir_breaks = c(0, 1/4, 1/3, 1/2, 1, 2, 3, 4, Inf)
  sir_colors = c(`[-Inf,0.25)`='#2166ac',
                 `[0.25,0.333)`='#4393c3',
                 `[0.333,0.5)`='#92c5de',
                 `[0.5,1)`='#d1e5f0',
                 `[1,2)`='#fddbc7',
                 `[2,3)`='#f4a582',
                 `[3,4)`='#d6604d',
                 `[4,Inf)`='#b2182b')
  sir_labels = c('<1/4', '1/4 to 1/3', '1/3 to 1/2', '1/2 to 1', '1 to 2', '2 to 3', '3 to 4', '>4')
})

# list containers for analysis artifacts
dat <- list()
fig <- list()

# Import data -----------------------------------------------------

dat$ukrgeo <- readRDS(paths$input$ukrgeo)
dat$thyroid <- readRDS(paths$input$thyroid)

# Create spatial outline of Ukraine -------------------------------

dat$ukr_outline_geo <-
  dat$ukrgeo %>%
  # avoid rendering artifacts
  st_buffer(0.0001) %>%
  st_make_valid() %>%
  summarise(id = 'ukr') %>%
  st_union()

# District centroids ----------------------------------------------

dat$centroids <-
  bind_cols(id_inc = dat$ukrgeo$id_inc, st_coordinates(st_centroid(dat$ukrgeo)))

# Prepare thyroid data --------------------------------------------

dat$thyroid <-
  dat$thyroid %>%
  # we exclude the regions where incidence is NA
  filter(!id_inc %in% c(14230, 44242)) %>%
  # add lon and lat coordinates
  left_join(dat$ukrgeo, by = 'id_inc') %>%
  select(-geometry)

# Estimate SIR ----------------------------------------------------

sir <- list()

# calculate ukrainian average
# age-sex-specific thyroid cancer incidence
sir$ukr_avg_inc <-
  dat$thyroid %>%
  group_by(age5, sex) %>%
  summarise(
    ukr_inc = sum(incidence_2001),
    ukr_pop = sum(population_2001),
    ukr_incrate = ukr_inc/ukr_pop
  ) %>%
  ungroup()

# calculate average dosage received by region and sex
# weighted by age-distribution
sir$avg_dosage_by_region <-
  dat$thyroid %>%
  group_by(id_inc, sex) %>%
  summarise(
    # age-structure weighted average dosage
    dos = sum(dose*population_2001)/sum(population_2001)
  ) %>%
  ungroup()

# expected incidences and crude SIRs by region, sex, and age
sir$sir_by_region_age_sex <-
  dat$thyroid %>%
  left_join(sir$ukr_avg_inc, by = c('age5', 'sex')) %>%
  mutate(
    # calculate expected incidence by region, sex and age,
    # given that ukrainian avg. thyroid cancer rates apply
    expected = population_2001*ukr_incrate,
    sir_est = incidence_2001/expected,
    sir_q025 = SIRCI(incidence_2001, expected, 0.025),
    sir_q975 = SIRCI(incidence_2001, expected, 0.975)
  )

# expected incidences and raw SIRs by region, and sex
sir$sir_by_region_sex <-
  sir$sir_by_region_age_sex %>%
  group_by(id_inc, sex) %>%
  # sum over age
  summarise(
    expected = sum(expected),
    observed = sum(incidence_2001),
    sir_est = observed/expected,
    sir_q025 = SIRCI(observed, expected, 0.025),
    sir_q975 = SIRCI(observed, expected, 0.975)
  ) %>%
  ungroup()

sir$ready_for_model <-
  sir$sir_by_region_sex %>%
  # add average dosage
  left_join(sir$avg_dosage_by_region) %>%
  # add region centroids
  left_join(dat$centroids)

# Model SIR via spatial smoothing NB regression -------------------

sir$predicted_sir <-
  sir$ready_for_model %>%
  group_by(sex) %>%
  group_modify(~{
    
    model_dat <- .x
    logexpected <- log(model_dat$expected)
    
    fit <-
      gam(
        formula =
          observed ~
          s(Y, X, bs = 'tp') + offset(log(expected)),
        family = nb(link = 'log'),
        data = model_dat
      )
    
    autocor <- testSpatialAutocorrelation(
      simulationOutput = simulateResiduals(fittedModel = fit),
      x = fit$model$X, y = fit$model$Y
    )[['p.value']]
    
    eta <- predict(
      fit, type = 'link', se.fit = TRUE, unconditional = TRUE
    )
    
    # eta, the linear predictor, is on the log incidence scale,
    # subtract log(expected) to get to the log SIR scale
    logsir <- eta$fit - logexpected
    # SD(X-c) = SD(X), thus the standard error also applies
    # when we put the linear predictor on the log SIR scale by
    # subtracting log(expected)
    logsirQ025 <- logsir - 1.97*eta$se.fit
    logsirQ975 <- logsir + 1.97*eta$se.fit
    
    tibble(
      predictions = list(tibble(
        .x,
        predSIRest = exp(logsir),
        predSIRQ025 = exp(logsirQ025),
        predSIRQ975 = exp(logsirQ975)
      )),
      autocor = autocor,
      fit = list(fit)
    )
    
  })

sir$predicted_sir

sir$predicted <-
  sir$predicted_sir %>%
  unnest(predictions) %>%
  left_join(dat$ukrgeo) %>%
  st_as_sf()

# Model SIR via spatial smoothing NB regression -------------------

sir$dosage_vs_sir <-
  sir$ready_for_model %>%
  group_by(sex) %>%
  group_modify(~{
    
    model_dat <- .x
    logexpected <- log(model_dat$expected)
    
    fit <-
      gam(
        formula =
          observed ~
          s(Y, X, bs = 'tp') + log2(dos) + offset(log(expected)),
        family = nb(link = 'log'),
        data = model_dat
      )
    
    autocor <- testSpatialAutocorrelation(
      simulationOutput = simulateResiduals(fittedModel = fit),
      x = fit$model$X, y = fit$model$Y
    )[['p.value']]
    
    eta <- predict(
      fit, type = 'link', se.fit = TRUE, unconditional = TRUE
    )
    
    # eta, the linear predictor, is on the log incidence scale,
    # subtract log(expected) to get to the log SIR scale
    logsir <- eta$fit - logexpected
    # SD(X-c) = SD(X), thus the standard error also applies
    # when we put the linear predictor on the log SIR scale by
    # subtracting log(expected)
    logsirQ025 <- logsir - 1.97*eta$se.fit
    logsirQ975 <- logsir + 1.97*eta$se.fit
    
    tibble(
      predictions = list(tibble(
        .x,
        predSIRest = exp(logsir),
        predSIRQ025 = exp(logsirQ025),
        predSIRQ975 = exp(logsirQ975)
      )),
      autocor = autocor,
      fit = list(fit)
    )
    
  })

sir$dosage_vs_sir

# Plot smooth SIRs ------------------------------------------------

fig$sirsmooth <- list()
fig$sirsmooth$data <-
  sir$predicted %>%
  # set regions with non significant SIRs to NA
  # so that they appear grey in the map
  mutate(
    sir_predicted =
      ifelse(predSIRQ975 < 1 | predSIRQ025 > 1, predSIRest, NA)
  ) %>% 
  # discretize the SIR for a discrete colors scale
  mutate(
    sir_predicted_discrete =
      cut(sir_predicted, breaks = cnst$sir_breaks, right = FALSE,
          include.lowest = TRUE)
  )

fig$sirsmooth$plot <-
  fig$sirsmooth$data %>%
  ggplot() +
  geom_sf(aes(fill = sir_predicted_discrete), size = 0.001) +
  geom_sf(data = dat$ukr_outline_geo, fill = NA, size = 0.02) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  facet_wrap(~sex) +
  scale_fill_manual(values = cnst$sir_colors, labels = cnst$sir_labels,
                    na.value = 'grey60') +
  labs(
    title = 'Thyroid cancer standardized incidence ratios Ukraine 2001 all ages',
    fill = 'SIR',
    y = NULL,
    x = NULL
  ) +
  guides(fill = guide_legend(reverse = TRUE)) +
  MyGGplotTheme(axis = '', axis_ticks = '') +
  theme(axis.text = element_blank())
fig$sirsmooth$plot

# Plot crude SIRs -------------------------------------------------

fig$sircrude <- list()
fig$sircrude$data <-
  sir$sir_by_region_sex %>%
  left_join(dat$ukrgeo) %>%
  st_as_sf() %>%
  # set regions with non significant SIRs to NA
  # so that they appear grey in the map
  mutate(
    # sir_predicted =
    #   ifelse(sir_q975 < 1 | sir_q025 > 1, sir_est, NA)
    sir_predicted = sir_est
  ) %>% 
  # discretize the SIR for a discrete colors scale
  mutate(
    sir_predicted_discrete =
      cut(sir_predicted, breaks = cnst$sir_breaks, right = FALSE,
          include.lowest = TRUE)
  )

fig$sircrude$plot <-
  fig$sircrude$data %>%
  ggplot() +
  geom_sf(aes(fill = sir_predicted_discrete), size = 0.001) +
  geom_sf(data = dat$ukr_outline_geo, fill = NA, size = 0.02) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  facet_wrap(~sex) +
  scale_fill_manual(values = cnst$sir_colors, labels = cnst$sir_labels,
                    na.value = 'grey60') +
  labs(
    title = 'Thyroid cancer standardized incidence ratios Ukraine 2001 all ages',
    fill = 'SIR',
    y = NULL,
    x = NULL
  ) +
  guides(fill = guide_legend(reverse = TRUE)) +
  MyGGplotTheme(axis = '', axis_ticks = '') +
  theme(axis.text = element_blank())
fig$sircrude$plot

# Plot dosage received --------------------------------------------

fig$dosage <- list()
fig$dosage$data <-
  sir$avg_dosage_by_region %>%
  left_join(dat$ukrgeo) %>%
  st_as_sf()
  
fig$dosage$plot <-
  fig$dosage$data %>%
  ggplot() +
  geom_sf(aes(fill = dos), size = 0.001) +
  geom_sf(data = dat$ukr_outline_geo, fill = NA, size = 0.02) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  facet_wrap(~sex) +
  scale_fill_viridis_c(
    trans = 'log2', breaks = c(20, 40, 80, 160, 320, 640, 1280)) +
  labs(
    title = 'District population average absorbed thyroid dose in 1986',
    fill = 'mGy',
    y = NULL,
    x = NULL
  ) +
  MyGGplotTheme(axis = '', axis_ticks = '') +
  theme(axis.text = element_blank())
fig$dosage$plot

# Correlation -----------------------------------------------------

cor <- list()

cor$lm <-
  lm(
    formula = predSIRest ~ log2dos*sex,
    data =
      sir$predicted %>%
      mutate(log2dos = log2(dos))
  )

library(multcomp)
glht(
  cor$lm,
  linfct = c(
  Male = "log2dos + log2dos:sexmale == 0",
  Female = "log2dos == 0"
  )
) %>% confint()

fig$lmdossir <-
  sir$predicted %>%
  group_by(sex) %>%
  ggplot(aes(y = predSIRest, x = dos)) +
  geom_hline(yintercept = 1, color = 'grey50', size = 1) +
  geom_point(size = 0.1) +
  geom_smooth(method = 'lm', color = 'darkgreen') +
  scale_x_continuous(trans = 'log2', n.breaks = 6) +
  scale_y_continuous(trans = 'log2', n.breaks = 10) +
  facet_wrap(~sex) +
  MyGGplotTheme(axis = '', grid = 'xy') +
  labs(
    x = 'District population average absorbed thyroid dose in 1986 [mGy]',
    y = 'District thyroid cancer SIR in 2001'
  )
fig$lmdossir

# Export ----------------------------------------------------------

ExportFigure(
  fig$sirsmooth$plot, path = paths$output$out, filename = '10-sirsmooth',
  device = 'pdf',
  width = config$figspec$width, scale = 1
)

ExportFigure(
  fig$sircrude$plot, path = paths$output$out, filename = '10-sircrude',
  device = 'pdf',
  width = config$figspec$width, scale = 1
)

ExportFigure(
  fig$dosage$plot, path = paths$output$out, filename = '10-dosage',
  device = 'pdf',
  width = config$figspec$width
)

ExportFigure(
  fig$lmdossir, path = paths$output$out, filename = '10-lmdossir',
  device = 'pdf',
  width = config$figspec$width
)
