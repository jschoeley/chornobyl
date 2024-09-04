# Estimate exposure rate ratio of absorbed dosage

# Init ------------------------------------------------------------

library(yaml)
library(tidyverse)
library(sf)
library(DHARMa)
library(INLA)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  glob = 'src/00-global_functions.R',
  config = 'cfg/config.yaml',
  maptemplates = 'out/10-maptemplates.rds',
  # thyroid cancer incidence 2001 and radiation exposure 1986
  modelinput = 'out/11-modelinput.rds'
)
paths$output <- list(
  dosage_rds = 'out/22-dosage.rds',
  car1_rds = 'out/22-car1.rds',
  car2_rds = 'out/22-car2.rds',
  nb1_rds = 'out/22-nb1.rds',
  nb2_rds = 'out/22-nb2.rds',
  out = 'out/'
)

# global configuration
config <- read_yaml(paths$input$config)

cnst <- list(
  nsim = 500,
  incidencescaler = 1e6,
  cilo = 0.025,
  cihi = 0.975
)

# global objects
source(paths$input$glob)

# list containers for analysis artifacts
dat <- list()

# Input data ------------------------------------------------------

dat$modelinput <- readRDS(paths$input$modelinput)
dat$maptemplates <- readRDS(paths$input$maptemplates)
dat$locations <-
  dat$modelinput$region_sex_age |>
  dplyr::select(region_id, X, Y) |>
  filter(
    !duplicated(region_id),
    unique(region_id %in% dat$modelinput$region_sex_age$region_id)
  )

# Plot dosage map -------------------------------------------------

dosage <- list()

dosage$national_average_dose <-
  dat$modelinput$region_sex |>
  summarise(
    dose = sum(population_2001*average_dose)/sum(population_2001)
  )

dosage$data <-
  dat$modelinput$avg_dosage_by_region_sex |>
  left_join(dat$maptemplates$ukrgeo, by = 'region_id') |>
  st_as_sf()

dosage$plot <- list()
map(c('female', 'male', 'total'), ~{
  dosage$plot[[.x]] <<-
    dosage$data |>
    filter(sex == .x) |>
    ggplot() +
    geom_sf(data = dat$maptemplates$background) +
    geom_sf(aes(fill = average_dose),
            linewidth = config$figspec$district_outline_width) +
    geom_sf(data = dat$maptemplates$outline, fill = NA,
            linewidth = config$figspec$national_outline_width) +
    geom_sf(
      data = dat$maptemplates$cities,
      size = config$figspec$cities_point_size, shape = 1
    ) +
    geom_sf_text(
      data = dat$maptemplates$cities,
      aes(label = city),
      family = 'roboto',
      size = config$figspec$cities_text_size,
      hjust = 0, vjust = 0,  position = position_nudge(0.21, -0.21),
      color = 'white'
    ) +
    geom_sf_text(
      data = dat$maptemplates$cities,
      aes(label = city),
      family = 'roboto',
      size = config$figspec$cities_text_size,
      hjust = 0, vjust = 0, position = position_nudge(0.20, -0.20)
    ) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    scale_fill_distiller(
      type = 'div', trans = 'log2',
      na.value = config$figspec$na_color,
      limits = c(
        0.1*dosage$national_average_dose$dose,
        10*dosage$national_average_dose$dose),
      oob = scales::squish,
      breaks = c(0.1*dosage$national_average_dose$dose,
                 0.2*dosage$national_average_dose$dose,
                 0.5*dosage$national_average_dose$dose,
                 dosage$national_average_dose$dose,
                 2*dosage$national_average_dose$dose,
                 5*dosage$national_average_dose$dose,
                 10*dosage$national_average_dose$dose),
      labels = c(
        '<x1/10',
        'x1/5',
        'x1/2',
        paste0(round(dosage$national_average_dose$dose, 2),
               'mGy (National avg.)'),
        'x2',
        'x5',
        '>x10'
      )
    ) +
    labs(
      fill = 'Relative absorbed thyroid dose',
      y = NULL,
      x = NULL
    ) +
    coord_sf(expand = FALSE) +
    MyGGplotTheme(axis = '', axis_ticks = '', panel_border = TRUE) +
    theme(axis.text = element_blank())
})
dosage$plot$total

# OLS dosage -> incidence rate ------------------------------------

ols <- list()

ols$dat <-
  dat$modelinput$national_incidence |>
  group_by(sex) |>
  mutate(stdpop = national_population_2001 / sum(national_population_2001)) |>
  right_join(dat$modelinput$region_sex_age) |>
  filter(sex != 'total') |>
  mutate(
    rate = incidence_2001 / population_2001,
    rate_std = rate * stdpop
  ) |>
  group_by(region_id, sex) |>
  summarise(
    rate_std = sum(rate_std)
  ) |>
  ungroup() |>
  left_join(dat$modelinput$avg_dosage_by_region_sex) |>
  mutate(sex = factor(sex, c('male', 'female')),
         dosGy = average_dose/1000)

ols$lm <-
  glm(
    formula = rate_std ~ dosGy*sex,
    data = ols$dat,
    family = poisson(link = 'identity')
  )
summary(ols$lm)

library(multcomp)
glht(
  ols$lm,
  linfct = c(
    Male = "dosGy + dosGy:sexfemale == 0",
    Female = "dosGy == 0"
  )
) |> confint()

ols$fig <-
  ols$dat |>
  ggplot(aes(y = rate_std*1e6, x = dosGy)) +
  geom_point(size = 0.1) +
  geom_smooth(method = 'lm', color = 'darkgreen', se = FALSE) +
  scale_y_continuous(trans = 'log2', breaks = unlist(map(c(0.01,0.1, 1, 10, 100), ~c(1,3,5)*.x))) +
  scale_x_continuous(trans = 'log2', breaks = unlist(map(c(0.01,0.1, 1), ~c(1,3,5)*.x))) +
  facet_wrap(~sex) +
  MyGGplotTheme(axis = 'xy', grid = 'xy') +
  labs(
    x = 'District population average absorbed thyroid dose in 1986 [Gy]',
    y = 'District age-standardized thyroid cancer incidence in 2001 [per million]'
  )
ols$fig

# Empirical thyroid cancer risk ratios over age -------------------

rrage <- list()

rrage$data <-
  dat$modelinput$region_sex_age |>
  group_by(sex, age) |>
  summarise(
    inc = sum(incidence_2001),
    pop = sum(population_2001)
  ) |>
  mutate(
    rate = inc/pop,
    rr = rate/rate[3]
  )

rrage$data |>
  ggplot(aes(x = age, y = rr, color = sex)) +
  geom_smooth(formula = y ~ poly(x, 2), method = 'lm', se = FALSE) +
  geom_point() +
  scale_y_log10() +
  MyGGplotTheme() +
  labs(
    y = 'Risk ratio vs. age 25', x = 'Age',
    title = 'Risk ratio of thyroid cancer incidence over age vs. age 25'
  )

# CAR -------------------------------------------------------------

car <- list()

# prepare data
car$data <-
  dat$modelinput$region_sex_age |>
  filter(sex != 'total') |>
  mutate(
    ageminus25 = age - 25,
    sex = factor(sex, levels = c('male', 'female')),
    dosGy = dose/1000,
    superregion_id = as.factor(substr(region_id, 1, 2)),
    # for INLAs CAR model (graph) factors in increasing interger order
    region_id = as.integer(as.factor(region_id))
  ) |>
  arrange(age, ageminus25, sex,
          superregion_id,
          region_id)

# get neighborhood matrix
car$W <- st_touches(
  filter(dat$maptemplates$ukrgeo,
         region_id %in% dat$locations$region_id),
  sparse = FALSE
)
#image(car$W)
diag(car$W) <- FALSE

# NB-CAR1 (overall) -----------------------------------------------

car1 <- within(list(), {
  
  # dimensions
  d <- list(nsim = cnst$nsim, ndat = nrow(car$data))
  
  # fit nb-CAR model
  fit <- inla(
    incidence_2001 ~
      1 + sex +
      poly(ageminus25, 2, raw = TRUE) +
      sex:poly(ageminus25, 2, raw = TRUE) +
      dosGy +
      #superregion_id +
      f(region_id, model = 'besag', graph = car$W) +
      offset(log(population_2001)),
    data = car$data,
    family = 'nbinomial',
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE)
  )
  summary(fit)
  
  # extract simulations of model parameters
  posterior_samples <- inla.posterior.sample(cnst$nsim, result = fit)
  theta <- within(list(), {
    # beta's
    b_intercept =
      inla.posterior.sample.eval(c('Intercept'), posterior_samples)
    b_female =
      inla.posterior.sample.eval(c('sexfemale'), posterior_samples)
    b_dosgy =
      inla.posterior.sample.eval(c('dosGy'), posterior_samples)
    b_age1 =
      inla.posterior.sample.eval(c('1'), posterior_samples)
    b_age2 =
      inla.posterior.sample.eval(c('2'), posterior_samples)
    b_femaleage1 =
      inla.posterior.sample.eval(c('sexfemale:1'), posterior_samples)
    b_femaleage2 =
      inla.posterior.sample.eval(c('sexfemale:2'), posterior_samples)
    # predicted mean response
    nb_mu =
      exp(inla.posterior.sample.eval(c('Predictor'), posterior_samples))
    # negative binominal parameters
    nb_size =
      fit$summary.hyperpar$mean[1]
    nb_var =
      nb_mu + nb_mu^2/nb_size
    nb_overdispersion =
      1/nb_size
    nb_p =
      nb_mu/nb_var
    nb_r =
      (nb_mu^2)/(nb_var-nb_mu)
  })
  
  # residual diagnostics
  # simulated count responses
  simulated_counts <- matrix(NA, nrow = d$ndat, ncol = d$nsim)
  for (i in 1:d$ndat) {
    simulated_counts[i,] <-
      rnbinom(d$nsim, size = theta$nb_r[i,], prob = theta$nb_p[i,])
  }
  # residual object
  dharm <- createDHARMa(
    simulated_counts, car$data$incidence_2001,
    fittedPredictedResponse = NULL, integerResponse = TRUE
  )
  # aggregate residuals to regions
  dharm_agg <- recalculateResiduals(dharm, group = car$data$region_id)
  residual_tests <- list(
    test_dispersion = testDispersion(dharm),
    test_zeroinflation = testZeroInflation(dharm),
    test_spatialautocor =
      testSpatialAutocorrelation(dharm_agg, dat$locations$X, dat$locations$Y)
  )
  
  # statistics of interest
  statistics_overall <-
    expand_grid(draw = 1:d$nsim) |>
    # expand beta's of interest
    mutate(
      b_intercept = theta$b_intercept,
      b_female = theta$b_female,
      b_age1 = theta$b_age1,
      b_age2 = theta$b_age2,
      b_femaleage1 = theta$b_femaleage1,
      b_femaleage2 = theta$b_femaleage2,
      b_dosgy_m = theta$b_dosgy
    ) |>
    mutate(
      intercept_m = b_intercept,
      intercept_f = b_intercept+b_female,
      # annual thyroid cancer incidence per million @ 0 dosage age 25
      baselineincidence_m = exp(intercept_m)*cnst$incidencescaler,
      baselineincidence_f = exp(intercept_f)*cnst$incidencescaler,
      # rate ratio of female to male baseline incidence
      baselineincidence_femalefactor = exp(b_female),
      # relative increase in thyroid cancer incidence for
      # 1 Gy increase in dosage
      exposureriskratio = exp(b_dosgy_m)
    ) |>
    # summarize mean + quantiles
    summarise(across(
      starts_with(c('intercept', 'baselineincidence', 'exposureriskratio', 'b_')),
      list(avg = mean, qlo = ~quantile(.x, cnst$cilo),
           qhi = ~quantile(.x, cnst$cihi))
    )) |>
    ungroup()
  
  # summarise statistics of interest
  statistics_by_dosage <-
    # evalutate parameters over levels of dosage
    expand_grid(dosgy = 0:5, draw = 1:d$nsim) |>
    # expand beta's of interest
    mutate(
      b_dosgy_m = rep(theta$b_dosgy, 6),
      b_female_dosgy = rep(theta$b_female_dosgy, 6)
    ) |>
    # derive measure of interest
    mutate(
      # risk ratio of thyroid cancer @ dosage compared to 0 dosage
      dosageriskratio = exp(b_dosgy_m*dosgy)
    ) |>
    group_by(dosgy) |>
    # summarize mean + quantiles
    summarise(across(
      starts_with(c('dosageriskratio')),
      list(avg = mean, qlo = ~quantile(.x, cnst$cilo),
           qhi = ~quantile(.x, cnst$cihi))
    )) |>
    ungroup()
  
  statistics_by_age <-
    # evalutate parameters over age
    expand_grid(age = 15:85, draw = 1:d$nsim) |>
    # expand beta's of interest
    mutate(
      b_age1 = rep(theta$b_age1, 71),
      b_age2 = rep(theta$b_age2, 71),
      b_femaleage1 = rep(theta$b_femaleage1, 71),
      b_femaleage2 = rep(theta$b_femaleage2, 71),
    ) |>
    # derive measure of interest
    mutate(
      # risk ratio of thyroid cancer incidence over age
      riskratio_m = exp(b_age1*(age-25) + b_age2*(age-25)^2),
      riskratio_f = exp((b_age1+b_femaleage1)*(age-25) +
                          (b_age2+b_femaleage2)*(age-25)^2)
    ) |>
    group_by(age) |>
    # summarize mean + quantiles
    summarise(across(
      starts_with(c('riskratio')),
      list(avg = mean, qlo = ~quantile(.x, cnst$cilo),
           qhi = ~quantile(.x, cnst$cihi))
    )) |>
    ungroup()
  
})

summary(car1$fit)
car1$statistics_overall |> t()
car1$statistics_by_age
car1$statistics_by_dosage
car1$theta$nb_overdispersion
car1$residual_tests

car1$statistics_by_age |>
  ggplot() +
  aes(x = age) +
  geom_ribbon(aes(ymin = riskratio_m_qlo, ymax = riskratio_m_qhi),
              fill = 'blue', alpha = 0.1) +
  geom_line(aes(y = riskratio_m_avg), color = 'blue') +
  geom_ribbon(aes(ymin = riskratio_f_qlo, ymax = riskratio_f_qhi),
              fill = 'red', alpha = 0.1) +
  geom_line(aes(y = riskratio_f_avg), color = 'red')

# NB-CAR2 (by sex) ------------------------------------------------

car2 <- within(list(), {
  
  # dimensions
  d <- list(nsim = cnst$nsim, ndat = nrow(car$data))
  
  # fit nb-CAR model
  fit <- inla(
    incidence_2001 ~
      1 + sex +
      poly(ageminus25, 2, raw = TRUE) +
      sex:poly(ageminus25, 2, raw = TRUE) +
      dosGy + sex:dosGy +
      f(region_id, model = 'besag', graph = car$W) +
      offset(log(population_2001)),
    data = car$data,
    family = 'nbinomial',
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE)
  )
  summary(fit)
  
  # extract simulations of model parameters
  posterior_samples <- inla.posterior.sample(cnst$nsim, result = fit)
  theta <- within(list(), {
    # beta's
    b_intercept =
      inla.posterior.sample.eval(c('Intercept'), posterior_samples)
    b_female =
      inla.posterior.sample.eval(c('sexfemale'), posterior_samples)
    b_dosgy =
      inla.posterior.sample.eval(c('dosGy'), posterior_samples)
    b_female_dosgy =
      inla.posterior.sample.eval(c('sexfemale:dosGy'), posterior_samples)
    b_age1 =
      inla.posterior.sample.eval(c('1'), posterior_samples)
    b_age2 =
      inla.posterior.sample.eval(c('2'), posterior_samples)
    b_femaleage1 =
      inla.posterior.sample.eval(c('sexfemale:1'), posterior_samples)
    b_femaleage2 =
      inla.posterior.sample.eval(c('sexfemale:2'), posterior_samples)
    # predicted mean response
    nb_mu =
      exp(inla.posterior.sample.eval(c('Predictor'), posterior_samples))
    # negative binominal parameters
    nb_size =
      fit$summary.hyperpar$mean[1]
    nb_var =
      nb_mu + nb_mu^2/nb_size
    nb_overdispersion =
      1/nb_size
    nb_p =
      nb_mu/nb_var
    nb_r =
      (nb_mu^2)/(nb_var-nb_mu)
  })
  
  # residual diagnostics
  # simulated count responses
  simulated_counts <- matrix(NA, nrow = d$ndat, ncol = d$nsim)
  for (i in 1:d$ndat) {
    simulated_counts[i,] <-
      rnbinom(d$nsim, size = theta$nb_r[i,], prob = theta$nb_p[i,])
  }
  # residual object
  dharm <- createDHARMa(
    simulated_counts, car$data$incidence_2001,
    fittedPredictedResponse = NULL, integerResponse = TRUE
  )
  # aggregate residuals to regions
  dharm_agg <- recalculateResiduals(dharm, group = car$data$region_id)
  residual_tests <- list(
    test_dispersion = testDispersion(dharm),
    test_zeroinflation = testZeroInflation(dharm),
    test_spatialautocor =
      testSpatialAutocorrelation(dharm_agg, dat$locations$X, dat$locations$Y)
  )
  
  # statistics of interest
  statistics_overall <-
    expand_grid(draw = 1:d$nsim) |>
    # expand beta's of interest
    mutate(
      b_intercept = theta$b_intercept,
      b_female = theta$b_female,
      b_age1 = theta$b_age1,
      b_age2 = theta$b_age2,
      b_femaleage1 = theta$b_femaleage1,
      b_femaleage2 = theta$b_femaleage2,
      b_dosgy_m = theta$b_dosgy,
      b_female_dosgy = theta$b_female_dosgy
    ) |>
    mutate(
      intercept_m = b_intercept,
      intercept_f = b_intercept+b_female,
      # annual thyroid cancer incidence per million @ 0 dosage age 25
      baselineincidence_m = exp(intercept_m)*cnst$incidencescaler,
      baselineincidence_f = exp(intercept_f)*cnst$incidencescaler,
      # rate ratio of female to male baseline incidence
      baselineincidence_femalefactor = exp(b_female),
      # relative increase in thyroid cancer incidence for
      # 1 Gy increase in dosage
      exposureriskratio_m = exp(b_dosgy_m),
      exposureriskratio_f = exp(b_dosgy_m+b_female_dosgy),
      # rate ratio of female to male exposure risk ratio
      exposureriskratio_femalefactor = exp(b_female_dosgy),
    ) |>
    # summarize mean + quantiles
    summarise(across(
      starts_with(c('intercept', 'baselineincidence', 'exposureriskratio', 'b_')),
      list(avg = mean, qlo = ~quantile(.x, cnst$cilo),
           qhi = ~quantile(.x, cnst$cihi))
    )) |>
    ungroup()
  
  # summarise statistics of interest
  statistics_by_dosage <-
    # evalutate parameters over levels of dosage
    expand_grid(dosgy = 0:5, draw = 1:d$nsim) |>
    # expand beta's of interest
    mutate(
      b_dosgy_m = rep(theta$b_dosgy, 6),
      b_female_dosgy = rep(theta$b_female_dosgy, 6)
    ) |>
    # derive measure of interest
    mutate(
      # risk ratio of thyroid cancer @ dosage compared to 0 dosage
      dosageriskratio_m = exp(b_dosgy_m*dosgy),
      dosageriskratio_f = exp((b_dosgy_m+b_female_dosgy)*dosgy)
    ) |>
    group_by(dosgy) |>
    # summarize mean + quantiles
    summarise(across(
      starts_with(c('dosageriskratio')),
      list(avg = mean, qlo = ~quantile(.x, cnst$cilo),
           qhi = ~quantile(.x, cnst$cihi))
    )) |>
    ungroup()
  
  statistics_by_age <-
    # evalutate parameters over age
    expand_grid(age = 15:85, draw = 1:d$nsim) |>
    # expand beta's of interest
    mutate(
      b_age1 = rep(theta$b_age1, 71),
      b_age2 = rep(theta$b_age2, 71),
      b_femaleage1 = rep(theta$b_femaleage1, 71),
      b_femaleage2 = rep(theta$b_femaleage2, 71),
    ) |>
    # derive measure of interest
    mutate(
      # risk ratio of thyroid cancer incidence over age
      riskratio_m = exp(b_age1*(age-25) + b_age2*(age-25)^2),
      riskratio_f = exp((b_age1+b_femaleage1)*(age-25) +
                          (b_age2+b_femaleage2)*(age-25)^2)
    ) |>
    group_by(age) |>
    # summarize mean + quantiles
    summarise(across(
      starts_with(c('riskratio')),
      list(avg = mean, qlo = ~quantile(.x, cnst$cilo),
           qhi = ~quantile(.x, cnst$cihi))
    )) |>
    ungroup()
  
})

summary(car2$fit)
car2$statistics_overall |> t()
car2$statistics_by_age
car2$statistics_by_dosage
car2$theta$nb_overdispersion
car2$residual_tests

car2$statistics_by_age |>
  ggplot() +
  aes(x = age) +
  geom_ribbon(aes(ymin = riskratio_m_qlo, ymax = riskratio_m_qhi),
              fill = 'blue', alpha = 0.1) +
  geom_line(aes(y = riskratio_m_avg), color = 'blue') +
  geom_ribbon(aes(ymin = riskratio_f_qlo, ymax = riskratio_f_qhi),
              fill = 'red', alpha = 0.1) +
  geom_line(aes(y = riskratio_f_avg), color = 'red')

# NB --------------------------------------------------------------

nb <- list()

# prepare data
nb$data <-
  dat$modelinput$region_sex_age |>
  filter(sex != 'total') |>
  mutate(
    sex = factor(sex, levels = c('male', 'female')),
    dosGy = dose/1000
  )

# NB1 (overall) ---------------------------------------------------

nb1 <- within(list(), {
  
  # dimensions
  d <- list(nsim = cnst$nsim, ndat = nrow(car$data))
  
  # fit nb-CAR model
  fit <- inla(
    incidence_2001 ~
      1 + sex +
      poly(ageminus25, 2, raw = TRUE) +
      sex:poly(ageminus25, 2, raw = TRUE) +
      dosGy +
      offset(log(population_2001)),
    data = car$data,
    family = 'nbinomial',
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE)
  )
  summary(fit)
  
  # extract simulations of model parameters
  posterior_samples <- inla.posterior.sample(cnst$nsim, result = fit)
  theta <- within(list(), {
    # beta's
    b_intercept =
      inla.posterior.sample.eval(c('Intercept'), posterior_samples)
    b_female =
      inla.posterior.sample.eval(c('sexfemale'), posterior_samples)
    b_dosgy =
      inla.posterior.sample.eval(c('dosGy'), posterior_samples)
    b_age1 =
      inla.posterior.sample.eval(c('1'), posterior_samples)
    b_age2 =
      inla.posterior.sample.eval(c('2'), posterior_samples)
    b_femaleage1 =
      inla.posterior.sample.eval(c('sexfemale:1'), posterior_samples)
    b_femaleage2 =
      inla.posterior.sample.eval(c('sexfemale:2'), posterior_samples)
    # predicted mean response
    nb_mu =
      exp(inla.posterior.sample.eval(c('Predictor'), posterior_samples))
    # negative binominal parameters
    nb_size =
      fit$summary.hyperpar$mean[1]
    nb_var =
      nb_mu + nb_mu^2/nb_size
    nb_overdispersion =
      1/nb_size
    nb_p =
      nb_mu/nb_var
    nb_r =
      (nb_mu^2)/(nb_var-nb_mu)
  })
  
  # residual diagnostics
  # simulated count responses
  simulated_counts <- matrix(NA, nrow = d$ndat, ncol = d$nsim)
  for (i in 1:d$ndat) {
    simulated_counts[i,] <-
      rnbinom(d$nsim, size = theta$nb_r[i,], prob = theta$nb_p[i,])
  }
  # residual object
  dharm <- createDHARMa(
    simulated_counts, car$data$incidence_2001,
    fittedPredictedResponse = NULL, integerResponse = TRUE
  )
  # aggregate residuals to regions
  dharm_agg <- recalculateResiduals(dharm, group = car$data$region_id)
  residual_tests <- list(
    test_dispersion = testDispersion(dharm),
    test_zeroinflation = testZeroInflation(dharm),
    test_spatialautocor =
      testSpatialAutocorrelation(dharm_agg, dat$locations$X, dat$locations$Y)
  )
  
  # statistics of interest
  statistics_overall <-
    expand_grid(draw = 1:d$nsim) |>
    # expand beta's of interest
    mutate(
      b_intercept = theta$b_intercept,
      b_female = theta$b_female,
      b_age1 = theta$b_age1,
      b_age2 = theta$b_age2,
      b_femaleage1 = theta$b_femaleage1,
      b_femaleage2 = theta$b_femaleage2,
      b_dosgy_m = theta$b_dosgy
    ) |>
    mutate(
      intercept_m = b_intercept,
      intercept_f = b_intercept+b_female,
      # annual thyroid cancer incidence per million @ 0 dosage age 25
      baselineincidence_m = exp(intercept_m)*cnst$incidencescaler,
      baselineincidence_f = exp(intercept_f)*cnst$incidencescaler,
      # rate ratio of female to male baseline incidence
      baselineincidence_femalefactor = exp(b_female),
      # relative increase in thyroid cancer incidence for
      # 1 Gy increase in dosage
      exposureriskratio = exp(b_dosgy_m)
    ) |>
    # summarize mean + quantiles
    summarise(across(
      starts_with(c('intercept', 'baselineincidence', 'exposureriskratio', 'b_')),
      list(avg = mean, qlo = ~quantile(.x, cnst$cilo),
           qhi = ~quantile(.x, cnst$cihi))
    )) |>
    ungroup()
  
  # summarise statistics of interest
  statistics_by_dosage <-
    # evalutate parameters over levels of dosage
    expand_grid(dosgy = 0:5, draw = 1:d$nsim) |>
    # expand beta's of interest
    mutate(
      b_dosgy_m = rep(theta$b_dosgy, 6),
      b_female_dosgy = rep(theta$b_female_dosgy, 6)
    ) |>
    # derive measure of interest
    mutate(
      # risk ratio of thyroid cancer @ dosage compared to 0 dosage
      dosageriskratio = exp(b_dosgy_m*dosgy)
    ) |>
    group_by(dosgy) |>
    # summarize mean + quantiles
    summarise(across(
      starts_with(c('dosageriskratio')),
      list(avg = mean, qlo = ~quantile(.x, cnst$cilo),
           qhi = ~quantile(.x, cnst$cihi))
    )) |>
    ungroup()
  
  statistics_by_age <-
    # evalutate parameters over age
    expand_grid(age = 15:85, draw = 1:d$nsim) |>
    # expand beta's of interest
    mutate(
      b_age1 = rep(theta$b_age1, 71),
      b_age2 = rep(theta$b_age2, 71),
      b_femaleage1 = rep(theta$b_femaleage1, 71),
      b_femaleage2 = rep(theta$b_femaleage2, 71),
    ) |>
    # derive measure of interest
    mutate(
      # risk ratio of thyroid cancer incidence over age
      riskratio_m = exp(b_age1*(age-25) + b_age2*(age-25)^2),
      riskratio_f = exp((b_age1+b_femaleage1)*(age-25) +
                          (b_age2+b_femaleage2)*(age-25)^2)
    ) |>
    group_by(age) |>
    # summarize mean + quantiles
    summarise(across(
      starts_with(c('riskratio')),
      list(avg = mean, qlo = ~quantile(.x, cnst$cilo),
           qhi = ~quantile(.x, cnst$cihi))
    )) |>
    ungroup()
  
})

summary(nb1$fit)
nb1$statistics_overall |> t() |> formatC(digits = 2, format = 'f')
nb1$statistics_by_age
nb1$statistics_by_dosage
nb1$theta$nb_overdispersion
nb1$residual_tests

nb1$statistics_by_age |>
  ggplot() +
  aes(x = age) +
  geom_ribbon(aes(ymin = riskratio_m_qlo, ymax = riskratio_m_qhi),
              fill = 'blue', alpha = 0.1) +
  geom_line(aes(y = riskratio_m_avg), color = 'blue') +
  geom_ribbon(aes(ymin = riskratio_f_qlo, ymax = riskratio_f_qhi),
              fill = 'red', alpha = 0.1) +
  geom_line(aes(y = riskratio_f_avg), color = 'red')

# NB2 (by sex) ----------------------------------------------------

nb2 <- within(list(), {
  
  # dimensions
  d <- list(nsim = cnst$nsim, ndat = nrow(car$data))
  
  # fit nb-CAR model
  fit <- inla(
    incidence_2001 ~
      1 + sex +
      poly(ageminus25, 2, raw = TRUE) +
      sex:poly(ageminus25, 2, raw = TRUE) +
      dosGy + sex:dosGy +
      offset(log(population_2001)),
    data = car$data,
    family = 'nbinomial',
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE)
  )
  summary(fit)
  
  # extract simulations of model parameters
  posterior_samples <- inla.posterior.sample(cnst$nsim, result = fit)
  theta <- within(list(), {
    # beta's
    b_intercept =
      inla.posterior.sample.eval(c('Intercept'), posterior_samples)
    b_female =
      inla.posterior.sample.eval(c('sexfemale'), posterior_samples)
    b_dosgy =
      inla.posterior.sample.eval(c('dosGy'), posterior_samples)
    b_female_dosgy =
      inla.posterior.sample.eval(c('sexfemale:dosGy'), posterior_samples)
    b_age1 =
      inla.posterior.sample.eval(c('1'), posterior_samples)
    b_age2 =
      inla.posterior.sample.eval(c('2'), posterior_samples)
    b_femaleage1 =
      inla.posterior.sample.eval(c('sexfemale:1'), posterior_samples)
    b_femaleage2 =
      inla.posterior.sample.eval(c('sexfemale:2'), posterior_samples)
    # predicted mean response
    nb_mu =
      exp(inla.posterior.sample.eval(c('Predictor'), posterior_samples))
    # negative binominal parameters
    nb_size =
      fit$summary.hyperpar$mean[1]
    nb_var =
      nb_mu + nb_mu^2/nb_size
    nb_overdispersion =
      1/nb_size
    nb_p =
      nb_mu/nb_var
    nb_r =
      (nb_mu^2)/(nb_var-nb_mu)
  })
  
  # residual diagnostics
  # simulated count responses
  simulated_counts <- matrix(NA, nrow = d$ndat, ncol = d$nsim)
  for (i in 1:d$ndat) {
    simulated_counts[i,] <-
      rnbinom(d$nsim, size = theta$nb_r[i,], prob = theta$nb_p[i,])
  }
  # residual object
  dharm <- createDHARMa(
    simulated_counts, car$data$incidence_2001,
    fittedPredictedResponse = NULL, integerResponse = TRUE
  )
  # aggregate residuals to regions
  dharm_agg <- recalculateResiduals(dharm, group = car$data$region_id)
  residual_tests <- list(
    test_dispersion = testDispersion(dharm),
    test_zeroinflation = testZeroInflation(dharm),
    test_spatialautocor =
      testSpatialAutocorrelation(dharm_agg, dat$locations$X, dat$locations$Y)
  )
  
  # statistics of interest
  statistics_overall <-
    expand_grid(draw = 1:d$nsim) |>
    # expand beta's of interest
    mutate(
      b_intercept = theta$b_intercept,
      b_female = theta$b_female,
      b_age1 = theta$b_age1,
      b_age2 = theta$b_age2,
      b_femaleage1 = theta$b_femaleage1,
      b_femaleage2 = theta$b_femaleage2,
      b_dosgy_m = theta$b_dosgy,
      b_female_dosgy = theta$b_female_dosgy
    ) |>
    mutate(
      intercept_m = b_intercept,
      intercept_f = b_intercept+b_female,
      # annual thyroid cancer incidence per million @ 0 dosage age 25
      baselineincidence_m = exp(intercept_m)*cnst$incidencescaler,
      baselineincidence_f = exp(intercept_f)*cnst$incidencescaler,
      # rate ratio of female to male baseline incidence
      baselineincidence_femalefactor = exp(b_female),
      # relative increase in thyroid cancer incidence for
      # 1 Gy increase in dosage
      exposureriskratio_m = exp(b_dosgy_m),
      exposureriskratio_f = exp(b_dosgy_m+b_female_dosgy),
      # rate ratio of female to male exposure risk ratio
      exposureriskratio_femalefactor = exp(b_female_dosgy),
    ) |>
    # summarize mean + quantiles
    summarise(across(
      starts_with(c('intercept', 'baselineincidence', 'exposureriskratio', 'b_')),
      list(avg = mean, qlo = ~quantile(.x, cnst$cilo),
           qhi = ~quantile(.x, cnst$cihi))
    )) |>
    ungroup()
  
  # summarise statistics of interest
  statistics_by_dosage <-
    # evalutate parameters over levels of dosage
    expand_grid(dosgy = 0:5, draw = 1:d$nsim) |>
    # expand beta's of interest
    mutate(
      b_dosgy_m = rep(theta$b_dosgy, 6),
      b_female_dosgy = rep(theta$b_female_dosgy, 6)
    ) |>
    # derive measure of interest
    mutate(
      # risk ratio of thyroid cancer @ dosage compared to 0 dosage
      dosageriskratio_m = exp(b_dosgy_m*dosgy),
      dosageriskratio_f = exp((b_dosgy_m+b_female_dosgy)*dosgy)
    ) |>
    group_by(dosgy) |>
    # summarize mean + quantiles
    summarise(across(
      starts_with(c('dosageriskratio')),
      list(avg = mean, qlo = ~quantile(.x, cnst$cilo),
           qhi = ~quantile(.x, cnst$cihi))
    )) |>
    ungroup()
  
  statistics_by_age <-
    # evalutate parameters over age
    expand_grid(age = 15:85, draw = 1:d$nsim) |>
    # expand beta's of interest
    mutate(
      b_age1 = rep(theta$b_age1, 71),
      b_age2 = rep(theta$b_age2, 71),
      b_femaleage1 = rep(theta$b_femaleage1, 71),
      b_femaleage2 = rep(theta$b_femaleage2, 71),
    ) |>
    # derive measure of interest
    mutate(
      # risk ratio of thyroid cancer incidence over age
      riskratio_m = exp(b_age1*(age-25) + b_age2*(age-25)^2),
      riskratio_f = exp((b_age1+b_femaleage1)*(age-25) +
                          (b_age2+b_femaleage2)*(age-25)^2)
    ) |>
    group_by(age) |>
    # summarize mean + quantiles
    summarise(across(
      starts_with(c('riskratio')),
      list(avg = mean, qlo = ~quantile(.x, cnst$cilo),
           qhi = ~quantile(.x, cnst$cihi))
    )) |>
    ungroup()
  
})

summary(nb2$fit)
nb2$statistics_overall |> t() |> formatC(digits = 2, format = 'f')
nb2$statistics_by_age
nb2$statistics_by_dosage
nb2$theta$nb_overdispersion
nb2$residual_tests

nb2$statistics_by_age |>
  ggplot() +
  aes(x = age) +
  geom_ribbon(aes(ymin = riskratio_m_qlo, ymax = riskratio_m_qhi),
              fill = 'blue', alpha = 0.1) +
  geom_line(aes(y = riskratio_m_avg), color = 'blue') +
  geom_ribbon(aes(ymin = riskratio_f_qlo, ymax = riskratio_f_qhi),
              fill = 'red', alpha = 0.1) +
  geom_line(aes(y = riskratio_f_avg), color = 'red')

# Plot dosage-risk-ratio ------------------------------------------

car2$dosageriskratio <-
  car2$statistics_by_dosage |>
  ggplot(aes(x = dosgy)) +
  geom_hline(yintercept = 1, color = 'black') +
  geom_ribbon(
    aes(
      ymin = dosageriskratio_qlo,
      ymax = dosageriskratio_qhi
    ),
    fill = NA, alpha = 0.1, color = 'black',
    lty = 2,
    data = car1$statistics_by_dosage
  ) +
  geom_ribbon(
    aes(
      ymin = dosageriskratio_m_qlo,
      ymax = dosageriskratio_m_qhi
    ),
    color = '#4295f5', alpha = 0.1, fill = NA,
    lty = 2
  ) +
  geom_ribbon(
    aes(
      ymin = dosageriskratio_f_qlo,
      ymax = dosageriskratio_f_qhi
    ),
    color = '#f5425a', alpha = 0.1, fill = NA,
    lty = 2
  ) +
  geom_line(
    aes(y = dosageriskratio_avg),
    color = 'black', size = 1,
    data = car1$statistics_by_dosage
  ) +
  geom_line(
    aes(y = dosageriskratio_m_avg),
    color = '#4295f5', size = 1
  ) +
  geom_line(
    aes(y = dosageriskratio_f_avg),
    color = '#f5425a', size = 1
  ) +
  scale_y_continuous(
    trans = 'log',
    breaks = c(0.1, 0.5, unlist(map(c(1,10, 100, 1000),
                                    ~c(1,2,3,4,5)*.x))),
  ) +
  labs(x = 'District population average absorbed thyroid dose in 1986 [Gy]', y = 'Exposure risk ratio') +
  coord_cartesian(xlim = c(0, 1.5), ylim = c(0.5, 10),
                  expand = FALSE) +
  MyGGplotTheme(axis = 'y', grid = 'xy')
car2$dosageriskratio

dat$modelinput$region_sex |>
  ggplot() +
  geom_histogram(
    aes(x = average_dose/1e3, weight = population_2001),
    breaks = c(0, 0.02, 0.04, 0.08, 0.1, 0.2, 0.5, 2)
  ) +
  coord_cartesian(xlim = c(0, 1.5),
                  expand = FALSE) +
  #scale_y_log10() +
  MyGGplotTheme(axis = 'xy', grid = 'xy')

# Export ----------------------------------------------------------

saveRDS(dosage, paths$output$dosage_rds)
saveRDS(car1, paths$output$car1_rds)
saveRDS(car2, paths$output$car2_rds)
saveRDS(nb1, paths$output$nb1_rds)
saveRDS(nb2, paths$output$nb2_rds)

ExportFigure(
  car2$dosageriskratio, path = paths$output$out, filename = '22-car2dosageriskratio',
  device = 'pdf',
  width = config$figspec$width, scale = 1.2
)

ExportFigure(
  ols$fig, path = paths$output$out, filename = '22-dosagevsincidence',
  device = 'pdf',
  width = config$figspec$width, scale = 1.2
)

ExportFigure(
  dosage$plot$total, path = paths$output$out,
  filename = '22-dosage_total',
  device = 'svg',
  width = config$figspec$width, scale = 1
)
