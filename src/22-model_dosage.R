# Estimate exposure rate ratio of absorbed dosage

# Init ------------------------------------------------------------

library(yaml)
library(tidyverse)
library(sf)
library(spdep)
library(DHARMa)
library(INLA)
library(inlabru)

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
  nb_rds = 'out/22-nb.rds',
  nbcar_rds = 'out/22-nbcar.rds',
  nbcar_errpergybyagesex_csv = 'out/22-errpergybyagesex.csv',
  nbcar_errpergybyage_csv = 'out/22-errpergybyage.csv',
  out = 'out/'
)

# global configuration
config <- read_yaml(paths$input$config)

cnst <- list(
  nsim = 10000,
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

# Spatial regression data -----------------------------------------

spatialdat <- list()

# prepare data
spatialdat$data <-
  dat$modelinput$region_sex_age |>
  filter(sex != 'total') |>
  mutate(
    ageatexposure = age - 15,
    ageatexposuresq = ageatexposure^2,
    sqrtageatexposure = sqrt(age - 15),
    log1pageatexposure = log1p(ageatexposure),
    ageatexposurefac = cut(
      ageatexposure, breaks = c(
        0, 5, 10, 20, Inf
      ), right = FALSE, labels = c(
        '0-4', '5-9', '10-19', '20+'
      ),
      ageatexposurefac2 = as.factor(data$ageatexposure)
    ),
    sex = factor(sex, levels = c('male', 'female')),
    sexsumto0 = case_when(sex == 'male' ~ -0.5, sex == 'female' ~ 0.5),
    dosGy = dose/1000,
    logdosGy = sqrt(dose/1000),
    dosGysq = dosGy^2,
    superregion_id = as.factor(substr(region_id, 1, 2)),
    # for INLAs CAR model (graph) factors in increasing integer order
    region_id_orig = region_id,
    region_id = as.integer(as.factor(region_id))
  ) |>
  arrange(age, ageatexposure, sex,
          superregion_id,
          region_id)

# get neighborhood matrix
spatialdat$nb <- poly2nb(
  filter(dat$maptemplates$ukrgeo, region_id %in% dat$locations$region_id),
  queen = TRUE,
  snap = 1e-4 # 10 meters
)
spatialdat$W <- listw2mat(nb2listw(spatialdat$nb, style = 'B', zero.policy = TRUE))

#image(nbcar$W)
# total number of neighbors
sum(spatialdat$W)
# average number of neighbors per district
sum(spatialdat$W)/nrow(spatialdat$W)
# maximum number of neighbors per district
max(rowSums(spatialdat$W))
# number of districts with 0 neighbors
sum(rowSums(spatialdat$W)==0)

# NB --------------------------------------------------------------

nb <- within(spatialdat, {

  set.seed(1986)
    
  # dimensions
  d <- list(nsim = cnst$nsim, ndat = nrow(data))
  
  # SPECIFY MODEL
  
  likelihood <- bru_obs(
    family = 'nbinomial',
    data = data,
    formula =
      incidence_2001 ~
      intercept +
      alpha_aae + alpha_aaesq +
      alpha_sex +
      log1p(err*exp(logerrmod_aae+logerrmod_sex)),
    E = population_2001
  )
  fit <- bru(
    components = ~
      # thyroid cancer risk baseline
      intercept(1) +
      alpha_aae(ageatexposure, model = 'linear') +
      alpha_aaesq(ageatexposuresq, model = 'linear') +
      alpha_sex(sexsumto0, model = 'linear') +
      # radiation exposure excess relative risk
      err(dosGy, model = 'linear') +
      logerrmod_aae(ageatexposure, model = 'linear') +
      logerrmod_sex(sexsumto0, model = 'linear'),
    likelihood,
    options = list(
      control.compute = list(dic = TRUE, waic = TRUE, cpo = FALSE),
      control.inla = list(
        int.strategy = 'eb'
      ),
      verbose = TRUE,
      bru_verbose = 4
    )
  )
  
  # SIMULATE FROM MODEL
  
  # err per gy
  errpergy <- predict(
    fit,
    expand.grid(
      sexsumto0 = unique(data$sexsumto0),
      ageatexposure = unique(data$ageatexposure),
      dosGy = 1
    ),
    ~ (err)*exp(logerrmod_aae+logerrmod_sex),
    n.samples = d$nsim
  )
  
  # coefficients of interest
  coefs <- predict(
    fit,
    NULL,
    ~ c(
      # baseline thyroid cancer rate at age of exposure 0 per million PY
      intercept = exp(intercept_latent)*1e6,
      # female baseline thyroid cancer rate at age of exposure 0 per million PY
      intercept_female = exp(intercept_latent + alpha_sex_latent/2)*1e6,
      # male baseline thyroid cancer rate at age of exposure 0 per million PY
      intercept_male = exp(intercept_latent - alpha_sex_latent/2)*1e6,
      # ratio of female to male baseline thyroid cancer rates
      alpha_sex = exp(alpha_sex_latent),
      # linear coefficient of age at exposure effect on log baseline thyroid cancer rate
      alpha_aae = alpha_aae_latent,
      # quadratic coefficient of squared age at exposure effect on log baseline thyroid cancer rate
      alpha_aaesq = alpha_aaesq_latent,
      # excess relative risk per 1 Gy at age of exposure 0
      err = err_latent,
      # female excess relative risk per 1 Gy at age of exposure 0
      err_female = err_latent*exp(logerrmod_sex_latent/2),
      # male excess relative risk per 1 Gy at age of exposure 0
      err_male = err_latent*exp(-logerrmod_sex_latent/2),
      # ratio of male to female excess relative risk per 1 Gy at age of exposure 0
      errmod_sex = exp(logerrmod_sex_latent),
      # multiplicative change in err over 5 year increase in age at exposure
      errmod_age = exp(logerrmod_aae_latent*5)
    ),
    n.samples = d$nsim
  )
  
  # predictive samples
  theta <- within(list(), {
    # predicted mean response
    nb_mu = generate(
      fit, newdata = data, formula = ~exp(
        intercept + alpha_aae + alpha_aaesq + alpha_sex +
          log1p(err*exp(logerrmod_aae+logerrmod_sex))
      )*population_2001,
      n.samples = d$nsim
    )
    # negative binominal parameters
    nb_size = fit$summary.hyperpar[1,1]
    nb_var = nb_mu + nb_mu^2/nb_size
    nb_overdispersion = 1/nb_size
    nb_p = nb_mu/nb_var
    nb_r = (nb_mu^2)/(nb_var-nb_mu)
  })
  
  # baseline thyroid cancer risk by age and sex
  baseline <- predict(
    fit,
    unique(data[,c('age', 'ageatexposure', 'ageatexposuresq', 'sexsumto0')]),
    ~ exp(intercept + alpha_aae + alpha_aaesq + alpha_sex),
    n.samples = d$nsim
  )
  
  # err over levels of dosage and by sex, age at exposure
  err <- predict(
    fit,
    expand_grid(
      unique(data[,c('ageatexposure', 'sexsumto0')]),
      data.frame(
        dosGy = c(0.1, 0.2, 0.3, 0.4, 0.5, 1, 2)
      )
    ),
    ~ (err)*exp(logerrmod_aae+logerrmod_sex),
    n.samples = d$nsim
  )
  
  # RESIDUAL DIAGNOSTICS
  
  # simulated count responses
  simulated_counts <- matrix(NA, nrow = d$ndat, ncol = d$nsim)
  for (i in 1:d$ndat) {
    simulated_counts[i,] <-
      rnbinom(d$nsim, size = theta$nb_r[i,], prob = theta$nb_p[i,])
  }
  simulated_counts <- simulated_counts[,!apply(simulated_counts, 2, anyNA)]
  # residual object
  dharm <- createDHARMa(
    simulated_counts, data$incidence_2001,
    fittedPredictedResponse = NULL, integerResponse = TRUE
  )
  # aggregate residuals to regions
  dharm_agg <- recalculateResiduals(dharm, group = data$region_id)
  residual_tests <- list(
    test_dispersion = testDispersion(dharm),
    test_zeroinflation = testZeroInflation(dharm),
    test_spatialautocor =
      testSpatialAutocorrelation(dharm_agg, dat$locations$X, dat$locations$Y)
  )
  
})

nb$textsummary <-
  list(
    date = date(),
    summary = summary(nb$fit),
    coefficients = round(nb$coefs, 3),
    overdispersion = nb$theta$nb_overdispersion,
    tests = nb$residual_tests
  )

capture.output(
  nb$textsummary, file = 'out/22-nb_textsummary.txt'
)

# NB-CAR ----------------------------------------------------------

nbcar <- within(spatialdat, {
  
  set.seed(1986)
  
  # dimensions
  d <- list(nsim = cnst$nsim, ndat = nrow(data))
  
  # SPECIFY MODEL
  
  likelihood <- bru_obs(
    family = 'nbinomial',
    data = data,
    formula =
      incidence_2001 ~
      intercept +
      alpha_aae + alpha_aaesq +
      alpha_sex +
      log1p(err*exp(logerrmod_aae+logerrmod_sex)) +
      zeta,
    E = population_2001
  )
  fit <- bru(
    components = ~
      # thyroid cancer risk baseline
      intercept(1) +
      alpha_aae(ageatexposure, model = 'linear') +
      alpha_aaesq(ageatexposuresq, model = 'linear') +
      alpha_sex(sexsumto0, model = 'linear') +
      # radiation exposure excess relative risk
      err(dosGy, model = 'linear') +
      logerrmod_aae(ageatexposure, model = 'linear') +
      logerrmod_sex(sexsumto0, model = 'linear') +
      # region control
      zeta(region_id, model = 'besag', graph = W),
    likelihood,
    options = list(
      control.compute = list(dic = TRUE, waic = TRUE, cpo = FALSE),
      control.inla = list(
        int.strategy = 'eb'
      ),
      verbose = TRUE,
      bru_verbose = 4
    )
  )
  
  # SIMULATE FROM MODEL
  
  # err per gy by age and sex
  errpergybyagesex <- predict(
    fit,
    expand.grid(
      sexsumto0 = unique(data$sexsumto0),
      ageatexposure = unique(data$ageatexposure),
      dosGy = 1
    ),
    ~ err*exp(logerrmod_aae+logerrmod_sex),
    n.samples = d$nsim
  )
  
  # err per gy by age
  errpergyage <- predict(
    fit,
    expand.grid(
      ageatexposure = unique(data$ageatexposure),
      dosGy = 1
    ),
    ~ err*exp(logerrmod_aae),
    n.samples = d$nsim
  )
  
  # coefficients of interest
  coefs <- predict(
    fit,
    NULL,
    ~ c(
      # baseline thyroid cancer rate at age of exposure 0 per million PY
      intercept = exp(intercept_latent)*1e6,
      # female baseline thyroid cancer rate at age of exposure 0 per million PY
      intercept_female = exp(intercept_latent + alpha_sex_latent/2)*1e6,
      # male baseline thyroid cancer rate at age of exposure 0 per million PY
      intercept_male = exp(intercept_latent - alpha_sex_latent/2)*1e6,
      # ratio of female to male baseline thyroid cancer rates
      alpha_sex = exp(alpha_sex_latent),
      # linear coefficient of age at exposure effect on log baseline thyroid cancer rate
      alpha_aae = alpha_aae_latent,
      # quadratic coefficient of squared age at exposure effect on log baseline thyroid cancer rate
      alpha_aaesq = alpha_aaesq_latent,
      # excess relative risk per 1 Gy at age of exposure 0
      err = err_latent,
      # female excess relative risk per 1 Gy at age of exposure 0
      err_female = err_latent*exp(logerrmod_sex_latent/2),
      # male excess relative risk per 1 Gy at age of exposure 0
      err_male = err_latent*exp(-logerrmod_sex_latent/2),
      # ratio of male to female excess relative risk per 1 Gy at age of exposure 0
      errmod_sex = exp(logerrmod_sex_latent),
      # multiplicative change in err over 5 year increase in age at exposure
      errmod_age = exp(logerrmod_aae_latent*5)
    ),
    n.samples = d$nsim
  )
  
  # predictive samples
  theta <- within(list(), {
    # predicted mean response
    nb_mu = generate(
      fit, newdata = data, formula = ~exp(
        intercept + alpha_aae + alpha_aaesq + alpha_sex +
          log1p(err*exp(logerrmod_aae+logerrmod_sex)) + zeta)*population_2001,
      n.samples = d$nsim
    )
    # negative binominal parameters
    nb_size = fit$summary.hyperpar[1,1]
    nb_var = nb_mu + nb_mu^2/nb_size
    nb_overdispersion = 1/nb_size
    nb_p = nb_mu/nb_var
    nb_r = (nb_mu^2)/(nb_var-nb_mu)
  })
  
  # baseline thyroid cancer risk by age and sex
  baseline <- predict(
    fit,
    unique(data[,c('age', 'ageatexposure', 'ageatexposuresq', 'sexsumto0')]),
    ~ exp(intercept + alpha_aae + alpha_aaesq + alpha_sex),
    n.samples = d$nsim
  )
  
  # err over levels of dosage and by sex, age at exposure
  err <- predict(
    fit,
    expand_grid(
      unique(data[,c('ageatexposure', 'sexsumto0')]),
      data.frame(
        dosGy = c(0.1, 0.2, 0.3, 0.4, 0.5, 1, 2)
      )
    ),
    ~ err*exp(logerrmod_aae+logerrmod_sex),
    n.samples = d$nsim
  )
  
  # region effects
  zeta <- predict(
    fit,
    unique(data[,c('region_id', 'region_id_orig')]),
    ~ exp(zeta),
    n.samples = d$nsim
  )
  
  # RESIDUAL DIAGNOSTICS
  
  # simulated count responses
  simulated_counts <- matrix(NA, nrow = d$ndat, ncol = d$nsim)
  for (i in 1:d$ndat) {
    simulated_counts[i,] <-
      rnbinom(d$nsim, size = theta$nb_r[i,], prob = theta$nb_p[i,])
  }
  simulated_counts <- simulated_counts[,!apply(simulated_counts, 2, anyNA)]
  # residual object
  dharm <- createDHARMa(
    simulated_counts, data$incidence_2001,
    fittedPredictedResponse = NULL, integerResponse = TRUE
  )
  # aggregate residuals to regions
  dharm_agg <- recalculateResiduals(dharm, group = data$region_id)
  residual_tests <- list(
    test_dispersion = testDispersion(dharm),
    test_zeroinflation = testZeroInflation(dharm),
    test_spatialautocor =
      testSpatialAutocorrelation(dharm_agg, dat$locations$X, dat$locations$Y)
  )
  
})

nbcar$textsummary <-
  list(
    date = date(),
    summary = summary(nbcar$fit),
    coefficients = round(nbcar$coefs, 3),
    overdispersion = nbcar$theta$nb_overdispersion,
    tests = nbcar$residual_tests
  )

capture.output(
  nbcar$textsummary, file = 'out/22-nbcar_textsummary.txt'
)

# Plot baseline ---------------------------------------------------

nbcar$plot <- list()

# baseline risk ratio over age by sex
nbcar$plot$baseline <-
  nbcar$baseline |>
  mutate(sex = ifelse(sexsumto0 == -0.5, 'Male', 'Female')) |>
  ggplot() +
  aes(x = ageatexposure) +
  geom_ribbon(aes(ymin = q0.025, ymax = q0.975, fill = sex), alpha = 0.1) +
  geom_line(aes(y = mean, color = sex)) +
  geom_point(aes(y = mean, color = sex)) +
  scale_y_continuous(labels = scales::label_comma(scale = 1e6), expand = c(0,0),
                     breaks = seq(0, 100, 10)/1e6) +
  scale_x_continuous(
    breaks = unique(nbcar$data$ageatexposure),
    labels = c(
      '0-4', '5-9', '10-14', '15-19', '20-24', '25-29', '30-34', '35-39',
      '40-44', '45-49', '50-54', '55-59', '60-64', '65-69', '70-74'
    )
  ) +
  labs(x = 'Age group', y = 'Baseline thyroid cancer incidence per million person-years') +
  MyGGplotTheme()

nbcar$plot$baseline

# Plot ERR --------------------------------------------------------

nbcar$plot$errpergybyagesex <-
  nbcar$errpergybyagesex |>
  mutate(sex = ifelse(sexsumto0 == -0.5, 'Male', 'Female')) |>
  ggplot() +
  aes(x = ageatexposure, y = median, color = sex, fill = sex, group = sex) +
  geom_errorbar(aes(ymax = q0.975, ymin = q0.025),
                position = position_dodge(width = 3), alpha = 0.5) +
  geom_point(position = position_dodge(width = 3)) +
  scale_x_continuous(breaks = unique(nbcar$data$ageatexposure), labels = c(
    '0-4', '5-9', '10-14', '15-19', '20-24', '25-29', '30-34', '35-39',
    '40-44', '45-49', '50-54', '55-59', '60-64', '65-69', '70-74'
  )) +
  scale_y_continuous(breaks = seq(0, 50, 10), expand = c(0.005, 0.005)) +
  labs(y = 'Excess relative risk per Gy', x = 'Age at exposure') +
  MyGGplotTheme()

nbcar$plot$errpergybyagesex

# Plot spatial estimates ------------------------------------------

zeta <- left_join(dat$maptemplates$ukrgeo, nbcar$zeta,
                  by = c(region_id = 'region_id_orig'))

nbcar$plot$zeta <-
  zeta |>
  ggplot() +
  geom_sf(data = dat$maptemplates$background) +
  geom_sf(aes(fill = mean),
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
  scale_fill_distiller(type = 'div', trans = 'log10',
                       na.value = config$figspec$na_color,
                       limits = c(1/3, 3),
                       oob = scales::squish,
                       breaks = c(1/3, 0.5, 1, 2, 3),
                       labels = c('<1/3', '1/2', '1', '2/1', '>3/1')
  ) +
  labs(
    fill = 'ICAR spatial effect',
    y = NULL,
    x = NULL
  ) +
  coord_sf(expand = FALSE) +
  MyGGplotTheme(axis = '', axis_ticks = '', panel_border = TRUE) +
  theme(axis.text = element_blank())

nbcar$plot$zeta

# Export ----------------------------------------------------------

saveRDS(dosage, paths$output$dosage_rds)
saveRDS(nb, paths$output$nb_rds)
saveRDS(nbcar, paths$output$nbcar_rds)

nbcar$errpergyage |>
  dplyr::select(ageatexposure, mean, sd, q0.025, q0.5, q0.975) |>
  mutate(across(where('is.numeric'), ~round(.x, 3))) |>
  write.csv(paths$output$nbcar_errpergybyage_csv, row.names = FALSE)

nbcar$errpergybyagesex |>
  mutate(sex = ifelse(sexsumto0 == -0.5, 'Male', 'Female')) |>
  dplyr::select(sex, ageatexposure, mean, sd, q0.025, q0.5, q0.975) |>
  mutate(across(where('is.numeric'), ~round(.x, 3))) |>
  write.csv(paths$output$nbcar_errpergybyagesex_csv, row.names = FALSE)

ExportFigure(
  nbcar$plot$baseline,
  path = paths$output$out, filename = '22-nbcar_baseline',
  device = 'pdf',
  width = config$figspec$width, scale = 1.2
)

ExportFigure(
  nbcar$plot$errpergybyagesex,
  path = paths$output$out, filename = '22-nbcar_errpergybyagesex',
  device = 'pdf',
  width = config$figspec$width, scale = 1.2
)

ExportFigure(
  nbcar$plot$zeta,
  path = paths$output$out, filename = '22-nbcar_zeta',
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
