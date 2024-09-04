# Calculate smooth regional SIR of thyroid cancer incidence across Ukraine

# Init ------------------------------------------------------------

library(yaml)
library(tidyverse)
library(mgcv)
library(sf)
# detect spatial autocorrelation via residual dignostics
library(DHARMa)
library(mgcViz)

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
  smoothsir_rds = 'out/21-smoothsir.rds',
  out = 'out/'
)

# global configuration
config <- read_yaml(paths$input$config)

# global objects
source(paths$input$glob)

# list containers for analysis artifacts
dat <- list()

# Input data ------------------------------------------------------

dat$modelinput <- readRDS(paths$input$modelinput)
dat$maptemplates <- readRDS(paths$input$maptemplates)
# retrieve region centroids for regions to be used in modeling
dat$locations <-
  dat$modelinput$region_sex_age |>
  select(region_id, X, Y) |>
  filter(!duplicated(region_id))

# Model SIR via spatial smoothing NB regression -------------------

# smooth thyroid cancer SIR by region and sex
dat$predicted_sir <-
  dat$modelinput$region_sex |>
  group_by(sex) |>
  group_modify(~{
    
    model_dat <- .x
    logexpected <- log(model_dat$incidence_expected)
    
    fit <-
      gam(
        formula =
          incidence_observed ~
          s(Y, X, bs = 'tp') + offset(log(incidence_expected)),
        family = nb(link = 'log'),
        data = model_dat
      )
    
    # test for spatial autocorrelation
    simresid <- simulateResiduals(fittedModel = fit)
    autocor <- testSpatialAutocorrelation(
      simulationOutput = simresid,
      x = fit$model$X, y = fit$model$Y
    )[['p.value']]
    # test for zero inflation
    zeroinflation <- testZeroInflation(simresid)[['p.value']]
    
    eta <- predict(
      fit, type = 'link', se.fit = TRUE, unconditional = TRUE
    )
    
    # eta, the linear predictor, is on the log incidence scale,
    # subtract log(expected) to get to the log SIR scale
    logsir <- c(eta[['fit']] - logexpected)
    # SD(X-c) = SD(X), thus the standard error also applies
    # when we put the linear predictor on the log SIR scale by
    # subtracting log(expected)
    logsirQ025 <- c(logsir - 1.97*eta[['se.fit']])
    logsirQ975 <- c(logsir + 1.97*eta[['se.fit']])
    
    tibble(
      predictions = list(tibble(
        .x,
        sir_smooth_est = exp(logsir),
        sir_smooth_qlo = exp(logsirQ025),
        sir_smooth_qhi = exp(logsirQ975)
      )),
      autocor = autocor,
      zeroinflation = zeroinflation,
      fit = list(fit)
    )
    
  }) |>
  ungroup()

# Plot smooth SIR -------------------------------------------------

smoothsir <- list()

smoothsir$data <-
  dat$predicted_sir |>
  unnest(predictions) |>
  # set regions with non significant SIRs to NA
  # so that they appear grey in the map
  mutate(
    sir_significant =
      ifelse(sir_smooth_qhi < 1 | sir_smooth_qlo > 1, TRUE, FALSE),
    sir_smooth_sig =
      ifelse(sir_significant, sir_smooth_est, NA)
  ) |> 
  left_join(dat$maptemplates$ukrgeo) |>
  st_as_sf()

smoothsir$plot <- list()
map(c('female', 'male', 'total'), ~{
  smoothsir$plot[[.x]] <<-
    smoothsir$data |>
    filter(sex == .x) |>
    ggplot() +
    geom_sf(data = dat$maptemplates$background) +
    geom_sf(aes(fill = sir_smooth_sig),
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
      fill = 'SIR',
      y = NULL,
      x = NULL
    ) +
    coord_sf(expand = FALSE) +
    MyGGplotTheme(axis = '', axis_ticks = '', panel_border = TRUE) +
    theme(axis.text = element_blank())
})
smoothsir$plot$total

# Export ----------------------------------------------------------

saveRDS(smoothsir, paths$output$smoothsir_rds)

ExportFigure(
  smoothsir$plot$total, path = paths$output$out,
  filename = '21-smoothsir_total',
  device = 'svg',
  width = config$figspec$width, scale = 1
)

ExportFigure(
  smoothsir$plot$female, path = paths$output$out,
  filename = '21-smoothsir_female',
  device = 'svg',
  width = config$figspec$width, scale = 1
)

ExportFigure(
  smoothsir$plot$male, path = paths$output$out,
  filename = '21-smoothsir_male',
  device = 'svg',
  width = config$figspec$width, scale = 1
)
