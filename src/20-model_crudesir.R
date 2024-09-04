# Calculate crude regional SIR of thyroid cancer incidence across Ukraine

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
  modelinput = 'out/11-modelinput.rds'
)
paths$output <- list(
  crudesir_rds = 'out/20-crudesir.rds',
  out = 'out'
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

# Functions -------------------------------------------------------

# Confidence intervals around crude standardized incidence rate
# based on delta method
SIRCI <- function(observed, expected, quant = 0.025) {
  Q <- qnorm(quant)
  EF <- exp(Q/sqrt(observed))
  CI <- observed/expected*EF
  CI <- ifelse(observed==0|expected==0, NA, CI)
}

# Calculate crude SIRs --------------------------------------------

dat$crudesir_by_age_sex <-
  dat$modelinput$expected_observed_by_region_sex |>
  mutate(
    sir_crude_est = incidence_observed/incidence_expected,
    sir_crude_lo = SIRCI(incidence_observed, incidence_expected,
                         config$quantiles$lo),
    sir_crude_hi = SIRCI(incidence_observed, incidence_expected,
                         config$quantiles$hi)
  )

# Plot crude SIRs -------------------------------------------------

crudesir <- list()

crudesir$cases_and_zeros <-
  dat$crudesir_by_age_sex  |> 
  group_by(sex) |>
  summarise(
    cases = sum(incidence_observed, na.rm = TRUE),
    zeros = sum(incidence_observed==0, na.rm = TRUE)
  )

crudesir$data <-
  dat$crudesir_by_age_sex |>
  left_join(dat$maptemplates$ukrgeo) |>
  st_as_sf()

crudesir$plot <- list()
map(c('female', 'male', 'total'), ~{
  crudesir$plot[[.x]] <<-
    crudesir$data |>
    filter(sex == .x) |>
    ggplot() +
    geom_sf(data = dat$maptemplates$background) +
    geom_sf(aes(fill = sir_crude_est),
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
crudesir$plot$total

# Export ----------------------------------------------------------

saveRDS(crudesir, paths$output$crudesir_rds)

ExportFigure(
  crudesir$plot$total, path = paths$output$out,
  filename = '20-crudesir_total',
  device = 'svg',
  width = config$figspec$width, scale = 1
)

ExportFigure(
  crudesir$plot$male, path = paths$output$out,
  filename = '20-crudesir_male',
  device = 'svg',
  width = config$figspec$width, scale = 1
)

ExportFigure(
  crudesir$plot$female, path = paths$output$out,
  filename = '20-crudesir_female',
  device = 'svg',
  width = config$figspec$width, scale = 1
)
