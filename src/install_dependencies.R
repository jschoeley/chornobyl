install.packages(
  'INLA',
  repos = c(getOption('repos'), INLA = 'https://inla.r-inla-download.org/R/stable'),
  dep = TRUE
)

pkgs <- c(
  'ggplot2',
  'geomtextpath',
  'rnaturalearth',
  'sf',
  'tidyverse',
  'yaml',
  'scales',
  'DHARMa',
  'mgcv',
  'mgcViz',
  'dplyr',
  'multcomp',
  'ggpubr',
  'ggrepel',
  'rgeoda',
  'spatialreg',
  'spdep',
  'inlabru',
  'spdep'
)

install.packages(pkgs, dependencies = TRUE)
