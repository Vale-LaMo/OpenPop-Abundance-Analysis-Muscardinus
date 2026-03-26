# densities

library(tidyverse)
library(readxl)

abundances <- readxl::read_xlsx("outputs/STAN_hierarchical_mod_inform.xlsx")
head(abundances)

extents <- data.frame(
  area_name = c("R3", "R2", "R1", "S3", "S2", "S1"),
  extent_ha = c(9.04, 12.09, 8.35, 12.27, 11.21, 9.40)
)

densities <- abundances |>
  left_join(extents, by = "area_name") |>
  mutate(
    dens_ha = mean / extent_ha,
    month = month(date),
    season = case_when(
      month %in% c(5, 6) ~ "spring",
      month %in% c(9, 10) ~ "autumn",
      .default = "other"
    )
  )
head(densities)

densities |>
  group_by(area_name) |>
  summarise(min_dens = min(dens_ha), max_dens = max(dens_ha))

densities |>
  group_by(area_name, season) |>
  summarise(min_dens = min(dens_ha), max_dens = max(dens_ha))

densities |>
  group_by(area_name) |>
  summarise(min_abund = min(mean), max_abund = max(mean))
