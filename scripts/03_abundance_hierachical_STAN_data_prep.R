## Data preparation for abundance estimation - hierarchical STAN model
{
  # load packages
  library(readxl)
  library(lubridate)
  library(cmdstanr)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(posterior)
}

{
  # function for time intervals (calculates difference in months)
  intervals_months <- function(dates) {
    diffs <- difftime(dates[-1], dates[-length(dates)], units = "days")
    as.numeric(diffs) / 30.44 # 1 month ≈ 30.44 days (mean Gregorian month)
  }
}

# every study area must share the same "master" timeline
# We need to align all areas to the same set of sampling dates.
# If an area was not sampled on a specific date, that column should be filled
# with zeros (or a placeholder if you distinguish "not sampled" from "zero captures").

# Define the Master List of unique sampling dates across ALL areas
master_dates <- sort(unique(as.Date(c(
  # "2019-07-15", "2019-08-05", "2019-08-20", "2019-10-10",
  # "2020-09-15",
  # "2021-09-15",
  "2022-05-05",
  "2022-05-20",
  "2022-06-05",
  "2022-06-20",
  "2022-07-05",
  "2022-07-20",
  "2022-08-05",
  "2022-08-20",
  "2022-09-05",
  "2022-09-20",
  "2022-10-15",
  #"2022-10-20", "2022-11-05", "2022-11-20", "2022-12-05",
  "2023-05-15",
  "2023-06-15",
  "2023-07-15",
  "2023-08-05",
  "2023-08-20",
  "2023-09-05",
  "2023-09-20", # "2023-10-10",
  "2024-06-15",
  "2024-07-15",
  "2024-08-15",
  "2024-09-15",
  "2024-10-15"
))))

# Calculate time intervals based on the Master Timeline
time_intervals <- intervals_months(master_dates)
K_total <- length(master_dates)


# List of study area names as they appear in your Excel sheets
areas <- c("S1", "S2", "S3", "R1", "R2", "R3")

# -----------------------------------------------------------
# file di riferimento "data/data_abundance.xlsx", sheet = areas[i])
# attenzione, deve esserci una colonna per ogni occasione per tutte le aree
# -----------------------------------------------------------

# Initialize containers
y_stacked <- matrix(nrow = 0, ncol = length(time_intervals) + 1)
area_idx_stacked <- c()

# Loop through each area to stack data
for (i in 1:length(areas)) {
  # 1. Load data
  temp_data <- readxl::read_xlsx(
    "data/data_abundance_art1.xlsx",
    sheet = areas[i]
  )
  temp_data[is.na(temp_data)] <- 0
  CH_raw <- as.matrix(temp_data[, -1]) # Remove non-capture columns (sex column)

  # 2. Augment per area
  n_obs <- nrow(CH_raw)
  M_area <- n_obs * 3 # Consistency in augmentation factor
  CH_aug <- rbind(CH_raw, matrix(0, nrow = M_area - n_obs, ncol = ncol(CH_raw)))

  # 3. Stack
  y_stacked <- rbind(y_stacked, CH_aug)
  area_idx_stacked <- c(area_idx_stacked, rep(i, M_area))
}

## ---- Data preparation for STAN ----
# Final Stan Data List
stan_data <- list(
  N_areas = length(areas), # Number of study area
  K = ncol(y_stacked), # Number of capture sessions
  M_total = nrow(y_stacked), # Total number of individuals in augmented dataset
  y = y_stacked, # Augmented dataset
  area_idx = area_idx_stacked, # ID number for each study area
  delta = time_intervals # Time intervals in months
)
