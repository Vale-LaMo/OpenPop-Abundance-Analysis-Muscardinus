## Population estimate with RMark, POPAN model
{
  # load packages
  library(readxl)
  library(lubridate)
  library(RMark)
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

## ---- Data preparation ----
lista_data_ch <- vector("list", length(areas))
{
  for (f in areas) {
    temp_data <- read_xlsx("data/data_abundance_art1.xlsx", sheet = f)
    temp_data[is.na(temp_data)] <- 0
    # Create capture histories excluding the first column (sex)
    ch_vec <- apply(temp_data[, -1], 1, paste0, collapse = "")
    lista_data_ch[[f]] <- data.frame(
      ch = as.character(ch_vec),
      area = f,
      freq = 1,
      stringsAsFactors = FALSE
    )
  }
  all_data_ch <- bind_rows(lista_data_ch)
  all_data_ch$area <- as.factor(all_data_ch$area)
  all_data_ch <- all_data_ch %>%
    mutate(
      valley = substr(area, 1, 1), # Estrae 'S' o 'R'
      altitude = as.factor(substr(area, 2, 2)) # Estrae '1', '2' o '3'
    )

  # check data
  str(all_data_ch)
  summary(all_data_ch$area) # how many individuals per area
}


## ---- POPAN analysis ----
{
  # Data preparation for POPAN
  dormouse.proc <- process.data(
    all_data_ch,
    model = "POPAN",
    groups = c("valley", "altitude"), # area
    time.intervals = time_intervals
  )
  dormouse.ddl <- make.design.data(dormouse.proc)
  dormouse.ddl$Phi$Season <- "active" # default value
  # identify hibernation intervals
  dormouse.ddl$Phi$Season[
    dormouse.ddl$Phi$par.index %in%
      c(12, 19, 34, 41, 56, 63, 78, 85, 100, 107, 122, 129)
  ] <- "hibernation"
  dormouse.ddl$Phi$Season <- as.factor(dormouse.ddl$Phi$Season)
  head(dormouse.ddl$Phi)
}

# Fit models
{
  model_null <- mark(
    processed,
    ddl, # Base model, simpler (no time varying parameters)
    model.parameters = list(
      Phi = list(formula = ~1),
      p = list(formula = ~1),
      pent = list(formula = ~1),
      N = list(formula = ~ altidude * valley)
    )
  )
  model_constant_p <- mark(
    processed,
    ddl,
    model.parameters = list(
      Phi = list(formula = ~Season),
      p = list(formula = ~1),
      pent = list(formula = ~time),
      N = list(formula = ~ altidude * valley)
    )
  )
  model_basic_demo <- mark(
    processed,
    ddl,
    model.parameters = list(
      Phi = list(formula = ~Season),
      p = list(formula = ~time),
      pent = list(formula = ~time),
      N = list(formula = ~ altidude * valley)
    )
  )
  model_common_det_additive_surv <- mark(
    processed,
    ddl,
    model.parameters = list(
      Phi = list(formula = ~ Season + altitude),
      p = list(formula = ~time),
      pent = list(formula = ~time),
      N = list(formula = ~ altidude * valley)
    )
  )
  model_common_det_altitudinal_surv <- mark(
    processed,
    ddl,
    model.parameters = list(
      Phi = list(formula = ~ Season * altitude),
      p = list(formula = ~time),
      pent = list(formula = ~time),
      N = list(formula = ~ altidude * valley)
    )
  )
  model_seasonal_surv <- mark(
    processed,
    ddl,
    model.parameters = list(
      Phi = list(formula = ~Season),
      p = list(formula = ~ time + altitude),
      pent = list(formula = ~time),
      N = list(formula = ~ altidude * valley)
    )
  )
  model_additive_altitudinal_surv <- mark(
    processed,
    ddl,
    model.parameters = list(
      Phi = list(formula = ~ Season + altitude),
      p = list(formula = ~ time + altitude),
      pent = list(formula = ~time),
      N = list(formula = ~ altidude * valley)
    )
  )
  model_full <- mark(
    processed,
    ddl,
    model.parameters = list(
      Phi = list(formula = ~ Season * altitude),
      p = list(formula = ~ time + altitude),
      pent = list(formula = ~time),
      N = list(formula = ~ altidude * valley)
    )
  )
}

# Model comparison
{
  mod_comparison <- collect.models()
  print(mod_comparison)
}

# PICK THE BEST MODEL: CHANGE MANUALLY!!
{
  # best_mod <- model_time_phi
  best_mod <- model
}

# Results
{
  # Model summary
  print(summary(best_mod))

  # # Interpretation of the parameters for the constant parameter model, see below for plots for more complex models
  # print(plogis(1.747))   # Phi: apparent average survival
  # plogis(-0.732)  # p: average capture probability
  # plogis(-0.312)  # pent: average probability of entry into the population
  # exp(1.488)      # N: average estimated N (log-scale → exp)

  {
    # abundance estimates across time
    abund_data <- best_mod$results$derived$N
    if (selected_area == "S3") {
      abund_data$date <- occasion_dates_S3
    }
    if (selected_area == "R3") {
      abund_data$date <- occasion_dates_R3
    }

    print(
      ggplot(abund_data, aes(x = date, y = estimate)) +
        geom_line(color = "steelblue", size = 1) +
        geom_point(size = 2, color = "steelblue") +
        geom_ribbon(
          aes(ymin = lcl, ymax = ucl),
          fill = "steelblue",
          alpha = 0.2
        ) +
        labs(
          title = "Estimated abundance (POPAN)",
          x = "Date",
          y = "Estimated N"
        ) +
        theme_minimal()
    )

    print(abund_data)
  }
}

# Optional
# writexl::write_xlsx(abund_data,
#                     paste("results/POPAN_abund_data", selected_area, ".xlsx",
#                           sep = ""))
