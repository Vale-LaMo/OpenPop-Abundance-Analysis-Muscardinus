## RMark POPAN
{
  # load packages
  library(readxl)
  library(lubridate)
  library(RMark)
  library(ggtext)
}

{
  # function for time intervals (calculates difference in months)
  intervals_months <- function(dates) {
    diffs <- difftime(dates[-1], dates[-length(dates)], units = "days")
    as.numeric(diffs) / 30.44 # 1 month ≈ 30.44 days (mean Gregorian month)
  }
}

## ---- Data preparation ----
# Study area S3
{
  # data preparation: study area S3
  {
    # capture histories
    data <- read_xlsx("data/data_abundance.xlsx", sheet = "S3_clean")
    # data <- read_xlsx("data/data_abundance.xlsx", sheet = "R3_clean")
    data[is.na(data)] <- 0
    CH <- as.matrix(data[, -1]) # remove Sex
  }

  {
    # dates and time intervals
    occasion_dates <- as.Date(c(
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
      "2023-05-15",
      "2023-06-15",
      "2023-07-15",
      "2023-08-05",
      "2023-08-20",
      "2023-09-05",
      "2023-09-20",
      "2024-06-15",
      "2024-07-15",
      "2024-08-15",
      "2024-09-15"
    ))

    time_intervals <- intervals_months(occasion_dates)
  }

  # dataframe for RMark
  {
    ch_strings <- apply(CH, 1, paste0, collapse = "")
    data_ch <- data.frame(ch = ch_strings, freq = 1)
    data_ch_S3 <- data_ch
    occasion_dates_S3 <- occasion_dates
    time_intervals_S3 <- time_intervals
  }
}

# Study area R3
{
  # data preparation: study area R3
  {
    # capture histories
    data <- read_xlsx("data/data_abundance.xlsx", sheet = "R3_clean")
    data[is.na(data)] <- 0
    CH <- as.matrix(data[, -1]) # togli la colonna Sex
  }

  {
    # dates and time intervals
    occasion_dates <- as.Date(c(
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
      "2023-05-15",
      "2023-06-15",
      "2023-07-15",
      "2023-08-05",
      "2023-08-20",
      "2023-09-05",
      "2023-09-20",
      "2024-06-15",
      "2024-07-15",
      "2024-08-15",
      "2024-09-15"
    ))

    time_intervals <- intervals_months(occasion_dates)
  }

  # dataframe for RMark
  {
    ch_strings <- apply(CH, 1, paste0, collapse = "")
    data_ch <- data.frame(ch = ch_strings, freq = 1)
    data_ch_R3 <- data_ch
    occasion_dates_R3 <- occasion_dates
    time_intervals_R3 <- time_intervals
  }
}


## ---- Select area (repeat for each area afterwards) ----
{
  # CHECK THE NAME OF THE STUDY AREA!! CHANGE MANUALLY
  {
    selected_area <- "R3"
  }

  {
    # pick the right data
    rm(data_ch, occasion_dates, time_intervals)
    if (selected_area == "S3") {
      data_ch <- data_ch_S3
      occasion_dates <- occasion_dates_S3
      time_intervals <- time_intervals_S3
    }
    if (selected_area == "R3") {
      data_ch <- data_ch_R3
      occasion_dates <- occasion_dates_R3
      time_intervals <- time_intervals_R3
    }
  }
}


## ---- POPAN analysis ----
{
  # Data preparation for POPAN
  processed <- process.data(
    data_ch,
    model = "POPAN",
    time.intervals = time_intervals
  )
  ddl <- make.design.data(processed)
}

# Fit models - only if Mark is installed
{
  model <- mark(
    processed,
    ddl, # Base model, simpler (no time varying parameters)
    model.parameters = list(
      Phi = list(formula = ~1),
      p = list(formula = ~1),
      pent = list(formula = ~1),
      N = list(formula = ~1)
    )
  )
  model_time <- mark(
    processed,
    ddl,
    model.parameters = list(
      Phi = list(formula = ~1), # constant survival
      p = list(formula = ~time), # time-dependent capture probability
      pent = list(formula = ~time), # time-dependent entry probability
      N = list(formula = ~1) # constant N - be careful, it is not the estimate of pop size!!
    )
  )
  model_time_phi <- mark(
    processed,
    ddl,
    model.parameters = list(
      Phi = list(formula = ~time), # time-dependent survival
      p = list(formula = ~time), # time-dependent capture probability
      pent = list(formula = ~time), # time-dependent entry probability
      N = list(formula = ~1) # constant N - be careful, it is not the estimate of pop size!!
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
