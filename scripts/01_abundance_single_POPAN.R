## Stime di popolazione con RMark, modello POPAN
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
    CH <- as.matrix(data[, -1]) # togli la colonna Sex
  }

  {
    # dates and time intervals
    occasion_dates <- as.Date(c(
      # "2019-07-15", "2019-08-05", "2019-08-20", "2019-10-10",
      # "2020-09-15",
      # "2021-09-15",
      # "2022-05-05",
      "2022-05-20",
      "2022-06-05",
      "2022-06-20",
      "2022-07-05",
      "2022-07-20",
      "2022-08-05",
      "2022-08-20",
      "2022-09-05",
      "2022-09-20",
      "2022-10-15", #"2022-10-20",
      # "2022-11-05", "2022-11-20", "2022-12-05",
      "2023-05-15",
      "2023-06-15",
      "2023-07-15",
      "2023-08-05",
      "2023-08-20",
      "2023-09-05",
      "2023-09-20",
      # "2023-10-10",
      "2024-06-15",
      "2024-07-15",
      "2024-08-15",
      "2024-09-15" #, "2024-10-15"
    ))

    time_intervals <- intervals_months(occasion_dates)
    # names(time_intervals) <- paste(head(occasion_dates, -1), tail(occasion_dates, -1), sep = " to ")
  }

  # dataframe for RMark
  {
    ch_strings <- apply(CH, 1, paste0, collapse = "")
    # data_ch <- data.frame(ch = ch_strings, sex = as.matrix(data[,1]), freq = 1)
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
      "2022-10-15", #"2022-10-20",
      # "2022-11-05", "2022-11-20", "2022-12-05",
      "2023-05-15",
      "2023-06-15",
      "2023-07-15",
      "2023-08-05",
      "2023-08-20",
      "2023-09-05",
      "2023-09-20",
      # "2023-10-10",
      "2024-06-15",
      "2024-07-15",
      "2024-08-15",
      "2024-09-15" #, "2024-10-15"
    ))

    time_intervals <- intervals_months(occasion_dates)
    # names(time_intervals) <- paste(head(occasion_dates, -1), tail(occasion_dates, -1), sep = " to ")
  }

  # dataframe for RMark
  {
    ch_strings <- apply(CH, 1, paste0, collapse = "")
    # data_ch <- data.frame(ch = ch_strings, sex = as.matrix(data[,1]), freq = 1)
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
  # processed <- process.data(data_ch, model="POPAN", groups = ("Sex"), time.intervals = time_intervals)
  processed <- process.data(
    data_ch,
    model = "POPAN",
    time.intervals = time_intervals
  )
  ddl <- make.design.data(processed)
}

# Fit models
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

## ---- plots for the paper: Figure S6 ----

library(RColorBrewer)
# Visualizza i codici hex per 8 colori della palette Dark2
brewer.pal(n = 6, name = "Dark2")
display.brewer.pal(n = 6, name = "Dark2")

R3_POPAN <- read_xlsx("outputs/POPAN_abund_dataR3.xlsx") |>
  mutate(area_name = "R1948")
S3_POPAN <- read_xlsx("outputs/POPAN_abund_dataS3.xlsx") |>
  mutate(area_name = "S1966")

plot_data <- bind_rows(R3_POPAN, S3_POPAN) |>
  mutate(
    Facet_Label = paste0(
      "<b>",
      area_name,
      "</b>" # Nome in grassetto e a capo
    )
  )
{
  p <- ggplot(plot_data, aes(x = date, y = estimate, group = Facet_Label)) +
    geom_ribbon(
      aes(ymin = lcl, ymax = ucl, fill = Facet_Label),
      alpha = 0.2
    ) +
    geom_line(aes(color = Facet_Label), linewidth = 1) +
    geom_point(aes(color = Facet_Label), size = 1.5) +
    scale_color_manual(values = c("#7570B3", "#E6AB02")) +
    # scale_color_brewer(palette = "Dark2") +
    # scale_fill_brewer(palette = "Dark2") +
    scale_fill_manual(values = c("#7570B3", "#E6AB02")) +
    facet_wrap(~Facet_Label, nrow = 2, scales = "free_y") + # scales = "fixed"
    # ggh4x::facetted_pos_scales(
    #   y = list(
    #     Facet_Label %in%
    #       grep(
    #         "S1966",
    #         unique(plot_data$Facet_Label),
    #         value = T
    #       ) ~ scale_y_continuous(limits = c(0, 35)),
    #     Facet_Label %in%
    #       grep(
    #         "R1948",
    #         unique(plot_data$Facet_Label),
    #         value = T
    #       ) ~ scale_y_continuous(limits = c(0, 85))
    #   )
    # ) +
    labs(
      x = "",
      y = "Estimated Population Size",
      # title = "Hierarchical Multi-Site Abundance Estimates",
      # subtitle = "Borrowing strength across study areas via global mortality rate"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "none", # Facet labels already identify the area
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      strip.text = element_markdown(lineheight = 1.2),
      strip.background = element_rect(fill = "grey90", color = NA)
    )

  print(p)
}
# ggsave(
#   filename = "figs/Figure_S6.tiff", # o .png, .pdf, .eps
#   plot = last_plot(), # esporta l'ultimo grafico visualizzato
#   device = "tiff", # formato file
#   compression = "lzw",
#   width = 180, # larghezza (es. 180mm è lo standard "full page width")
#   height = 120, # altezza in mm
#   units = "mm", # unità di misura
#   dpi = 300, # risoluzione richiesta
#   bg = "white" # sfondo bianco (evita trasparenze indesiderate)
# )
