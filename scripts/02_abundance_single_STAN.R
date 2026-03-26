## Abundance estimation with STAN (single study areas)
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
}

# Study area R3
{
  # data preparation: study area R3
  {
    # capture histories
    data <- read_xlsx("data/data_abundance.xlsx", sheet = "R3_clean")
    data[is.na(data)] <- 0
    CH <- as.matrix(data[, -1]) # remove Sex column
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
}


## ---- Data preparation for STAN ----
{
  # build augmented matrix

  # Number of capture sessions
  K <- ncol(CH)

  # Number of observed individuals
  n <- nrow(CH)

  # Define M: total number of individuals in augmented dataset
  M <- n * 3 # You can later try 60, 80, etc., to check robustness
  # A common rule of thumb is to set M to about 2–3 times the number of observed individuals
  # to start. You can increase this later if the posterior mass for N is too close to M.

  # Create M-n rows of all zeros
  data_aug <- rbind(CH, matrix(0, nrow = M - n, ncol = K))

  # Check dimensions
  dim(data_aug) # Should be M x K

  # sex <- data$Sex  # assuming it's coded as "M" and "F"
  # sex_num <- ifelse(sex == "M", 1, 0)
  # sex_aug <- c(sex_num, rep(0, M - length(sex_num)))
}

{
  # Stan data
  stan_data <- list(
    M = M,
    K = K,
    y = data_aug,
    delta = time_intervals #,
    # sex = sex_aug
  )
}


## ---- Analysis with STAN ----

# Fit model
{
  mod <- cmdstan_model("models/abundance_single.stan")

  fit <- mod$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 2000,
    iter_sampling = 3000,
    seed = 123,
    refresh = 500,
    save_cmdstan_config = TRUE
  )
}

## ---- Results and diagnostics ----
{
  # model summary
  mod_summary <- fit$summary()
  print(mod_summary)
  # posterior_draws <- fit$draws()
}

fit$cmdstan_diagnose()

# Abundance estimates
{
  # Extract draws for N (iterations x chains x sessions)
  N_draws <- fit$draws("N")
  # Combine chains & iterations
  N_mat <- posterior::as_draws_matrix(N_draws)
  # Convert to tidy format: one row per draw per session
  df <- as.data.frame(N_mat)
  colnames(df) <- paste0("N_session_", 1:ncol(df))

  df_long <- df %>%
    pivot_longer(cols = everything(), names_to = "session", values_to = "N") %>%
    mutate(session_num = as.integer(gsub("N_session_", "", session)))

  # Summarize mean and 95% credible intervals per session
  {
    summary_stats <- df_long %>%
      group_by(session_num) %>%
      summarise(
        mean = mean(N),
        lower_95 = quantile(N, 0.025),
        upper_95 = quantile(N, 0.975)
      )
    summary_stats <- summary_stats %>%
      mutate(date = occasion_dates[session_num])
    print(summary_stats)
  }

  # Plot
  {
    library(ggplot2)
    print(
      ggplot(summary_stats, aes(x = date, y = mean)) +
        geom_line(color = "blue") +
        geom_ribbon(
          aes(ymin = lower_95, ymax = upper_95),
          alpha = 0.3,
          fill = "lightblue"
        ) +
        labs(
          x = "Date",
          y = "Estimated abundance (N)",
          title = "Estimated Abundance Over Time with 95% Credible Interval"
        ) +
        theme_minimal()
    )
  }
}

# Optional
# writexl::write_xlsx(summary_stats, "outputs/STAN_abund_dataS3.xlsx")
# writexl::write_xlsx(summary_stats, "outputs/STAN_abund_dataR3.xlsx")


## ---- plots for the paper: Figure S6 ----

library(RColorBrewer)
# Visualizza i codici hex per 8 colori della palette Dark2
brewer.pal(n = 6, name = "Dark2")
display.brewer.pal(n = 6, name = "Dark2")

R3_STAN <- read_xlsx("outputs/STAN_abund_dataR3.xlsx") |>
  mutate(area_name = "R1948")
S3_STAN <- read_xlsx("outputs/STAN_abund_dataS3.xlsx") |>
  mutate(area_name = "S1966")

plot_data <- bind_rows(R3_STAN, S3_STAN) |>
  mutate(
    Facet_Label = paste0(
      "<b>",
      area_name,
      "</b>" # Nome in grassetto e a capo
    )
  )
{
  p <- ggplot(plot_data, aes(x = date, y = mean, group = Facet_Label)) +
    geom_ribbon(
      aes(ymin = lower_95, ymax = upper_95, fill = Facet_Label),
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
#   filename = "figs/Figure_S7.tiff", # o .png, .pdf, .eps
#   plot = last_plot(), # esporta l'ultimo grafico visualizzato
#   device = "tiff", # formato file
#   compression = "lzw",
#   width = 180, # larghezza (es. 180mm è lo standard "full page width")
#   height = 120, # altezza in mm
#   units = "mm", # unità di misura
#   dpi = 300, # risoluzione richiesta
#   bg = "white" # sfondo bianco (evita trasparenze indesiderate)
# )
