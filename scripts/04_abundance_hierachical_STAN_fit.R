## Abundance estimation - fit hierarchical STAN model

rm(list = ls())

# Run the data preparation script
source("scripts/03_abundance_hierachical_STAN_data_prep.R")

# load  additional packages
{
  # library()
}

## ---- Analysis with STAN ----

# Fit model
{
  mod <- cmdstan_model("models/abundance_hierarchical.stan")

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
saveRDS(fit, file = "outputs/fit_hierarchical.rds")

# Abundance estimates
{
  N_draws <- fit$draws("N") # Extract draws for N (iterations x chains x sessions)
  N_mat <- posterior::as_draws_matrix(N_draws) # Combine chains & iterations

  # Convert to tidy format: one row per draw per session
  # In a hierarchical model, Stan names matrix columns as N[area,session]
  df_long <- as.data.frame(N_mat) %>%
    pivot_longer(
      cols = starts_with("N["),
      names_to = "parameter",
      values_to = "N"
    ) %>%
    # Extract area index and session index from the string "N[area,session]"
    mutate(
      indices = gsub("N\\[|\\]", "", parameter),
      area_num = as.integer(sub(",.*", "", indices)),
      session_num = as.integer(sub(".*,", "", indices))
    ) %>%
    # Optional: map area numbers back to names
    mutate(area_name = areas[area_num])

  # Summarize statistics per area and per session
  {
    summary_stats <- df_long %>%
      group_by(area_name, session_num) %>%
      summarise(
        mean = mean(N),
        lower_95 = quantile(N, 0.025),
        upper_95 = quantile(N, 0.975),
        .groups = "drop"
      ) %>%
      # Add the dates back in (using the master_dates list)
      mutate(date = master_dates[session_num])

    print(summary_stats)
  }

  # Multi-Site Plot (Using facet_wrap to see all areas at once)
  {
    p <- ggplot(summary_stats, aes(x = date, y = mean, group = area_name)) +
      geom_ribbon(
        aes(ymin = lower_95, ymax = upper_95, fill = area_name),
        alpha = 0.2
      ) +
      geom_line(aes(color = area_name), linewidth = 1) +
      facet_wrap(~area_name, scales = "free_y") + # 'free_y' is helpful if areas have very different sizes
      labs(
        x = "Date",
        y = "Estimated abundance (N)",
        title = "Hierarchical Multi-Site Abundance Estimates",
        subtitle = "Borrowing strength across study areas via global mortality rate"
      ) +
      theme_minimal() +
      theme(
        legend.position = "none", # Facet labels already identify the area
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
      )

    print(p)
  }
}

# Optional
# writexl::write_xlsx(summary_stats, "outputs/STAN_hierarchical.xlsx")

# Filtriamo solo i parametri N per un Caterpillar Plot delle abbondanze medie
n_summary <- summary_diag %>%
  filter(grepl("N\\[", variable)) %>%
  # Estraiamo gli indici per ordinare o filtrare se necessario
  mutate(
    indices = gsub("N\\[|\\]", "", variable),
    area_num = as.integer(sub(",.*", "", indices)),
    session_num = as.integer(sub(".*,", "", indices)),
    area_name = areas[area_num],
    date = master_dates[session_num]
  )

# Plot: Caterpillar plot per area (ultime 5 sessioni come esempio)
ggplot(
  n_summary, # %>% filter(session_num > (max(session_num) - 23)),
  aes(x = date, y = mean, color = area_name)
) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = q5, ymax = q95), width = 0.2) +
  facet_wrap(~area_name, scales = "free_y") +
  labs(
    title = "Caterpillar Plot: Stime di Abbondanza (N)",
    subtitle = "Media posteriore e intervalli di credibilità al 90% (q5-q95)",
    x = "Sessione",
    y = "Abbondanza stimata"
  ) +
  theme_light() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )


##---- grafico per la sopravvivenza media ----
# Estraiamo i parametri phi (che variano in base a delta e mu)
phi_summary <- summary_diag %>%
  filter(grepl("phi", variable)) %>%
  mutate(
    interval_num = as.integer(gsub("phi\\[|\\]", "", variable)),
    # Associamo l'intervallo alla data (inizio dell'intervallo)
    date = master_dates[interval_num]
  )

# Grafico della sopravvivenza nel tempo
ggplot(phi_summary, aes(x = date, y = mean)) +
  geom_line(color = "darkgreen", size = 1) +
  geom_point(color = "darkgreen") +
  geom_ribbon(aes(ymin = q5, ymax = q95), fill = "darkgreen", alpha = 0.2) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Sopravvivenza Apparente (phi) nel tempo",
    subtitle = "Basata sul tasso di mortalità globale (mu) e gli intervalli temporali (delta)",
    x = "Data di inizio intervallo",
    y = "Probabilità di sopravvivenza"
  ) +
  theme_minimal()
# sopravvivenza nel tempo per aree?

# creare ancora grafici detectability e mu
