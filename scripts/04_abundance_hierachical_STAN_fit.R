## Abundance estimation - fit hierarchical STAN model

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
    save_cmdstan_config=TRUE
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
  N_draws <- fit$draws("N")                       # Extract draws for N (iterations x chains x sessions)
  N_mat <- posterior::as_draws_matrix(N_draws)    # Combine chains & iterations

  # Convert to tidy format: one row per draw per session
  # In a hierarchical model, Stan names matrix columns as N[area,session]
  df_long <- as.data.frame(N_mat) %>%
    pivot_longer(cols = starts_with("N["), 
                 names_to = "parameter", 
                 values_to = "N") %>%
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
      geom_ribbon(aes(ymin = lower_95, ymax = upper_95, fill = area_name), alpha = 0.2) +
      geom_line(aes(color = area_name), linewidth = 1) +
      facet_wrap(~area_name, scales = "free_y") + # 'free_y' is helpful if areas have very different sizes
      labs(
        x = "Date",
        y = "Estimated abundance (N)",
        title = "Hierarchical Multi-Site Abundance Estimates",
        subtitle = "Borrowing strength across study areas via global mortality rate"
      ) +
      theme_minimal() +
      theme(legend.position = "none") # Facet labels already identify the area
  
    print(p)
  }
  
}  

# Optional
# writexl::write_xlsx(summary_stats, "outputs/STAN_hierarchical.xlsx")
