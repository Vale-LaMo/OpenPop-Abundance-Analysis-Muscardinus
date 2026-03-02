## Abundance estimation - fit hierarchical STAN model

rm(list = ls())

# Run the data preparation script
source("scripts/03_abundance_hierachical_STAN_data_prep.R")

# # load  additional packages
# {
#   # library()
# }

# load model
{
  fit <- readRDS("outputs/rds/fit_mod_inform.rds")
  # fit <- readRDS("outputs/rds/fit_mod_wide.rds")
}

##---- Abundance estimates ----
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

    print("Estimated N per area and session - summary")
    print(summary_stats)
  }

  summary_diag <- fit$summary(
    variables = NULL, # NULL estrae tutti i parametri inclusi N, phi, p, beta, mu
    "mean",
    "median",
    "sd",
    "quantile2",
    "rhat",
    "ess_bulk"
  )
  print("Parameters summary")
  print(summary_diag)

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
      title = "Caterpillar plot: Hierarchical Multi-Site Abundance Estimates",
      subtitle = "Posterior means and 90% credibility intervals",
      x = "Session",
      y = "Estimated abundance (N)"
    ) +
    theme_light() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) -> caterpillar_plot
  print(caterpillar_plot)

  ##---- Super-population ----
  {
    n_super_summary <- summary_diag %>%
      filter(grepl("N_super", variable)) %>%
      mutate(
        area_num = as.integer(gsub("N_super\\[|\\]", "", variable)),
        area_name = areas[area_num]
      )
    print("Super Population per area")
    print(n_super_summary)

    ggplot(n_super_summary, aes(x = area_name, y = mean, color = area_name)) +
      geom_point(size = 4) +
      geom_errorbar(aes(ymin = q5, ymax = q95), width = 0.2) +
      labs(
        title = "Super-population Estimate (N_super) per Area",
        subtitle = "Total individuals ever present (observed + estimated)",
        y = "Number of individuals",
        x = "Study Area"
      ) +
      theme_light() -> superP_plot
    print(superP_plot)
  }

  # Multi-Site Plot (Using facet_wrap to see all areas at once)
  n_super_summary_plot <- n_super_summary |>
    dplyr::rename(mean_SuperP = mean, q5_SuperP = q5, q95_SuperP = q95)
  plot_data <- summary_stats |>
    left_join(n_super_summary_plot, by = "area_name") |>
    mutate(
      Area_name = case_when(
        area_name == "S1" ~ "S1521",
        area_name == "S2" ~ "S1868",
        area_name == "S3" ~ "S1966",
        area_name == "R1" ~ "R1200",
        area_name == "R2" ~ "R1840",
        area_name == "R3" ~ "R1948"
      ),
      Facet_Label = paste0(
        "<b>",
        Area_name,
        "</b><br>", # Nome in grassetto e a capo
        "<span style='font-weight:normal; font-size:8pt;'>", # Inizio stile normale e più piccolo
        "Super-N: ",
        round(mean_SuperP, 0),
        " (CrI: ",
        round(q5_SuperP, 0),
        "-",
        round(q95_SuperP, 0),
        ")",
        "</span>"
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
      scale_color_brewer(palette = "Dark2") +
      scale_fill_brewer(palette = "Dark2") +
      facet_wrap(~Facet_Label, ncol = 3, scales = "free_y") + # scales = "fixed"
      ggh4x::facetted_pos_scales(
        y = list(
          Facet_Label %in%
            grep(
              "S1521",
              unique(plot_data$Facet_Label),
              value = T
            ) ~ scale_y_continuous(limits = c(0, 25)),
          Facet_Label %in%
            grep(
              "S1868",
              unique(plot_data$Facet_Label),
              value = T
            ) ~ scale_y_continuous(limits = c(0, 25)),
          Facet_Label %in%
            grep(
              "S1966",
              unique(plot_data$Facet_Label),
              value = T
            ) ~ scale_y_continuous(limits = c(0, 25)),
          Facet_Label %in%
            grep(
              "R1200",
              unique(plot_data$Facet_Label),
              value = T
            ) ~ scale_y_continuous(limits = c(0, 25)),
          Facet_Label %in%
            grep(
              "R1840",
              unique(plot_data$Facet_Label),
              value = T
            ) ~ scale_y_continuous(limits = c(0, 25)),
          Facet_Label %in%
            grep(
              "R1948",
              unique(plot_data$Facet_Label),
              value = T
            ) ~ scale_y_continuous(limits = c(0, 65))
        )
      ) +
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
  #   filename = "figs/Abundance_Plot_Bayes.tiff", # o .png, .pdf, .eps
  #   plot = last_plot(),                          # esporta l'ultimo grafico visualizzato
  #   device = "tiff",                             # formato file
  #   compression = "lzw",
  #   width = 180,                                 # larghezza (es. 180mm è lo standard "full page width")
  #   height = 120,                                # altezza in mm
  #   units = "mm",                                # unità di misura
  #   dpi = 300,                                   # risoluzione richiesta
  #   bg = "white"                                 # sfondo bianco (evita trasparenze indesiderate)
  # )
}

# Optional
# writexl::write_xlsx(summary_stats, "outputs/STAN_hierarchical_mod_inform.xlsx")

##---- detection probability ----
{
  summary_p <- fit$summary(
    variables = "alpha_p", # NULL estrae tutti i parametri inclusi N, phi, p, beta, mu
    "mean",
    "median",
    "sd",
    "quantile2",
    "rhat",
    "ess_bulk",
    "ess_tail"
  )
  print("Parameters summary - detectability")
  print(summary_p)

  p_summary <- summary_diag %>%
    filter(grepl("p_eff", variable)) %>%
    mutate(
      indices = gsub("p_eff\\[|\\]", "", variable),
      area_num = as.integer(sub(",.*", "", indices)),
      session_num = as.integer(sub(".*,", "", indices)),
      area_name = areas[area_num],
      date = master_dates[session_num]
    )
  print("Parameters summary (probs) - detectability")
  print(p_summary)

  ggplot(p_summary, aes(x = date, y = mean, color = area_name)) +
    geom_line() +
    geom_ribbon(
      aes(ymin = q5, ymax = q95, fill = area_name),
      alpha = 0.1,
      color = NA
    ) +
    facet_wrap(~area_name) +
    labs(
      title = "Detection Probability (p) per Area",
      y = "Probability",
      x = "Date"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) -> detection_plot
  print(detection_plot)

  # tabelle detection probability
  {
    # 1. Estrazione dei draws per p_eff
    p_draws <- as_draws_matrix(fit$draws("p_eff"))

    # 2. Definizione dei raggruppamenti (Intervalli e Aree)
    # Usa gli stessi indici definiti per phi
    idx_summer <- c(1:11, 13:18, 19:23)
    idx_winter <- c(12, 19)
    num_areas <- length(unique(p_summary$area_num))
    num_sessions <- length(unique(p_summary$session_num))

    # Funzione per mappare gli indici della matrice p_eff[area, sessione]
    get_p_draws <- function(a, s_indices) {
      # Costruisce i nomi dei parametri p_eff[a,s] per l'area 'a' e le sessioni 's_indices'
      cols <- paste0("p_eff[", a, ",", s_indices, "]")
      # Fa la media riga per riga (per ogni draw)
      return(rowMeans(p_draws[, cols, drop = FALSE]))
    }

    # 3. Calcolo delle statistiche aggregate
    res_list <- list()

    for (a in 1:num_areas) {
      # Global per Area (tutte le sessioni)
      draws_area <- get_p_draws(a, 1:num_sessions)

      # Summer per Area
      draws_summer <- get_p_draws(a, idx_summer)

      # Winter per Area
      draws_winter <- get_p_draws(a, idx_winter)

      # Salvataggio risultati
      res_list[[length(res_list) + 1]] <- data.frame(
        Area = areas[a],
        Season = "Global",
        Mean = mean(draws_area),
        Q5 = quantile(draws_area, 0.05),
        Q95 = quantile(draws_area, 0.95)
      )
      res_list[[length(res_list) + 1]] <- data.frame(
        Area = areas[a],
        Season = "Summer",
        Mean = mean(draws_summer),
        Q5 = quantile(draws_summer, 0.05),
        Q95 = quantile(draws_summer, 0.95)
      )
      res_list[[length(res_list) + 1]] <- data.frame(
        Area = areas[a],
        Season = "Winter",
        Mean = mean(draws_winter),
        Q5 = quantile(draws_winter, 0.05),
        Q95 = quantile(draws_winter, 0.95)
      )
    }

    p_aggregated_results <- bind_rows(res_list)
    print("Detection probability - global and seasonal")
    print(p_aggregated_results)
    # print(filter(p_aggregated_results, Season == "Global"))
    # print(filter(p_aggregated_results, Season == "Summer"))
    # print(filter(p_aggregated_results, Season == "Winter"))
  }
}


##---- entry probability ----
{
  beta_summary <- summary_diag %>%
    filter(grepl("beta", variable)) %>%
    mutate(
      indices = gsub("beta\\[|\\]", "", variable),
      area_num = as.integer(sub(",.*", "", indices)),
      session_num = as.integer(sub(".*,", "", indices)),
      area_name = areas[area_num],
      date = master_dates[session_num]
    )
  print("Parameters summary (probs) - entry")
  print(beta_summary)

  ggplot(beta_summary, aes(x = date, y = mean, fill = area_name)) +
    geom_col() + # Le entrate sono impulsi, il grafico a barre (col) è spesso più chiaro
    facet_wrap(~area_name) +
    labs(
      title = "Recruitment Pulses (beta)",
      y = "Entry Probability",
      x = "Date"
    ) +
    theme_minimal() +
    theme(legend.position = "none") -> recruitment_pulses
  print(recruitment_pulses)

  # calcolo del reclutamento reale per sessione, per confronto con Mark POPAN
  {
    # 1. Estrazione dei draws per beta e N_super
    # Assicurati che i nomi corrispondano a quelli nel tuo oggetto 'fit'
    beta_draws <- as_draws_matrix(fit$draws("beta"))
    n_super_draws <- as_draws_matrix(fit$draws("N_super"))

    n_iterations <- nrow(beta_draws)
    n_areas <- 6
    n_sessions <- 23

    # 2. Funzione per calcolare B per tutti i draw di un'area
    get_B_uncertainty <- function(area_idx) {
      # Matrice per salvare i risultati (iterazioni x sessioni)
      B_matrix <- matrix(NA, nrow = n_iterations, ncol = n_sessions)

      # Estraiamo i draw di N_super per questa area
      # (Verifica il nome esatto della colonna, potrebbe essere "N_super[area_idx]")
      ns_col <- paste0("N_super[", area_idx, "]")
      area_ns_draws <- n_super_draws[, ns_col]

      for (draw in 1:n_iterations) {
        # Pool iniziale per questo draw
        pool_curr <- area_ns_draws[draw]

        for (k in 1:n_sessions) {
          # Prendiamo il beta specifico per area, sessione e draw
          beta_col <- paste0("beta[", area_idx, ",", k, "]")
          b_val <- beta_draws[draw, beta_col]

          # Calcolo del reclutamento netto (B)
          B_matrix[draw, k] <- pool_curr * b_val

          # Aggiornamento del pool rimanente per la sessione successiva
          pool_curr <- pool_curr * (1 - b_val)
        }
      }

      # Calcolo dei riassunti statistici per ogni sessione
      B_summary <- as.data.frame(t(apply(B_matrix, 2, function(x) {
        c(mean = mean(x), q5 = quantile(x, 0.05), q95 = quantile(x, 0.95))
      })))

      B_summary$session_num <- 1:n_sessions
      B_summary$area_num <- area_idx
      return(B_summary)
    }

    # 3. Esecuzione per tutte le aree
    B_final_with_error <- bind_rows(lapply(1:n_areas, get_B_uncertainty)) %>%
      mutate(area_name = areas[area_num])

    # 4. Confronto finale (Media del reclutamento post-iniziale)
    B_comparison_paper <- B_final_with_error %>%
      filter(session_num > 1) %>% # Escludiamo lo stock iniziale per confrontare con MARK
      group_by(area_name) %>%
      summarise(
        B_mean = mean(mean),
        B_q5 = mean(`q5.5%`),
        B_q95 = mean(`q95.95%`)
      )

    print("Net recruitment (B) for comparison with Mark")
    print(B_comparison_paper)
  }

  # estrazione delle stime di beta
  {
    # 1. Estrazione dei draws di beta
    beta_draws <- as_draws_matrix(fit$draws("beta"))

    # 2. Parametri dello studio
    num_areas <- 6
    # Gli intervalli totali sono solitamente N_sessioni - 1 (nel tuo caso 22)
    n_intervalli_tot <- 22
    n_intervalli_summer <- length(idx_summer) # quanti intervalli estivi hai (es. 20)

    # 3. Funzione per il calcolo del beta medio
    calculate_beta_avg <- function(area_idx, interval_indices) {
      # Identifica le colonne beta[area, sessione] corrette
      cols <- paste0("beta[", area_idx, ",", interval_indices, "]")
      # Filtra solo le colonne esistenti nel modello
      cols <- cols[cols %in% colnames(beta_draws)]

      # Calcola la media per ogni draw (riga)
      draws_avg <- rowMeans(beta_draws[, cols])

      # Restituisce il summary
      return(data.frame(
        Mean = mean(draws_avg),
        Q5 = quantile(draws_avg, 0.05),
        Q95 = quantile(draws_avg, 0.95)
      ))
    }

    # 4. Esempio: Calcolo per Area R3 (R1948)
    beta_avg_R3_summer <- calculate_beta_avg(6, idx_summer)
    beta_avg_R3_global <- calculate_beta_avg(6, 1:n_intervalli_tot)
  }

  print("Average betas (entry probs from Super pop) per area")
  print(calculate_beta_avg(1, 1:n_intervalli_tot))
  print(calculate_beta_avg(2, 1:n_intervalli_tot))
  print(calculate_beta_avg(3, 1:n_intervalli_tot))
  print(calculate_beta_avg(4, 1:n_intervalli_tot))
  print(calculate_beta_avg(5, 1:n_intervalli_tot))
  print(calculate_beta_avg(6, 1:n_intervalli_tot))
}

##---- survival ----
{
  mu_summary <- summary_diag %>%
    filter(variable %in% c("mu_active", "mu_hibernation"))

  ggplot(mu_summary, aes(x = variable, y = mean, fill = variable)) +
    geom_bar(stat = "identity", width = 0.5) +
    geom_errorbar(aes(ymin = q5, ymax = q95), width = 0.1) +
    labs(
      title = "Instantaneous Mortality Rates (mu)",
      subtitle = "Comparison between active season and hibernation",
      y = "Mortality Rate (mu)",
      x = ""
    ) +
    scale_fill_manual(
      values = c("mu_active" = "orange", "mu_hibernation" = "darkblue")
    ) +
    theme_minimal() -> istantaneous_seasonal_mu
  print(istantaneous_seasonal_mu)

  # Estraiamo i parametri phi (che variano in base a delta e mu)
  phi_summary <- summary_diag %>%
    filter(grepl("phi", variable)) %>%
    mutate(
      interval_num = as.integer(gsub("phi\\[|\\]", "", variable)),
      # Associamo l'intervallo alla data (inizio dell'intervallo)
      date = master_dates[interval_num]
    )
  print("Parameters summary (probs) - survival")
  print(phi_summary)

  # Grafico della sopravvivenza nel tempo
  ggplot(phi_summary, aes(x = date, y = mean)) +
    geom_line(color = "darkgreen", size = 1) +
    geom_point(color = "darkgreen") +
    geom_ribbon(aes(ymin = q5, ymax = q95), fill = "darkgreen", alpha = 0.2) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = "Apparent survival",
      subtitle = "Based on global mortality rate (mu) and temporal intervals",
      x = "",
      y = "Survival probability"
    ) +
    theme_minimal() -> survival_plot
  print(survival_plot)

  # table for summer and winter survival
  {
    # 1. Definisci gli indici (come fatto prima)
    idx_summer <- c(1:10, 12:17, 19:22)
    idx_winter <- c(11, 18)

    # 2. Estrai i draws completi per phi
    phi_draws <- fit$draws("phi") # Formato draws_array

    # 3. Calcola la media stagionale PER OGNI CAMPIONE della posterior
    # Questo preserva la corretta incertezza statistica
    summer_draws_avg <- rowMeans(as.matrix(phi_draws[,, idx_summer]))
    winter_draws_avg <- rowMeans(as.matrix(phi_draws[,, idx_winter]))

    # 4. Crea la tabella riassuntiva finale
    seasonal_summary <- data.frame(
      Season = c("Summer", "Winter"),
      Mean = c(mean(summer_draws_avg), mean(winter_draws_avg)),
      Q5 = c(
        quantile(summer_draws_avg, 0.05),
        quantile(winter_draws_avg, 0.05)
      ),
      Q95 = c(
        quantile(summer_draws_avg, 0.95),
        quantile(winter_draws_avg, 0.95)
      )
    )

    print("Seasonal survival summary")
    print(seasonal_summary)
  }

  # table for global survival
  {
    # 2. Estrai i draws completi per phi
    phi_draws <- fit$draws("phi") # Formato draws_array

    # 3. Calcola la media stagionale PER OGNI CAMPIONE della posterior
    # Questo preserva la corretta incertezza statistica
    draws_avg <- rowMeans(as.matrix(phi_draws))

    # 4. Crea la tabella riassuntiva finale
    global_phi_summary <- data.frame(
      Mean = c(mean(draws_avg)),
      Q5 = c(quantile(draws_avg, 0.05)),
      Q95 = c(quantile(draws_avg, 0.95))
    )

    print("Global survival summary")
    print(global_phi_summary)
  }
}
