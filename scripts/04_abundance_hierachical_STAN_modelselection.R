# STAN model selection

library(cmdstanr)
library(loo)
library(dplyr)
library(gt)

# Crea le cartelle se non esistono
if (!dir.exists("outputs")) {
  dir.exists("outputs")
}
if (!dir.exists("outputs/rds")) {
  dir.create("outputs/rds")
}

# Definisci le configurazioni delle prior da testare
scenari <- list(
  mod_base = list(
    p_prior_mean = 0,
    p_prior_sd = 1.5,
    mu_active_mean = 0,
    mu_active_sd = 1.5,
    mu_hibernation_mean = 0,
    mu_hibernation_sd = 1.5
  ),
  mod_inform = list(
    p_prior_mean = -1,
    p_prior_sd = 1,
    mu_active_mean = 0,
    mu_active_sd = 1.5,
    mu_hibernation_mean = 0,
    mu_hibernation_sd = 1.5
  ),
  mod_regular = list(
    p_prior_mean = 0,
    p_prior_sd = 1,
    mu_active_mean = 0,
    mu_active_sd = 1,
    mu_hibernation_mean = 0,
    mu_hibernation_sd = 1
  ),
  mod_strong = list(
    p_prior_mean = 0,
    p_prior_sd = 0.5,
    mu_active_mean = 0,
    mu_active_sd = 0.5,
    mu_hibernation_mean = 0,
    mu_hibernation_sd = 0.5
  ),
  mod_wide = list(
    p_prior_mean = 0,
    p_prior_sd = 2,
    mu_active_mean = 0,
    mu_active_sd = 2,
    mu_hibernation_mean = 0,
    mu_hibernation_sd = 2
  )
)

source("scripts/03_abundance_hierachical_STAN_data_prep.R")

# Carica il modello (compila una volta sola!)
stan_file <- "models/abundance_hierarchical.stan"
mod <- cmdstan_model(stan_file)

# Ciclo di esecuzione
risultati_loo <- list()

for (nome in names(scenari)) {
  message("\n>>> RUNNING SCENARIO: ", nome)

  # Uniamo i dati originali con i parametri della prior corrente
  stan_data_current <- c(stan_data, scenari[[nome]])

  # Lancio del campionamento
  fit <- mod$sample(
    data = stan_data_current,
    seed = 123,
    chains = 4,
    parallel_chains = 2,
    iter_warmup = 2000,
    iter_sampling = 3000
  )

  # SALVATAGGIO PERMANENTE (fondamentale per evitare l'errore dei file temporanei)
  file_rds <- paste0("outputs/rds/fit_", nome, ".rds")
  fit$save_object(file = file_rds)

  # Calcolo e salvataggio del LOO
  risultati_loo[[nome]] <- fit$loo()
}

# CONFRONTO FINALE
confronto <- loo_compare(risultati_loo)
print(confronto, simplify = FALSE)
# Convertiamo il confronto in un data frame
confronto_df <- as.data.frame(confronto)
# Aggiungiamo i nomi dei modelli come prima colonna
confronto_df <- cbind(model = rownames(confronto_df), confronto_df)
# writexl::write_xlsx(confronto_df, "outputs/loo_comparison.xlsx")
# saveRDS(confronto, "outputs/loo_comparison.rds")
# saveRDS(risultati_loo, "outputs/loo_results.rds")

print(risultati_loo$mod_base)
print(risultati_loo$mod_inform)
print(risultati_loo$mod_regular)
print(risultati_loo$mod_strong)
print(risultati_loo$mod_wide)

# Preparazione dati per tabella gt
confronto <- readRDS("outputs/loo_comparison.rds")
confronto_df <- as.data.frame(confronto)
confronto_df <- cbind(model = rownames(confronto_df), confronto_df)
confronto_df$model <- rownames(confronto_df)
confronto_df <- confronto_df[, c(
  "model",
  "elpd_diff",
  "se_diff",
  "elpd_loo",
  "se_elpd_loo",
  "p_loo",
  "looic"
)]

confronto_df %>%
  gt() %>%
  # tab_header(
  #   title = md("**Model Comparison via LOO-CV**"),
  #   subtitle = "Sensitivity analysis of prior distributions"
  # ) %>%
  # Riduzione della dimensione del font e spaziatura
  tab_options(
    table.font.size = px(10),
    data_row.padding = px(2.5),
    column_labels.font.size = px(10),
    column_labels.font.weight = "bold"
  ) %>%
  fmt_number(
    columns = where(is.numeric),
    decimals = 1
  ) %>%
  cols_label(
    model = "Model Scenario",
    elpd_diff = "Δ ELPD",
    se_diff = "SE (diff)",
    elpd_loo = "ELPD",
    p_loo = "p_LOO",
    looic = "LOOIC"
  ) %>%
  # Evidenzia il modello vincitore
  tab_style(
    style = cell_fill(color = "lightgrey", alpha = 0.5),
    locations = cells_body(rows = 1)
  ) %>%
  gtsave("outputs/loo_comparison.png")
