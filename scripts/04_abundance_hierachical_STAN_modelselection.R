# STAN model selection

library(cmdstanr)
library(loo)
library(dplyr)
library(gt)

if (!dir.exists("outputs")) {
  dir.exists("outputs")
}
if (!dir.exists("outputs/STAN")) {
  dir.create("outputs/STAN")
}
if (!dir.exists("outputs/STAN/rds")) {
  dir.create("outputs/STAN/rds")
}

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

stan_file <- "models/abundance_hierarchical.stan"
mod <- cmdstan_model(stan_file)

# PLEASE NOTE: REFITTING THE MODELS TAKES A BIT OF TIME
risultati_loo <- list()

for (nome in names(scenari)) {
  message("\n>>> RUNNING SCENARIO: ", nome)

  stan_data_current <- c(stan_data, scenari[[nome]])

  fit <- mod$sample(
    data = stan_data_current,
    seed = 123,
    chains = 4,
    parallel_chains = 2,
    iter_warmup = 2000,
    iter_sampling = 3000
  )

  file_rds <- paste0("outputs/STAN/rds/fit_", nome, ".rds")
  fit$save_object(file = file_rds)

  risultati_loo[[nome]] <- fit$loo()
}

# comparison
confronto <- loo_compare(risultati_loo)
print(confronto, simplify = FALSE)
confronto_df <- as.data.frame(confronto)
confronto_df <- cbind(model = rownames(confronto_df), confronto_df)
# writexl::write_xlsx(confronto_df, "outputs/loo_comparison.xlsx")
# saveRDS(confronto, "outputs/loo_comparison.rds")
# saveRDS(risultati_loo, "outputs/loo_results.rds")

print(risultati_loo$mod_base)
print(risultati_loo$mod_inform)
print(risultati_loo$mod_regular)
print(risultati_loo$mod_strong)
print(risultati_loo$mod_wide)

# data prep for gt table
confronto <- readRDS("outputs/STAN/loo_comparison.rds")
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
  tab_style(
    style = cell_fill(color = "lightgrey", alpha = 0.5),
    locations = cells_body(rows = 1)
  ) %>%
  gtsave("outputs/STAN/loo_comparison.png")
