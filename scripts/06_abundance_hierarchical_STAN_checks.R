{
  library(cmdstanr)
  library(posterior)
  library(bayesplot)
  library(corrplot)
  library(vegan)
  library(sf)
  library(patchwork)
  library(tidybayes)
  library(bayestestR)
  library(overlapping)
}

rm(list = ls())
source("scripts/03_abundance_hierachical_STAN_data_prep.R")
fit <- readRDS("outputs/STAN/rds/fit_mod_inform.rds")


# ---- Prior posterior overlap for survival ----

n_samples <- 12000 # Total number of samples (3000 iterations * 4 chains)
# temporal intervals
# delta[1] ~15 gg, delta[11] winter (~210 gg)
d_breve <- 15 
d_lungo <- 210

# generate phi prior transforming the mu prior ~ HalfNormal(0,1)
prior_mu <- abs(rnorm(n_samples, 0, 1))
prior_phi_breve <- exp(-prior_mu * (d_breve/30.44)) # normalised on a monthly basis (cf. intervals)
prior_phi_lungo <- exp(-prior_mu * (d_lungo/30.44))
post_phi_breve <- as.vector(fit$draws("phi[1]"))
post_phi_lungo <- as.vector(fit$draws("phi[18]"))

get_mode <- function(x) {
  dens <- density(x)
  return(dens$x[which.max(dens$y)])
}
picco_breve <- get_mode(post_phi_breve)
picco_lungo <- get_mode(post_phi_lungo)
cat("Picco sopravvivenza breve (~15gg):", round(picco_breve, 3), "\n")
cat("Picco sopravvivenza invernale (~7 mesi):", round(picco_lungo, 3), "\n")

post_N <- as.vector(fit$draws("N[6,11]"))
df_phi_short <- data.frame(
  Value = c(prior_phi_breve, post_phi_breve),
  Type = rep(c("Prior", "Posterior"), each = length(post_N))
)
ggplot(df_phi_short, aes(x = Value, fill = Type)) +
  geom_density(alpha = 0.5) +
  labs(#title = "Overlap for Survival (phi) - short interval",
       subtitle = "Short summer interval",
       x = "Survival", y = "Density") +
  theme_minimal() -> df_phi_short

df_phi_long <- data.frame(
  Value = c(prior_phi_lungo, post_phi_lungo),
  Type = rep(c("Prior", "Posterior"), each = length(post_N))
)

ggplot(df_phi_long, aes(x = Value, fill = Type)) +
  geom_density(alpha = 0.5) +
  labs(#title = "Overlap for Survival (phi) - long interval",
       subtitle = "Long winter interval)",
       x = "Survival", y = "Density") +
  theme_minimal() -> df_phi_long

idx_summer <- c(1:10, 12:17, 19:22)
idx_winter <- c(11, 18)
mu_prior_samples <- abs(rnorm(n_samples, 0, 1))
prior_phi_summer_all <- as.vector(sapply(c(15), function(d) exp(-mu_prior_samples * (d/30.44))))
prior_phi_winter_all <- as.vector(sapply(c(210), function(d) exp(-mu_prior_samples * (d/30.44))))
post_phi_summer_all <- as.vector(fit$draws(paste0("phi[", idx_summer, "]")))
post_phi_winter_all <- as.vector(fit$draws(paste0("phi[", idx_winter, "]")))

get_mode <- function(x) {
  dens <- density(x)
  return(dens$x[which.max(dens$y)])
}

picco_summer <- get_mode(post_phi_summer_all)
picco_winter <- get_mode(post_phi_winter_all)
ov_summer <- overlapping::overlap(list(prior_phi_summer_all, post_phi_summer_all))$OV
ov_winter <- overlapping::overlap(list(prior_phi_winter_all, post_phi_winter_all))$OV

# --- paper fig ---
cat("--- Risultati Aggregati Sopravvivenza (phi) ---\n")
cat("ESTATE - Picco:", round(picco_summer, 3), "| Overlap Prior-Post:", round(ov_summer, 3), "\n")
cat("INVERNO - Picco:", round(picco_winter, 3), "| Overlap Prior-Post:", round(ov_winter, 3), "\n")

df_summer <- data.frame(
  Value = c(
    sample(prior_phi_summer_all, length(post_phi_summer_all), replace = TRUE), 
    post_phi_summer_all
  ),
  Type = rep(c("Prior", "Posterior"), each = length(post_phi_summer_all))
)

df_winter <- data.frame(
  Value = c(
    sample(prior_phi_winter_all, length(post_phi_winter_all), replace = TRUE), 
    post_phi_winter_all
  ),
  Type = rep(c("Prior", "Posterior"), each = length(post_phi_winter_all))
)

p1 <- ggplot(df_summer, aes(x = Value, fill = Type)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = picco_summer, linetype = "dashed") +
  labs(subtitle = paste0("Summer Aggregate (OV: ", round(ov_summer, 2), ")"),
       x = "Survival (phi)", y = "Density") +
  theme_minimal()

p2 <- ggplot(df_winter, aes(x = Value, fill = Type)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = picco_winter, linetype = "dashed") +
  labs(subtitle = paste0("Winter Aggregate (OV: ", round(ov_winter, 2), ")"),
       x = "Survival (phi)", y = "Density") +
  theme_minimal()

p1 + p2 + 
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Aggregate Prior-Posterior Overlap for Seasonal Survival")

ggsave(
    filename = "figs/Prior_Post_overlap_phi.tiff", # o .png, .pdf, .eps
    plot = last_plot(),                          # esporta l'ultimo grafico visualizzato
    device = "tiff",                             # formato file
    compression = "lzw",
    width = 180,                                 # larghezza (es. 180mm è lo standard "full page width")
    height = 120,                                # altezza in mm
    units = "mm",                                # unità di misura
    dpi = 300,                                   # risoluzione richiesta
    bg = "white"                                 # sfondo bianco (evita trasparenze indesiderate)
  )



# ---- Sensitivity checks ----

fit_mod_base <- readRDS("outputs/STAN/rds/fit_mod_base.rds")
fit_mod_inform <- readRDS("outputs/STAN/rds/fit_mod_inform.rds")
fit_mod_regular <- readRDS("outputs/STAN/rds/fit_mod_regular.rds")
fit_mod_strong <- readRDS("outputs/STAN/rds/fit_mod_strong.rds")
fit_mod_wide <- readRDS("outputs/STAN/rds/fit_mod_wide.rds")
modelli <- list(
  "Base" = fit_mod_base,
  "Informative" = fit_mod_inform,
  "Regular" = fit_mod_regular,
  "Strong" = fit_mod_strong,
  "Wide" = fit_mod_wide
)

df_forest <- lapply(names(modelli), function(nome) {
  modelli[[nome]]$summary(c("mu_active", "mu_hibernation", "alpha_p", "beta")) %>%
    mutate(Prior = nome)
}) %>% bind_rows()

df_forest <- df_forest %>%
  mutate(variable = gsub("\\[.*\\]", "", variable)) # Rimuove gli indici [1,1] per raggruppare i tipi di parametro

df_forest_clean <- df_forest %>%
  mutate(param_type = gsub("\\[.*\\]", "", variable)) %>%
  group_by(Prior, param_type) %>%
  summarise(
    mean_global = mean(mean),
    q5_global = mean(q5),
    q95_global = mean(q95),
    .groups = "drop"
  )

# PAPER FIG
ggplot(df_forest_clean, aes(x = mean_global, y = Prior, color = Prior)) +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = q5_global, xmax = q95_global), height = 0.3, size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~param_type, scales = "free_x") +
  labs(
    title = "Prior Sensitivity Check (Global Trends)",
    subtitle = "Aggregated posterior means across different prior specifications",
    x = "Posterior Mean (Logit scale for alpha_p, mortality scale for mu)",
    y = NULL
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "none", strip.text = element_text(face = "bold", size = 12))

ggsave(
    filename = "figs/Prior_Sensitivity_checks.tiff", # o .png, .pdf, .eps
    plot = last_plot(),                          # esporta l'ultimo grafico visualizzato
    device = "tiff",                             # formato file
    compression = "lzw",
    width = 180,                                 # larghezza (es. 180mm è lo standard "full page width")
    height = 120,                                # altezza in mm
    units = "mm",                                # unità di misura
    dpi = 300,                                   # risoluzione richiesta
    bg = "white"                                 # sfondo bianco (evita trasparenze indesiderate)
  )

