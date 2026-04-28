## Population estimate with RMark, POPAN model
{
  # load packages
  library(readxl)
  library(lubridate)
  library(RMark)
  library(tidyverse)
}

{
  # function for time intervals (calculates difference in months)
  intervals_months <- function(dates) {
    diffs <- difftime(dates[-1], dates[-length(dates)], units = "days")
    as.numeric(diffs) / 30.44 # 1 month ≈ 30.44 days (mean Gregorian month)
  }
}

# Define the Master List of unique sampling dates across ALL areas
master_dates <- sort(unique(as.Date(c(
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
      valley = substr(area, 1, 1), # Extract 'S' o 'R'
      altitude = as.factor(substr(area, 2, 2)) # Extract '1', '2' o '3'
    )

  # check data
  str(all_data_ch)
  summary(all_data_ch$area) # how many individuals per area
}

# Data preparation for POPAN
{
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
  dormouse.ddl$pent$Season <- "active" # default value
  # identify hibernation intervals
  dormouse.ddl$pent$Season[
    dormouse.ddl$pent$par.index %in%
      c(12, 19, 34, 41, 56, 63, 78, 85, 100, 107, 122, 129)
  ] <- "hibernation"
  dormouse.ddl$pent$Season <- as.factor(dormouse.ddl$pent$Season)
  head(dormouse.ddl$Phi)
}


## ---- POPAN analysis ----

### ---- Fit models (optional) ----
# Please note: The models have already been fitted and the output files 
# (.res, .vcv, etc.) are saved in the directory: outputs/POPAN/models
# You can skip this block and load the existing results to save time.
# Re-run only if you need to refit the models from scratch (and on Windows laptop with Mark).
# {
#   model1 <- mark(
#     # full
#     dormouse.proc,
#     dormouse.ddl,
#     model.parameters = list(
#       Phi = list(formula = ~ Season * valley + altitude),
#       p = list(formula = ~time),
#       pent = list(formula = ~time),
#       N = list(formula = ~group)
#     ),
#     model.name = "full",
#     filename = paste0("POPAN_full")
#   )
#   model2 <- mark(
#     dormouse.proc,
#     dormouse.ddl,
#     model.parameters = list(
#       Phi = list(formula = ~ Season * valley + altitude),
#       p = list(formula = ~time),
#       pent = list(formula = ~Season),
#       N = list(formula = ~group)
#     ),
#     model.name = "full_pent_seasonal",
#     filename = paste0("POPAN_full_pent_seasonal")
#   )
#   model3 <- mark(
#     dormouse.proc,
#     dormouse.ddl,
#     model.parameters = list(
#       Phi = list(formula = ~ Season * valley + altitude),
#       p = list(formula = ~time),
#       pent = list(formula = ~1),
#       N = list(formula = ~group)
#     ),
#     model.name = "full_pent_constant",
#     filename = paste0("POPAN_full_pent_constant")
#   )
#   model4 <- mark(
#     dormouse.proc,
#     dormouse.ddl,
#     model.parameters = list(
#       Phi = list(formula = ~ Season * valley),
#       p = list(formula = ~time),
#       pent = list(formula = ~1),
#       N = list(formula = ~group)
#     ),
#     model.name = "Phi_season_valley",
#     filename = paste0("POPAN_Phi_season_valley")
#   )
#   model5 <- mark(
#     dormouse.proc,
#     dormouse.ddl,
#     model.parameters = list(
#       Phi = list(formula = ~ Season * altitude),
#       p = list(formula = ~time),
#       pent = list(formula = ~1),
#       N = list(formula = ~group)
#     ),
#     model.name = "Phi_season_altitude",
#     filename = paste0("POPAN_Phi_season_altitude")
#   )
#   model_null <- mark(
#     dormouse.proc,
#     dormouse.ddl,
#     model.parameters = list(
#       Phi = list(formula = ~1),
#       p = list(formula = ~1),
#       pent = list(formula = ~1),
#       N = list(formula = ~group)
#     ),
#     model.name = "null",
#     filename = paste0("POPAN_null")
#   )
#   model_basic_demo <- mark(
#     dormouse.proc,
#     dormouse.ddl,
#     model.parameters = list(
#       Phi = list(formula = ~Season),
#       p = list(formula = ~1),
#       pent = list(formula = ~time),
#       N = list(formula = ~group)
#     ),
#     model.name = "basic_demo",
#     filename = paste0("POPAN_basic_demo")
#   )
#   model_constant_p <- mark(
#     dormouse.proc,
#     dormouse.ddl,
#     model.parameters = list(
#       Phi = list(formula = ~Season),
#       p = list(formula = ~1),
#       pent = list(formula = ~time),
#       N = list(formula = ~group)
#     ),
#     model.name = "constant_p",
#     filename = paste0("POPAN_constant_p")
#   )
# }
# # Model comparison
# {
#   mod_comparison <- collect.models()
#   print(mod_comparison$model.table)
# }
# # save model comparison table
# # writexl::write_xlsx(mod_comparison$model.table,
# #                     "outputs/POPAN/POPAN_mod_comparison_table.xlsx")
# saveRDS(model4, "outputs/POPAN/models/POPAN_model4.rds")
# saveRDS(model5, "outputs/POPAN/models/POPAN_model5.rds")
# saveRDS(model3, "outputs/POPAN/models/POPAN_model3.rds")


##---- Results ----
{
  mod_comparison_table <- read_xlsx("outputs/POPAN/POPAN_mod_comparison_table.xlsx")
  print(mod_comparison_table)
}

# PICK THE BEST MODEL
{
  best_mod <- readRDS("outputs/POPAN/models/POPAN_model4.rds")
}

# helper function to plot abundance results
plot_popan_abundance <- function(
  chosen_model,
  model_description = "Single Best Model",
  is_averaged = FALSE
) {
  abund_data <- chosen_model$results$derived$N

  # Mapping Group <-> Area
  mapping <- arrange(distinct(dormouse.proc$data[, c("group", "area")]), group)
  area_names <- mapping$area
  n_areas <- length(area_names)
  n_occasions <- K_total

  # Data frame of results
  if (nrow(abund_data) != (n_areas * n_occasions)) {
    stop(
      "Error: The number of estimates in the model does not match Areas * Occasions."
    )
  }

  abund_clean <- data.frame(
    Area = rep(area_names, each = n_occasions),
    Occasion = rep(1:n_occasions, times = n_areas),
    Date = rep(master_dates, times = n_areas),
    Estimate = abund_data$estimate,
    SE = abund_data$se,
    LCL = abund_data$lcl,
    UCL = abund_data$ucl
  )

  # Faceted plot
  p <- ggplot(
    abund_clean,
    aes(x = Date, y = Estimate, group = Area, color = Area, fill = Area)
  ) +
    geom_ribbon(aes(ymin = LCL, ymax = UCL), color = NA, alpha = 0.2) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    facet_wrap(~Area, ncol = 3, scales = "free_y") +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(
      title = "Estimated Abundance of Muscardinus avellanarius",
      subtitle = paste("Model:", model_description),
      x = "Survey Date",
      y = "Estimated Population Size"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold")
    )

  return(list(plot = p, data = abund_clean))
}

# Results
{
  # Model summary
  print(summary(best_mod))
  res_best <- plot_popan_abundance(best_mod)
  print(res_best$plot)
  print(res_best$data)
}

# comparison with other models
model5 <- readRDS("outputs/POPAN/models/POPAN_model5.rds")
second_best <- plot_popan_abundance(
  model5,
  model_description = "Second best model"
)
print(second_best$plot)
model3 <- readRDS("outputs/POPAN/models/POPAN_model3.rds")
third_best <- plot_popan_abundance(
  model3,
  model_description = "Third best model"
)
print(third_best$plot)


##---- Model averaging ----
# Please note: The models averaging has already been performed.
# If you haven't re-run the models, you can skip this block and load the existing results to save time.
# Re-run only if you refit the models from scratch (and on Windows laptop with Mark).
# {
#   # define range of models for phi
#   Phi.full = list(formula = ~ Season * valley + altitude)
#   Phi.altitude = list(formula = ~ Season * altitude)
#   Phi.valley = list(formula = ~ Season * valley)
#   # for p
#   p.time = list(formula = ~time)
#   pent.const = list(formula = ~1)
#   N.group = list(formula = ~group)

#   # Run all pairings of models
#   dormouse.model.list = create.model.list("POPAN")
#   dormouse.results = mark.wrapper(
#     dormouse.model.list,
#     data = dormouse.proc,
#     ddl = dormouse.ddl,
#     delete = TRUE
#   )
#   dormouse.results
# }
# # saveRDS(dormouse.results, "outputs/POPAN_dormouse_results.rds")

dormouse.results <- readRDS("outputs/POPAN/POPAN_dormouse_results.rds")

# averaged N estimates per area
{
  N.estimates = model.average(
    dormouse.results,
    "derived",
    parameter = "N",
    vcv = TRUE
  )
  # calculate unique individuals (Mt) for each area
  Mt_per_sito <- aggregate(
    rep(1, nrow(all_data_ch)),
    by = list(Area = all_data_ch$area),
    FUN = sum
  )
  names(Mt_per_sito) <- c("Area", "Mt")

  # Order Mt according to model levels
  Mt_per_sito <- Mt_per_sito[
    match(levels(all_data_ch$area), Mt_per_sito$Area),
  ]

  # Extract model averaging results
  av_data <- N.estimates$estimates

  # Create data frame
  # RMark orders results by group, then by occasion
  mapping <- arrange(distinct(dormouse.proc$data[, c("group", "area")]), group)
  area_names <- mapping$area

  abund_avg <- N.estimates$estimates %>%
    left_join(Mt_per_sito, join_by("group" == "Area")) %>%
    mutate(
      Estimate = estimate + Mt,
      LCL = lcl + Mt,
      UCL = ucl + Mt
    )
  print(abund_avg)
}

# averaged trends
w4 <- dormouse.results$model.table$weight[1]
w5 <- dormouse.results$model.table$weight[2]
w3 <- dormouse.results$model.table$weight[3]

# sum weight to normalise them (sum to 1)
sum_weights <- w4 + w5 + w3
W4 <- w4 / sum_weights
W5 <- w5 / sum_weights
W3 <- w3 / sum_weights

# extract 'estimate' from the three models
N4 <- best_mod$results$derived$N$estimate
N5 <- model5$results$derived$N$estimate
N3 <- model3$results$derived$N$estimate

N_averaged_estimates <- (N4 * W4) + (N5 * W5) + (N3 * W3)

n_occasions = 23
abund_final <- data.frame(
  Area = rep(area_names, each = n_occasions),
  Occasion = rep(1:23, times = 6),
  Date = rep(master_dates, times = 6),
  Estimate = N_averaged_estimates
)

# SE_avg = sqrt( sum( W_i * (SE_i^2 + (N_i - N_avg)^2) ) )
SE4 <- best_mod$results$derived$N$se
SE5 <- model5$results$derived$N$se
SE3 <- model3$results$derived$N$se

abund_final$SE <- sqrt(
  W4 *
    (SE4^2 + (N4 - N_averaged_estimates)^2) +
    W5 * (SE5^2 + (N5 - N_averaged_estimates)^2) +
    W3 * (SE3^2 + (N3 - N_averaged_estimates)^2)
)

# approximate confidence intervale
abund_final$LCL <- abund_final$Estimate - (1.96 * abund_final$SE)
abund_final$UCL <- abund_final$Estimate + (1.96 * abund_final$SE)
# writexl::write_xlsx(abund_final, "outputs/POPAN/POPAN_abund_final.xlsx")


