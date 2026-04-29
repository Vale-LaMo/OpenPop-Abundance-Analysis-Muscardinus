# OpenPop-Abundance-Analysis-Muscardinus

This repository contains the data and the main R scripts used to estimate the abundance of the hazel dormouse (*Muscardinus avellanarius*) in two Alpine valleys. The analysis compares a classic Frequentist approach (**POPAN** via `RMark`) with a **Hierarchical Bayesian Model** implemented in **Stan**.

## Repository Structure

- `data/`: Raw capture histories in Excel format.
- `models/`: Source code for the Stan models (`.stan`).
- `scripts/`: R scripts for data processing and analysis, numbered by execution order.
- `outputs/`: Model outputs, including saved Mark fits, Stan fits, and summary tables (`.xlsx`).

## Instructions

To obtain the results of the POPAN analyses, run `scripts/01_abundance_grouped_POPAN.R`. For the POPAN analyses ran separately for the two study areas R1948 and S1966, run `scripts/01_abundance_single_POPAN.R`.    
For the tables showing detailed results of the POPAN analyses, see `scripts/07_supplementary_mark_tables.R`: tables are generated using the Mark output saved in `outputs/POPAN/models` (`*.out` files).

To obtain the results of the Bayesian analyses, run `scripts/04_abundance_hierachical_STAN_modelselection.R` and `scripts/05_abundance_hierachical_STAN_results.R`. For the Bayesian analyses ran separately for the two study areas R1948 and S1966, run `scripts/02_abundance_single_STAN.R`. Plots for the prior sensitivity checks and prior-posterior overlap are produced by `scripts/06_abundance_hierarchical_STAN_checks.R`. 

## Software Requirements

To reproduce the analyses, you will need:
1. **R** (version >= 4.0)
2. **Program MARK**: Required for the POPAN analyses ([download here](http://www.phidot.org/software/mark/)).
3. **Stan**: Required for Bayesian models via the `cmdstanr` package ([download here](https://mc-stan.org/)).

## Contact
Valentina La Morgia - [valentina.lamorgia@isprambiente.it]      
Institute for Environmental Protection and Research (ISPRA), Italy
