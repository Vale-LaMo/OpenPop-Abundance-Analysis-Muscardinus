# OpenPop-Abundance-Analysis-Muscardinus

This repository contains the data and R scripts used to estimate the abundance of the hazel dormouse (*Muscardinus avellanarius*). The analysis compares a classic Frequentist approach (**POPAN** via `RMark`) with a **Hierarchical Bayesian Model** implemented in **Stan**.

## Repository Structure

- `data/`: Raw capture histories in Excel format.
- `models/`: Source code for the Stan models (`.stan`).
- `scripts/`: R scripts for data processing and analysis, numbered by execution order.
- `outputs/`: Model outputs, including saved Mark fits, Stan fits, and summary tables (`.xlsx`).

## Software Requirements

To reproduce the analyses, you will need:
1. **R** (version >= 4.0)
2. **Program MARK**: Required for the POPAN analyses ([download here](http://www.phidot.org/software/mark/)).
3. **Stan**: Required for Bayesian models via the `cmdstanr` package ([download here](https://mc-stan.org/)).

## Contact
Valentina La Morgia - [valentina.lamorgia@isprambiente.it]
Institute for Environmental Protection and Research (ISPRA), Italy
