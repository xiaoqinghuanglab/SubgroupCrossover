# Subgroup Analysis in Crossover Studies

This repository contains the R simulation code and SAS example application code accompanying the manuscript on subgroup analysis in crossover studies.

## Files

- `simulation.R`: Main program for the Monte Carlo simulation study.
- `functions.R`: Functions for fitting the CS, UN, and ANOVA models, including subgroup and naive models.
- `Truth_Table.xlsx`: Parameter settings used in the simulation scenarios.
- `NACC_example_code.sas`: Example SAS code for cohort selection, model fitting, and figure generation in the NACC application.

## Simulation

The simulation evaluates model performance across different sample sizes, covariance structures, sequence and subgroup allocation ratios, and treatment/carryover interaction scenarios.

Performance measures include bias, coverage probability, and power.

The R code requires the following packages:

`dplyr`, `tidyr`, `parallel`, `readxl`, `openxlsx`, `MASS`, and `mmrm`.

The simulation is designed for command-line execution using a task index:

```bash
Rscript simulation.R 1
```

File paths and computing-resource settings in `simulation.R` may need to be modified for use in a different computing environment.

## NACC Example

`NACC_example_code.sas` provides example code for the real-world application using data from the National Alzheimer's Coordinating Center (NACC).

The NACC data and the complete data-cleaning pipeline are not included in this repository. Therefore, this code is provided as an illustration of the cohort construction, statistical modeling, and figure generation used in the application rather than as a fully reproducible NACC analysis.

## Citation

If you use this code, please cite the accompanying manuscript when available.
