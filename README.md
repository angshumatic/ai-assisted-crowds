# AI-Assisted Crowds

This repository contains the self-contained R script and data files needed to
reproduce the main manuscript analyses and figures for:

**Enhancing Crowd Accuracy in AI-Assisted Judgmental Forecasting: Theory and Experimental Evidence**

## Repository Structure

- `main.R`: main analysis script
- `data/`: cleaned experimental data
- `preregistration/`: preregistrations, postmortems, and postmortem scripts
- `figures/`: generated manuscript figures; ignored by Git
- `AI Advice Generation.pdf`: study material

## Required Files

Run `main.R` from the repository root. The script expects the following
files:

- `main.R`
- `data/qualtricsdata_baseline.csv`
- `data/qualtricsdata_tr_worse.csv`
- `data/qualtricsdata_tr_better.csv`

The three CSV files contain the cleaned Qualtrics response data for the
Baseline, Worse Advice, and Better Advice experimental conditions.

## Outputs

The script creates the main analysis objects in the R session, including:

- `responses_all`
- `participant_summary`
- `tasks_all`
- `desc_stats` (Table 1)
- `forecast_rmse` (Table 2)
- `forecast_ttests` (Table 2)

It also creates the manuscript figures:

- `gamma_contourplot.pdf` (Figure 3a)
- `benefit_contourplot.pdf` (Figure 3b)
- `scatterplot.pdf` (Figure 4a)
- `kappa_hist.pdf` (Figure 4b)
- `gamma_hist.pdf` (Figure 4c)

Cite as: Pal, Angshuman, Asa B. Palley, Ville A. Satopää. Enhancing Crowd Accuracy in AI-Assisted Judgmental Forecasting: Theory and Experimental Evidence. Available at SSRN 6045634.