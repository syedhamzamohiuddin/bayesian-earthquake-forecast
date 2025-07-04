# Bayesian Forecasting of Global Earthquakes (1918–2020)

This project uses Bayesian autoregressive models to forecast the annual count of worldwide earthquakes (magnitude > 7) using data from the USGS. The analysis focuses on a standalone AR(3) model with posterior inference, model order selection, and forecasting. The report also includes theoretical comparison with a mixture of AR(3) models.

## 🧠 Methods
- Model: AR(3) with Normal-Inverse Gamma conjugate prior
- Inference: Posterior sampling with 5000 draws in R
- Model order selection using PACF, AIC, and BIC
- Prior sensitivity analysis on hyperparameters

> 🔍 Note: While the accompanying report includes discussion of a mixture of AR(3) models and their DIC-based comparison, the current codebase implements only the standalone AR(3) analysis and forecasting pipeline.

## 📈 Results
- AR(3) model selected via PACF and AIC/BIC
- Posterior inference shown to be robust to prior hyperparameters
- Forecasts for 2021–2024: **12**, **10**, **10**, **11** earthquakes

## 📂 Project Structure
- `R/`: Scripts for model order selection, prior sensitivity analysis, and prediction
- `data/`: Earthquake count data (taken from USGS website)
- `results/`: Saved plots and intermediate outputs (PACF, AIC/BIC)
- `report/`: PDF capstone report written in Overleaf

## 📌 Key Takeaways
- Demonstrates practical application of Bayesian time series modeling in R
- Implements MCMC-based posterior estimation with uncertainty quantification
- Uses real-world geophysical data in a forecasting context

## 📚 Tools Used
- R (version 4.x)
- Key packages: `mvtnorm`, `MASS`, `stats`, `base`

Install required packages with:

```r
install.packages(c("mvtnorm", "MASS"))

## 🧩 Planned Improvements
- Add inline comments to R scripts for better readability
- Reproduce missing mixture model code (optional)
- Convert scripts into R Markdown notebooks (optional)
