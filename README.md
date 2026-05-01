# Bayesian Method for Causal Inference in Spatially-Correlated Time Series

## Overview
Assessing the causal impact of advertising campaigns across multiple retail locations is notoriously difficult due to spatial correlations between stores and low signal-to-noise ratios in sales data. This project introduces and implements a novel Bayesian approach utilizing a multivariate structural time series model to accurately identify and quantify causal effects under these challenging conditions. 

Unlike traditional methods that rely on point estimates, this approach leverages the full posterior distributions of latent variables, comparing observed data against synthetic counterfactuals to capture full uncertainty in the causal estimates.

## Key Features & Methodology

* **Multivariate Structural Time Series:** Models trends, seasonality, and regression effects while accounting for spatial dependencies between test stores using a G-Wishart prior on the precision matrix.
* **Synthetic Control Construction via DAEMVS:** Employs an adapted Deterministic Annealing Expectation-Maximization Variable Selection (DAEMVS) algorithm. This temperature-based approach efficiently selects the most relevant control stores from high-dimensional data, avoiding local optima.
* **Two-Stage Estimation Algorithm:**
  1. **EMVS Step:** Estimates sparse regression coefficients to de-trend the series.
  2. **MCMC & Kalman Filtering:** Samples parameters from posterior distributions (time-varying states, stationary constraints, and covariance matrices) using Markov Chain Monte Carlo methods, the Kalman filter, and Koopman's simulation smoother.
* **Stationarity Constraints:** Imposes a strict stationarity constraint on the local linear trend. This prevents prediction intervals from becoming overly wide during the intervention period, solving a common issue where true causal impacts are masked by expanding variance.
* **KS-Distance Causal Estimand:** Replaces traditional scalar estimands with the **one-sided Kolmogorov-Smirnov (KS) distance**. By comparing the empirical cumulative distribution functions (CDFs) of the observed trend against simulated counterfactual trends, the model rigorously determines statistical significance based on computed posterior predictive thresholds.

## Technology Stack
* **Language:** R
* **Core Libraries:** `Matrix`, `MCMCpack`, `mvtnorm`, `MASS`, `BDgraph`, `matrixStats`, `parallel`
* **Visualization:** `ggplot2`, `gridExtra`

## Results & Simulation Performance
Extensive simulation studies across various datasets demonstrate the superiority of the KS-distance estimand. In comparative testing:
* **Non-Stationary Models:** Produced overly expansive prediction intervals, failing to detect true causal impacts.
* **Stationary Models (Proposed):** Maintained moderate prediction intervals, allowing the KS-distance metric to accurately distinguish impacted datasets from non-impacted ones and avoid false positives.

## Repository Structure
* `main.R`: Contains the core implementation including the Kalman filter, Koopman's filter, EMVS variable selection, and the multivariate MCMC algorithm.
* `causal_impact.R`: Contains the causal impact wrapper and threshold generation logic using the one-sided KS distance.
* `simulation.R`: Code to generate spatially correlated synthetic datasets with varying degrees of causal impact and visualize the results.
