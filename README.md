# Random Perturbation Subsampling for Composite Quantile Regression with Massive Data

This repository contains the code used for the simulation and real data studies presented in the paper *"Perturbation Subsampling Estimation for Diverging-Dimensional Composite Quantile Regression."*

## Structure

### 1. Simulation Code (`simulation code`)

The folder `simulation code` contains all the code for the simulations in the paper. It includes three subfolders corresponding to the simulation steps in the paper:

- **1 sensitivity analysis**: Code related to the selection of appropriate parameter settings in Section 4.1.
- **2 comparison of several subsampling methods**: Code for comparing different subsampling methods in Section 4.2.
- **3 running time**: Code for calculating the running time in Section 4.2.

### 2. Real Data Code (`real data code`)

The folder `real data code` contains the code for the empirical analysis, organized into two subfolders:

- **physicochemical properties**: Code for calculating the Prediction Errors (PEs) for the physicochemical properties dataset, corresponding to the empirical results in Figure 8.
- **superconductivty**: Code for calculating the Prediction Errors (PEs) for the superconductivty dataset, corresponding to the empirical results in Figure 8.

Additionally, the empirical data files are named `CASP.csv` and `superconductivty_train.csv` respectively, which support the analysis results.

## Usage

Each folder contains code and related data files for reproducing the results presented in the paper. Please refer to the specific folder corresponding to the simulation or empirical analysis you are interested in.
