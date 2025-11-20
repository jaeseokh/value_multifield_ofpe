# Project Structure Overview

This repository contains the full workflow for the **Value of Multi-Field OFPE** project.  
It includes reproducible code for:

- Tree-ensemble prediction models (2022 version)
- Updated Bayesian–hierarchical yield–N modeling (2025)
- APSIM-based prior construction and Bayesian updating
- Complete data-processing pipeline for soils, elevation, Daymet weather, and experimental inputs

The project is organized as follows.

---

## 1. Code Folder

The `Code/` directory contains all modeling scripts and processing workflows, divided into two major components.

---

### 1.1 **Main (Current Workflow)**  
Located at: `Code/main/`

This folder includes the **active 2025 pipeline**:

#### **`0_set_up.Rmd`**
Initializes the environment, loads libraries, sets project paths, and sources helper functions.

#### **`1.data_processing.Rmd`**
Full ETL pipeline:

- Merge raw OFPE plot-level datasets  
- Extract SSURGO soils (area-weighted)  
- Extract DEM and topographic variables  
- Generate **Daymet timing-specific weather features**  
- Assemble field-year analysis-ready datasets  

#### **`2a.tree_model_analysis.RMD`**
Analysis script for the updated 2022 tree-ensemble framework:  
Runs XGBoost and Random Forest models and evaluates predictability across unseen fields.

#### **`2b.tree_model_results.RMD`**
Summaries, cross-field diagnostics, and limitations of tree-based approaches.

#### **`3a.apsim_priors_ii.RMD`**
Generates APSIM-based prior moments (mean, variance, skewness) for yield given N and soil/weather combinations.

#### **`3b.bayesian_update_ii.RMD`**
Bayesian updating pipeline:

- Incorporates OFPE data  
- Produces posterior N-response curves  
- Computes posterior predictive yield distributions  

#### **`3c.tree_vs_bayesian.RMD`**
Model comparison between tree-ensemble and Bayesian posterior predictions.

#### **`Anonymize_Naming.RMD`**
Maps field/producer identifiers to anonymized labels for publication.

---

### 1.2 **SRC (Function Library)**  
Located at: `Code/src/`

Reusable modular R functions:

#### **`functions_for_processing.R`**
Core processing utilities:
- Daymet timing-specific feature construction  
- SSURGO soil extraction  
- DEM/topographic calculations  
- Geometry (sf) handling  

#### **`functions_tree_model.R`**
Wrappers for Random Forest, XGBoost, tuning, and prediction utilities.

#### **`functions_bayes_priors.R`**
Tools to extract APSIM simulation moments and build Bayesian priors.

#### **`functions_bayes_model.R`**
Posterior computation utilities:
- Likelihood  
- Prior transformation  
- Posterior integration  

All functions are automatically sourced by the scripts in `Code/main/`.

---

### 1.3 **Sandbox (Archive of Older Versions)**  
Located at: `Code/Sandbox/`

Contains the full 2022 pipeline for documentation and reproducibility.

#### **`Main(ver2022)/`**
Complete earlier workflow:
- Setup
- Data processing
- Analysis
- Figures and tables

#### **`Functions(ver2022)/`**
Archived versions of older function scripts.

---

## 2. Book Folder  
Located at: `book/`

Contains the Quarto manuscript:

- **Chapter 1** – Updated introduction & methodological motivation  
- **Part I** – 2022 tree-ensemble prediction  
- **Part II** – 2025 Bayesian–hierarchical analysis  
  (priors → update → posterior)

The book renders to HTML/PDF for dissemination.

---

## 3. Data Folder  
Located at: `Data/`

Contains:

- Raw OFPE experimental data  
- Cleaned/processed datasets  
- Extracted weather, soil, and DEM variables  
- APSIM simulation outputs  

**Note:**  
Data are not publicly distributed.  
For access, contact **jaeseok2@illinois.edu**.

---
