## Code structure **

### **Functions Folder**

This folder contains reproducible R scripts with modular functions:

- `functions_for_analysis.R`: Core functions for data analysis.
- `functions_for_process.R`: Functions used for data cleaning and processing.
- `functions_for_table_figure.R`: Functions to generate tables and visualizations.
- `unpack_trial_info.R`: Helper script to process field trial-level information.

### **Main Folder**

This folder contains scripts for the main data pipeline and analysis workflow:

- **`0_Set_up_preparation.R`**: Prepares the workspace, loads libraries, and sets up the environment.
- **`1.Data_Process.Rmd`**: Processes raw data, including cleaning and transformation.
- **`2.Data_Analysis.Rmd`**: Performs statistical analysis, modeling, and generates outputs.

## **How to reproduce R scripts**

1. Start with `0_Set_up_preparation.R` to initialize the environment.
2. Run `1.Data_Process.Rmd` to clean and process the data.
3. Use `2.Data_Analysis.Rmd` for analysis and result generation.

Please contact me at jaeseok2@illinois.edu if you need access to the data, figures, or tables.
