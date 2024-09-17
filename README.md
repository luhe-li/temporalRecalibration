# Audiovisual temporal recalibration

This project develops models of audiovisual temporal recalibration. It contains scripts for experiments, model fitting, simulation, GUI of major models, and plotting related to all the recalibration models in the [paper](https://elifesciences.org/reviewed-preprints/97765).

## Folder Structure

The project is organized into the following main folders:

- **atheoretical_models**: Includes atheoretical models of recalibration.
- **data**: Storage for all datasets used in the project.
- **figures**: Codes that generate all the figures in the paper.
- **recalibration_models_VBMC**: Includes major models of recalibration. Subfolders contain different model scripts such as:
  - `caulnf_asym`, `caulnf_sym`, `criteria_asym`, `fixed_asym`, `fixed_sym`, `heu_asym`, `heu_sym`: Specific model implementations.
  - **Scripts**:
    - `fit_local.m` and `fit_recal_model_VBMC.m`: Major scripts for fitting models, calling functions stored in the model subfolders.
    - `model_recovery.m`, `param_recovery.m`: Scripts for model recovery and parameter recovery.
    - `plot_confusion_matrix.m`, `plot_model_prediction.m`, `plot_model_recovery.m`: Generate different plots to inspect modeling results and save them to the respective folders.
    - `sim_*`: Scripts for simulating different aspects of the model and saving the results.
- **utils**: Utility scripts and helper functions.

## Main Scripts

- **Fitting Scripts**:
  - `fit_local.m`, `fit_recal_model_VBMC.m`: The main scripts for fitting recalibration models. These scripts call functions from different model subfolders (e.g., `cauInf_asym`, `fixed_sym`).
  
- **Model Recovery and Parameter Recovery**:
  - `model_recovery.m`, `model_recovery_s1.m`, `model_recovery_s2.m`: Scripts used for model recovery experiments.
  - `param_recovery.m`: Used to recover parameters from the model.

- **Plotting Scripts**:
  - `plot_*`: Scripts for plotting results, saving the figures under the corresponding folders in `figures`. Examples include `plot_confusion_matrix.m`, `plot_model_prediction.m`, and `plot_model_recovery.m`.

- **Simulation Scripts**:
  - `sim_*`: These scripts, such as `sim_cauInf_asym_reduced.m`, simulate various aspects of the model for analysis and testing.

## Usage

1. **Model Fitting**: Use `fit_local.m` to run model fitting locally, which calls 
the function `fit_recal_model_VBMC.m`. Alternatively, you can also run it on hpc by directly calling this function. These scripts will call functions from the subfolders to process models. It requires the access to the data, utility functions and VBMC.

2. **Model Recovery**: Run model recovery of all models using the `model_recovery*.m` script series. Also requires using the data and models stored in their respective subfolders. Use `plot_model_recovery.m` to visualize and inspect results.

3. **Parameter Recovery**: Use `param_recovery.m` to test the accuracy of parameter estimation of the causal-inference model. Use `plot_parameter_recovery.m` to visualize and inspect results.

4. **GUI**: Use the `modelGUI.mlapp` under each model subfolders to interact with parameters and see its impact on the model prediction.


## Prerequisites

- **Psychtoolbox**: Required for running the experiment codes.
- **VBMC**: Necessary for model fitting. You can find it [here](https://github.com/acerbilab/vbmc).

## Contact

For any questions or issues, please contact luhe.li@nyu.edu.