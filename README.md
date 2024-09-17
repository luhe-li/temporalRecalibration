# Audiovisual temporal recalibration

This project develops models of audiovisual temporal recalibration, containing scripts for experiments, model fitting, GUI simulation, and plotting related to all the recalibration models in the [paper](https://elifesciences.org/reviewed-preprints/97765).

## Folder Structure

The project is organized into the following main folders:

- **atheoretical_models**: Includes atheoretical models of recalibration.
- **data**: Storage for all datasets used in the project.
- **figures**: Codes that generate all the figures in the paper.
- **recalibration_models**: Includes major models of recalibration. Subfolders contain different model scripts such as:
  - `caulnf_asym`, `caulnf_sym`, `criteria_asym`, `fixed_asym`, `fixed_sym`, `heu_asym`, `heu_sym`: Specific model implementations. `asym` corresponds to the models with modality-specific precision, and `sym` refers to the models with modality-independent precision.
- **utils**: Utility scripts and helper functions.

## Usage

1. **Model Fitting**: Use `fit_local.m` to run model fitting locally, which calls 
the function `fit_recal_model_VBMC.m`. Alternatively, you can also run it on hpc by directly calling this function. These scripts will call functions from the subfolders to process respective models. It requires the access to the data, utility functions and VBMC.

2. **Model Recovery**: Run model recovery of all models using the `model_recovery*.m` script series. Also requires using the data and models stored in their respective subfolders. Use `plot_model_recovery.m` to visualize and inspect results.

3. **Parameter Recovery**: Use `param_recovery.m` to test the accuracy of parameter estimation of the causal-inference model. Use `plot_parameter_recovery.m` to visualize and inspect results.

4. **GUI**: Use the `modelGUI.mlapp` under each model subfolders to interact with parameters and see its impact on the model prediction.

## Prerequisites

- **Psychtoolbox**: Required for running the experiment codes.
- **VBMC**: Necessary for model fitting. You can find it [here](https://github.com/acerbilab/vbmc).

## Contact

For any questions or issues, please contact luhe.li@nyu.edu.
