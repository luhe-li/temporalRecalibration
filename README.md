# Audiovisual temporal recalibration

This project develops computational models of audiovisual temporal recalibration. It includes scripts for:  
- **Experiments**  
- **Model fitting** (Bayesian optimization and parameter estimation)  
- **GUI-based simulations**  
- **Plotting utilities** for reproducing figures in the [associated paper](https://elifesciences.org/reviewed-preprints/97765).  

## Installation & Setup

### 1. Download dependencies
- Install [Psychtoolbox](https://psychtoolbox.org/) for experiment scripts.  
- Install [VBMC](https://github.com/acerbilab/vbmc) (MATLAB toolbox for Bayesian model fitting). 

### 2. Download precomputed results
- Download `fit_results.zip` from [OSF](https://osf.io/8s7qv/) and unzip it into the project root.  

### 3. Add paths in MATLAB
- Add the project folders to your MATLAB path using `addpath(genpath('temporalRecalibration'))`.  

## Usage

### 1. Model Fitting
- **Local fitting**: Run the `run_local.m` script under `recalibration_models/` to fit all recalibration models or a specific one.
- **Cluster fitting**: Modify script to run `fit_recal_model.m` to submit jobs via your cluster’s scheduling system.
- Recalibraton models:
   - `caulnf_asym/` Causal inference model, modality-specific temporal precision
   - `caulnf_sym/` Causal inference model, modality-independent temporal precision
   - `caulnf_asym_update/` A variant of causal inference model, modality-specific temporal precision, mentioned in Appendix 12
   - `heu_asym/` Asynchrony-contingent model, modality-specific temporal precision
   - `heu_sym/` Asynchrony-contingent model, modality-independent temporal precision
   - `trigger_asym/` Asynchrony-correction model, modality-specific temporal precision
   - `trigger_sym/` Asynchrony-correction model, modality-independent temporal precision

### 2. Model Recovery
Run the `model_recovery_s1.m`, `model_recovery_s2.m`, and `model_recovery_s3.m` scripts sequentially. Alternatively, modify and combine them to create a full pipeline. 

### 3. Parameter Recovery
Use the `param_recovery.m` script to run parameter recovery for a specific model. Modify the `currModelStr` variable to test other models.  

### 4. GUI Simulation
Launch the `modelGUI.mlapp` from the respective model subfolder if it exists (e.g., `recalibration_models/caulnf_asym`). Ensure your MATLAB working directory matches the model folder.

### 5. Regenerate Figures
Run the plotting scripts in the `figures` folder (e.g., `fig2.m`) to regenerate figures from the paper. Requires the `fit_results` folder in the root directory.  

## Troubleshooting
- **Path errors**: Ensure all folders (especially `fit_results` and `data`) are in the root directory.  
- **Missing dependencies**: Confirm Psychtoolbox and VBMC are installed and added to MATLAB’s path.  
- **GUI issues**: Use MATLAB ≥2021b for compatibility with `.mlapp` files.  

## Contact
For questions or bugs, contact [Luhe Li](mailto:luhe.li@nyu.edu).  
*Include error logs and MATLAB version when reporting issues.*  