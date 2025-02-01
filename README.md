# Audiovisual temporal recalibration

This project develops computational models of audiovisual temporal recalibration. It includes scripts for:  
- **Experiments**  
- **Model fitting** (Bayesian optimization and parameter estimation)  
- **GUI-based simulations**  
- **Plotting utilities** for reproducing figures in the [associated paper](https://elifesciences.org/reviewed-preprints/97765).  

---

## Folder Structure
root/
├── atheoretical_models/ # Models without assumptions about recalibration mechanisms
├── data/ # Raw experimental datasets (Matlab formats)
├── fit_results/ # Precomputed fitting results (download from [OSF](https://osf.io/8s7qv/))
├── figures/ # Scripts to regenerate all paper figures (Figs 1–6, Appendix 1-13)
├── recalibration_models/ # Core model implementations
│ ├── caulnf_asym/ # Causal inference model, modality-specific temporal precision
│ ├── caulnf_sym/ # Causal inference model, modality-independent temporal precision
│ ├── caulnf_asym_update/ # A variant of causal inference model, modality-specific temporal precision, mentioned in Appendix 12
│ ├── heu_asym/ # Asynchrony-contingent model, modality-specific temporal precision
│ ├── heu_sym/ # Asynchrony-contingent model, modality-independent temporal precision
│ ├── trigger_asym/ # Asynchrony-correction model, modality-specific temporal precision
│ └── trigger_sym/ # Asynchrony-correction model, modality-independent temporal precision
└── utils/ # Helper functions

---

<<<<<<< Updated upstream
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
=======
## Installation & Setup

1. **Download dependencies**:  
   - Install [Psychtoolbox](https://psychtoolbox.org/) for experiment scripts.  
   - Install [VBMC](https://github.com/acerbilab/vbmc) (MATLAB toolbox for Bayesian model fitting).  
   - Clone this repository:  
     ```bash  
     git clone https://github.com/your-repo/audiovisual-recalibration.git  
     ```  

2. **Download precomputed results**:  
   - Download `fit_results.zip` from [OSF](https://osf.io/8s7qv/) and unzip it into the project root.  

3. **Add paths in MATLAB**:  
   - Add the project folders to your MATLAB path using `addpath(genpath('temporalRecalibration'))`.  

---

## Usage

### 1. Model Fitting
- **Local fitting**: Run the `run_local.m` script for all recalibration models or a specific one.
- **Cluster fitting**: Modify script to run `fit_recal_model.m` to submit jobs via your cluster’s scheduling system.

### 2. Model Recovery
Run the `model_recovery_s1.m`, `model_recovery_s2.m`, and `model_recovery_s3.m` scripts sequentially. Alternatively, modify and combine them to create a full pipeline. 

### 3. Parameter Recovery
Use the `param_recovery.m` script to run parameter recovery for a specific model. Modify the `currModelStr` variable to test other models.  

### 4. GUI Simulation
Launch the `modelGUI.mlapp` from the respective model subfolder if it exists (e.g., `recalibration_models/caulnf_asym`). Ensure your MATLAB working directory matches the model folder.

### 5. Regenerate Figures
Run the plotting scripts in the `figures` folder (e.g., `fig2.m`) to regenerate figures from the paper. Requires the `fit_results` folder in the root directory.  
>>>>>>> Stashed changes

---

## Troubleshooting
- **Path errors**: Ensure all folders (especially `fit_results` and `data`) are in the root directory.  
- **Missing dependencies**: Confirm Psychtoolbox and VBMC are installed and added to MATLAB’s path.  
- **GUI issues**: Use MATLAB ≥2021b for compatibility with `.mlapp` files.  

---

## Contact
For questions or bugs, contact [Luhe Li](mailto:luhe.li@nyu.edu).  
*Include error logs and MATLAB version when reporting issues.*  