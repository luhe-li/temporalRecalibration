# Audiovisual temporal recalibration

This project develops models of audiovisual temporal recalibration, containing scripts for experiments, model fitting, GUI simulation, and plotting related to all the recalibration models in the [paper](https://elifesciences.org/reviewed-preprints/97765).

## Folder Structure

The project is organized into the following main folders:

- **atheoretical_models**: Define and compare four recalibration models that do not assume the recalibration mechanism. 
- **data**: Storage for all raw datasets from experiments.
- **fit_results**: Download the existing fitting results from [OSF](https://elifesciences.org/reviewed-preprints/97765) and place it under the main directory. You can generate all the figures without rerunning the analyses by using these results files.
- **figures**: Codes that generate all the figures in the paper.
- **recalibration_models**: Includes major models of recalibration. Subfolders contain different model scripts including:`caulnf_asym`, `caulnf_sym`,`heu_asym`, `heu_sym`,`trigger_asym`,`trigger_sym`.
  - The first part is the specific recalibration implementaion. 
  - The second part is the type of temporal precision: `asym` corresponds to the models with modality-specific precision, and `sym` refers to the models with modality-independent precision.
- **utils**: Utility scripts and helper functions. `latexTable` and `colorGradient` are borrowed for illustration purposes.

## Usage

1. **Model Fitting**: Use `fit_local.m` to run model fitting locally, which calls 
the function `fit_recal_model.m`. Alternatively, you can also run it on the cluster by directly calling this function. These scripts will call functions from the subfolders to process respective models. It requires the access to the data, utility functions and VBMC.

2. **Model Recovery**: Run model recovery of all models using the `model_recovery*.m` script series step by step. Also requires using the data and models stored in their respective subfolders. These three steps can also be grouped as one script.

3. **Parameter Recovery**: Use `param_recovery.m` to run model recovery of the causal-inference model. You can do parameter recovery on other models by simply changing the `currModelStr`.

4. **GUI**: If available, use the `modelGUI.mlapp` under each model subfolders to interact with parameters and see its impact on the model prediction. Make sure your path is in the same directory as GUI.

5. **Regenerate figures**: You can regenerate figures in the main text and appendices of the paper once the `fit_results` are under the main directory. 

## Prerequisites

- **Psychtoolbox**: Required for running the experiment codes.
- **VBMC**: Necessary for model fitting. You can find it [here](https://github.com/acerbilab/vbmc).

## Contact

For any questions or issues, please contact luhe.li@nyu.edu.
