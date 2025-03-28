# Mouse-cerebral-Arteriolar-vasodynamics_Model

## About
A computational model of mouse cerebral arteriolar vasodynamics under resting state and functional hyperemia.

## Overview
This repository contains MATLAB scripts for the in-silico simulation of mouse cerebral arterioles, as described in our article, **"Systems Biology Analysis of Vasodynamics in Mouse Cerebral Arterioles During Resting State and Functional Hyperemia."**

## Abstract
In our earlier work ([link to previous repository or publication](https://github.com/your_previous_repo)), we developed a coarsely segmented model of mouse cerebral vasculature, featuring a closed circulatory system and key morphological characteristics of mouse cerebrovasculature. In this repository, we extend that framework to simulate hemodynamic and vasodynamic interactions within the vessel network under both active and resting brain states.

## Project Description
This repository includes hemodynamic-vasodynamic simulations:
- **Micro-scale (cellular-level)**: Penetrating arteriole (PA) model simulations.
- **Macro-scale**: Simulations in a network of coupled arteriolar segments.

## Repository Contents
This repository contains MATLAB scripts and associated data files for:
- **Micro-scale simulations** (cellular-level penetrating arteriole model)
- **Macro-scale simulations** (network of macro-scale coupled arteriolar segments)

## Usage Instructions

### Software Requirements
- MATLAB version: Tested on MATLAB 2020b (should be compatible with other MATLAB versions)
- No additional MATLAB toolbox dependencies required.
- Simulations are computationally intensive; powerful computing resources are strongly recommended.

### Folder Descriptions

#### `Micro-scale/`
Scripts for simulations at the cellular-level PA model:
- **Vasomotion**: Spontaneous oscillations without neuronal input.
- **Functional Hyperemia (FH)**: Simulations with neurovascular coupling pathways initiated after 150 seconds (once vasomotion oscillations dampen), simulating impulse or step neuronal activities.

#### `Macro-scale/`
Scripts for simulations in a network of macro-scale coupled arteriolar segments. The algorithms corresponding to these simulations are fully described in the associated article.

## Acknowledgments

- **Hemodynamic Analysis Code**
  - Adapted from the C program developed by Prof. Timothy W. Secomb’s research group ([Secomb's Lab Website](http://www.secomb.org)).
  - We sincerely thank Prof. Secomb’s group for making these valuable computational tools publicly available.

- **SMC Electrophysiology Model**
  - Adapted primarily from Karlin et al. (2015), providing detailed equations and parameterization.
  - We appreciate their comprehensive descriptions and availability.

- **SMC Contraction Mechanism Model**
  - Adapted from work by Tim David’s research group.
  - Available on their GitHub repository ([OO-NVU Model Repository](https://github.com/BlueFern/OO-NVU/releases/tag/v2.1)).
  - We acknowledge their generosity in providing public access to this model.

- **In-vivo LFP and Vasodynamics Data**
  - Data used in simulations sourced from [Author et al., 20XX](https://doi.org/your-doi-link).
  - We appreciate their provision of accurate concurrent electrical and optical measurements in awake mouse brains.

## Citation
If you use this repository or its components in your research, please cite our paper:

> Esfandi et al., Systems Biology Analysis of Vasodynamics in Mouse Cerebral Arterioles During Resting State and Functional Hyperemia, submitted to eLife (2024).

---

For questions, contributions, or issues, please use the GitHub issue tracker or contact the authors directly.
