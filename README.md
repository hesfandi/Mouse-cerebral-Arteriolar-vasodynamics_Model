# Mouse-cerebral-Arteriolar-vasodynamics_Model
A computational model of mouse cerebral arteriolar vasodynamics under resting state and functional hyperemia.
Overview

This repository contains MATLAB scripts for the in-silico simulation of mouse cerebral arterioles, as described in our article, "Systems Biology Analysis of Vasodynamics in Mouse Cerebral Arterioles During Resting State and Functional Hyperemia."

Abstract

In our earlier work ([Link to previous repository or publication]), we developed a coarsely segmented model of mouse cerebral vasculature, featuring a closed circulatory system and key morphological characteristics of mouse cerebrovasculature. In this repository, we extend that framework to simulate hemodynamic and vasodynamic interactions within the vessel network under both active and resting brain states.

Project Description

This repository includes hemodynamic-vasodynamic simulations:

Micro-scale (cellular-level): Penetrating arteriole (PA) model simulations.

Macro-scale: Simulations in a network of coupled arteriolar segments.

Repository Contents

This repository contains MATLAB scripts and associated data files for:

Micro-scale (cellular-level) penetrating arteriole simulations

Macro-scale coupled arteriolar network simulations

Usage Instructions

Software Requirements

MATLAB version: Tested on MATLAB 2020b (should be compatible with other MATLAB versions)

No additional MATLAB toolbox dependencies are required

Due to the computational intensity of simulations, running scripts on powerful computing resources is strongly recommended.

Folder Descriptions

Micro-scale/

This folder contains scripts to simulate:

Vasomotion: Spontaneous vessel diameter oscillations with no neuronal input.

Functional Hyperemia (FH): Simulation scripts where neurovascular coupling pathways are initiated after 150 seconds (once vasomotion-related oscillations dampen), simulating impulse or step neuronal activities.

Macro-scale/

This folder contains scripts for simulations in a network of macro-scale coupled arteriolar segments. Detailed algorithms corresponding to these simulations are fully described in the associated article.

Acknowledgments

Hemodynamic analysis code:

Adapted from the C program developed by Prof. Timothy W. Secomb’s research group ([Link to Timothy W. Secomb's website or repository]).

We sincerely thank Prof. Secomb's group for making these valuable computational tools publicly available.

SMC electrophysiology model:

Primarily adapted from Karlin et al. (2015), providing detailed equations and parametrization.

We thank the authors for their comprehensive model description and availability.

SMC contraction mechanism model:

Adapted from the comprehensive work by Tim David's research group.

Publicly available on their GitHub repository (https://github.com/BlueFern/OO-NVU/releases/tag/v2.1).

We acknowledge and thank them for providing access to this valuable resource.

In-vivo local field potential (LFP) and vasodynamics data:

Data used in simulations was sourced from [citation placeholder].

We deeply appreciate their provision of accurate concurrent electrical and optical measurements in awake mouse brains.

Citation

If you use this repository or its components in your research, please cite our paper:

Esfandi et al., Systems Biology Analysis of Vasodynamics in Mouse Cerebral Arterioles During Resting State and Functional Hyperemia, submitted to eLife (2024).

