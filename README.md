# Mouse-Cerebral-Arteriolar-Vasodynamics_Model

## About  
A computational model of mouse cerebral arteriolar vasodynamics under resting state and functional hyperemia.

## Overview  
This repository contains MATLAB scripts for the in-silico simulation of vasodynamics in mouse cerebral arterioles, as described in our article, **"Systems Biology Analysis of Vasodynamics in Mouse Cerebral Arterioles During Resting State and Functional Hyperemia."**

## Abstract  
In our earlier work ([Cerebrovascular_Model](https://github.com/hesfandi/Cerebrovascular_Model)), we developed a coarsely segmented model of mouse cerebral vasculature, featuring a closed circulatory system and key morphological characteristics of mouse cerebrovasculature. In this repository, we extend that framework to simulate hemodynamic and vasodynamic interactions within the vessel network under both active and resting brain states.

## Repository Contents  
This repository contains MATLAB scripts for hemodynamic-vasodynamic simulations:  
- **Micro-scale**: Simulating vasodynamics in a cellular-level model of a penetrating arteriole (PA)  
- **Macro-scale**: Simulating vasodynamics in a network of coupled arteriolar segments  
 

### Software Requirements  
- MATLAB version: Tested on MATLAB 2020b (should be compatible with other MATLAB versions)  
- No additional MATLAB toolbox dependencies required  
- Simulations typically require 10–30 minutes on an up-to-date personal computer  

### Folder Descriptions  

#### `Micro-scale/`  
- **Vasomotion**: Scripts beginning with `Run_Vaso` simulate dampened oscillations in a cellular-level PA under varied mechanotransduction dynamics  
- **Functional Hyperemia (FH)**:  
  - Scripts beginning with `Run_IFH` simulate instantaneous neuronal activity (impulse response)  
  - Scripts beginning with `Run_SFH` simulate sustained neuronal activity (step response)  

#### `Macro-scale/`  
- `Run_Macro_Vaso` simulates undampened oscillations in a macro-scale model of coupled arteriolar segments  
- Scripts beginning with `Run_Macro_NVC` simulate modulation of the myogenic response by in-vivo recorded LFP signals to predict NVC-mediated arteriolar vasodynamics  

## Acknowledgments  

- **Hemodynamic Analysis Code**  
  - Adapted from the C program developed by Prof. Timothy W. Secomb’s research group ([Secomb's Lab Website](https://sites.arizona.edu/secomb/network-hemodynamics/))  
  - We sincerely thank Prof. Secomb’s group for making these valuable computational tools publicly available  

- **SMC Electrophysiology Model**  
  - Adapted primarily from Karlin et al. (2015), providing detailed equations and parameterization  
  - We appreciate their comprehensive modeling approach  

- **SMC Contraction Mechanism Model**  
  - Adapted from work by Tim David’s research group  
  - Available from their GitHub repository ([OO-NVU Model Repository](https://github.com/ModelDBRepository/237604))  

- **In-vivo LFP and Vasodynamics Data**  
  - Data used in simulations sourced from ([Mateo et al., 2017](https://www.cell.com/neuron/fulltext/S0896-6273(17)30980-7))  
  - We appreciate their provision of accurate concurrent electrical and optical measurements in awake mouse brains  

## Citation  
If you use this repository or its components in your research, please cite our paper:

> Esfandi et al., *Systems Biology Analysis of Vasodynamics in Mouse Cerebral Arterioles During Resting State and Functional Hyperemia*, submitted to eLife (2024)

---

For questions, contributions, or issues, please use the GitHub issue tracker or contact the authors directly.
