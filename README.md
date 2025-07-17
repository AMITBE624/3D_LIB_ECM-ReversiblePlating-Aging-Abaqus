# 3D_LIB_ECM-ReversiblePlating-Aging-Abaqus
## Overview
This repository contains Abaqus input files and user subroutines developed to simulate electro-chemo-mechanical (ECM) models for lithium-ion batteries. The models focus on analyzing reversible lithium plating, swelling force evolution, and aging behavior under mechanical constraint, as presented in the related publication mentioned below.

## Contents
- `input_files/` – Abaqus input files for constrained pouch cell simulations.
- `subroutines/` – User-defined subroutines for ECM model.
- `data/` – Experimental and simulation data.
- `figures/` – Visualizations related to the simulation results.
- `docs/` – Additional documentation and paper references.

## Features
- Coupled electrochemical and mechanical modeling.
- Simulation of lithium plating and stripping.
- Swelling force prediction during battery cycling.
- Aging characterization under mechanical constraints.

## Getting Started
1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/3D_LIB_ECM-ReversiblePlating-Aging-Abaqus.git
2. Compile and run user subroutines according to Abaqus documentation.
   ```bash
   @call "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2019.5.281\windows\bin\compilervars" intel64
   abaqus -j 3D_ECHMECH_Model_v12.inp -user usersublplfld.f -cpus 20 -int
4. Analyze output data using your preferred tools.

## Requirements

- Abaqus CAE (version 2025 or later)
- Fortran compiler compatible with Abaqus
- Python or MATLAB (optional, for data analysis)

## Citation
If you use this repository in your research, please cite:

Author(s). (Year). Electro-Chemo-Mechanical Coupling Unlocks New Insights into Lithium Plating and Aging Upon Mechanical Stress. Journal/Conference Name, Volume(Issue), pages. DOI: [insert DOI here]

Example:
Bhowmick, A., et al. (2025). Electro-Chemo-Mechanical Coupling Unlocks New Insights into Lithium Plating and Aging Upon Mechanical Stress. -------,---, ---. DOI: 10.xxxx/xxxxx

## Contact Information
- Research Group: 
Energy Mechanics and Sustainability Laboratory (EMSLab)

University of Delaware

Newark, Delaware, 19716
- Principal Investigator:
Prof. Jun Xu

Email: junxu@udel.edu
- Postdoctoral Researcher / Maintainer:
Dr. Amit Bhowmick

Email: amitbhowmick555@gmail.com

## Sponsor:
Startup funding of the University of Delaware and the research project with Dassault Systèmes
- Contact information (Abaqus)

Dr. Youngwon HAHN

Email: Youngwon.HAHN@3ds.com

