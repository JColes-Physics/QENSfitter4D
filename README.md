# QENS fitter 4D 
## Description  

A plugin for the software package ([NeXpy](https://nexpy.github.io/nexpy/)). This plugin can be loaded into the window menu and used for semi-automated fitting of 4D-QENS data collected at time-of-flight direct geometry instruments. This software has been built on the foundations developed in one of the first ([4D-QENS papers](https://doi.org/10.1103/b93d-755s)).

As of Tues Jul 21st 2026, this software is still in beta early developement. Developement is limited to available datasets for testing this software. Therefore, expect possible gaps in capability while necessary features are discovered. Further collaboration and developements are welcome and encouraged from the greater single crystal QENS community. Any interest in contributing or collaboration can reach out to Jared Coles (jared.coles@chalmers.se)

## Installation  

Clone the repo and install dependencies:

```bash
git clone https://github.com/JColes-Physics/QENSfitter4D.git
cd QENSfitter4D
pip install .
```

## Usage

This plugin allows for analysis of 4D QENS data in the NeXus file format (signal shape taking the form of [E,L,K,H]). 

Single Q-point fits can be performed using the "Initialize Fitting" option under the QENS tab in the pluggin bar. By selecting "Test Fits" one can select initial fitting parameters which can be stored in an NXprocess. Plotting options are available to verify the quality of fits given a set of initial conditions.

These initial fitting conditions, once stored, can then be passed to an automated mutli-threaded algorithm which raster-scans through the data fitting points independantly. This can be done by selecting the "Raster-Scan Fits" option on the "Initialize Fitting" window. From here you can select the Q-range and Q-steps you want to automatically fit over.

Initial values for linewidths of one of the lorentzian functions can be supplied in the form of an array. Current work is being done to implement a tool to generate these models with current plans focussed on developing for the Chudley Elliott model in 3D. This will be built into the "Second Voigt Model" option under the QENS tab in the pluggin bar and will have the option to save the chosen model under any entry of your choosing for ease of loading. Any additional models are welcome to be suggested, preferably with example implementations.