# QENS Fitter 4D 
## Description  

A plugin for the software package ([NeXpy](https://nexpy.github.io/nexpy/)). This plugin can be loaded into the window menu and used for semi-automated fitting of 4D-QENS data collected at time of flight instruments.



## Installation  

Clone the repo and install dependencies:

```bash
git clone https://github.com/yourusername/project-name.git
cd project-name
pip install .
```

## Usage

This plugin allows for analysis of 4D QENS data in the NeXus file format (signal shape taking the form of [E,L,K,H]). The files are fit by allowing users to first initialized initial conditions for fitting by fitting individual Q-points. 

These initial fitting conditions can then be passed to an automated mutli-threaded algorithm which raster-scans through the data fitting points independantly. 

Initial values for linewidths of one of the lorentzian functions can be supplied in the form of an array. Current work is being done to implement a tool to generate these models with current plans focussed on developing for the Chudley Elliott model in 3D. Any additional models are welcome to be suggested, preferably with example implementations.