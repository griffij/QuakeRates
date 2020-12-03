# QuakeRates
Code for analysis of long-term earthquake occurrence rates

This code is designed to analyse earthquake inter-event times from palaeoseismology data. Inputs are expected to be probability distributions of earthquake occurrence times, such as the output from OxCal estimation of event dates, or mean event dates with uncertainties. The code is developed with the aim of producing long-term time-dependent earthquake rate models for individual faults.

If you use this code or the associated data in a publication, please cite this paper: Griffin, J. D., M. W. Stirling and T. Wang (2020). Periodicity and Clustering in the Long‐Term Earthquake Record, Geophys. Res. Lett., 47 (22). https://doi.org/10.1029/2020GL089272

To cite the code directly, use the DOI: https://doi.org/10.5281/zenodo.4001031

The main script for running the analysis is `./analysis/fault_B_M_analysis.py`. This script is used to generate Monte Carlo samples for earthquake chronologies and calculate coefficient of variation (COV), burstiness (B) and memory coefficients (M) and their respective uncertainties from the data.

Data is stored in the directory `./data`. Data is either .txt files containing event times and uncertainties, or `.csv` files containing OxCal outputs. Each data file is associated with a parameter file in `./params` that describes how the event timing data is given and some basic information abou the fault.
