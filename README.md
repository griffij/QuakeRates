# QuakeRates
Code for analysis of long-term earthquake occurrence rates

https://doi.org/10.5281/zenodo.4001031

This code is initially designed to analyse earthquake inter-event times from palaeoseismology data. Inputs are expected to be probability distributions of earthquake occurrence times, such as the output from OxCal estimation of event dates. The code is developed with the aim of producing long-term time-dependent earthquake rate models for individual faults.

The main script for running the analysis is ./analysis/sample_COV.py. This script is used to generate Monte-Carlo samples for earthquake chronologies and calculate coefficient of variation, burstiness and memory coefficients and their respective uncertainties from the data.

Data is stored in the directory ./data. Data is either .txt files containing event times and uncertainties, or .csv files containing OxCal outputs. Each data file is associated with a parameter file in ./params that describes how the event timing data is given and some basic information abou the fault.
