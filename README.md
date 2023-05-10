# PredictingMotifEvo
Repository for the manuscript "Predicting evolutionary outcomes through the probability of accessing sequence variants"

Contains the relevant code used for the simulations and analyses in the manuscript.

## Evolutionary simulations
The simulations can be run stepwise through QuasiEvolveV3-Selection-DualFitnessPeaks.ipynb or from the commandline using QuasiEvolveV3-SelectionCMD.py.

A user defined nucleotide sequence gets replicated in accordance with influenza replication and influenza mutation rates to simulate sequence evolution and the emergence of potential new motifs. Sequences are checked against user defined motif regular expressions and saved in pandas dataframes.

## Numerical evolutionary calculations and plots
A numerical analysis to determine growth and population sizes of idealized populations can be performed using Numerical-Accessibility.ipynb. Here user defined relative fitness and accessibility are compared, and the sizes and average emergence times of new hypothetical populations are plotted.

## Network analysis
From the output from the evolutionary simulations, the Network-population-sampling scripts create output files that can be uploaded into cytoscape to analyse as networks and create network visualisations. It takes source-target nucleotide sequences fomr simulations as input and calculates probability weights and outputs as csv for cytoscape or other network use.

## Accessibility profiles
Visualisation tool to create satellite plots for codon accessibility based on the probabilities calculated based on the probability function in CommonMotifProb.py

## CommonMotifProb and Pyvolve_Sim
Contain general probability functions imported by all other scripts and several examples of probability calculations for codons and short motifs respectively.
