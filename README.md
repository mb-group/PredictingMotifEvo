# PredictingMotifEvo
Repository for the manuscript "Predicting evolutionary outcomes through the probability of accessing sequence variants"

Contains the relevant code used for the simulations and analyses in the manuscript.

## Evolutionary simulations
The simulations can be run from the commandline using QuasiEvolveV3-SelectionCMD.py with -o outputfile.csv as argument.

A user defined nucleotide sequence gets replicated in accordance with influenza replication and influenza mutation rates to simulate sequence evolution and the emergence of potential new motifs. Sequences are checked against user defined motif regular expressions and saved in pandas dataframes.

The output is a csv containing columns for Founder, Gen(eration), Population size,Motif1 (count), Motif2(count), and genotype of the selected founders for the subsequent generation.

## Numerical evolutionary calculations and plots
A numerical analysis to determine growth and population sizes of idealized populations can be performed using Numerical-Accessibility.ipynb. Here user defined relative fitness and accessibility are compared, and the sizes and average emergence times of new hypothetical populations are plotted. Default values are already in the notebook, and can be changed to explore alternative parameters.

## Network analysis
From the output from the evolutionary simulations, the Network-population-sampling script creates output files that can be uploaded into cytoscape to analyse as networks and create network visualisations. It takes source-target nucleotide sequences from simulations as input and calculates probability weights and outputs as csv for cytoscape or other network use. Currently the script generates outputs for cytoscape for the data for supplementary figure 4 based on the data in Data/ .

## Accessibility profiles
Visualisation tool to create satellite plots for codon accessibility based on the probabilities calculated based on the probability function in CommonMotifProb.py

## Historical codon usage analysis
The scipt Analyse_Codons_#3.1 can be run to create the plot from Supplementary data, and the numerical analysis of which codons are used at conserved S/T residues in phosphorylation sites. It uses supplied data from Residue_Dataframes which in turn are based on multiple sequence alignments. The statistical analysis can subsequently be performed using Analyse_ST_codonUse-Statsphospho.ipynb.

## Numerical codon and motif probabilitiy calculations
The script probability_calc.ipynb can be used to explore the calculations used to determine evolution and mutation probabilities between codons or setw of codons that make up motifs, using different mutation probability models.
