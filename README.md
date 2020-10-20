# Quantum Annealing Repository

## What is this?
This repository contains quantum anealling implementations for D-Wave Ocean API maintained by CJS Inc.
They can be run in D-Wave Leap' IDE or locally, according to your preference. 

## Available examples
### Nurse Schedule Problem: Replication of Ikeda *et al.* DOI: 10.1038/s41598-019-49172-3
Reproduction report: https://scigen.report/reviews/get?doi=10.1038/s41598-019-49172-3#aOZZCJO
 - `Nurse Shift.py` performs forward (standard) annealing to find solutions and saves on pickle files.
 - `Reverse Nurse Shift.py` uses solutions and embeddings from `Nurse Shift.py` to perform reverse annealing, i.e., starting from a given initial state increase quantum fluctuations and anneal the state again locally around the given initial state. 
 - `results_analysis.py` prints the first sample of a result (ideally, it will be a ground state) and computes the fraction of the solutions that were in the ground state. It also computes the mean Hamming distance between consecutive samples and their standard deviation. 
#### Usage
 1. Using D-Wave IDE, which can be done by prefixing the URL to this repository with `https://ide.dwavesys.io/#` (requires D-Wave Leap account), first run `Nurse Shift.py` and then `Reverse Nurse Shift.py`. This will produce results for forwar annealing first, then reverse annealling. Before runnning, check the parameters are what you want (e.g., `topology` and `numSamples`).
 2. Adjust the variable `filename` in `results_analysis.py` according to the output you want to see the results.


### Notice
We do **not** have any direct association with D-Wave and their products. Questions about D-Wave services and products should be directed directly to them.
