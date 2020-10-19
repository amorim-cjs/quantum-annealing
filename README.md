# Quantum Annealing Repository

## What is this?
This repository contains quantum anealling implementations for D-Wave Ocean API maintained by CJS Inc.
They can be run in D-Wave Leap' IDE or locally, according to your preference. 

## Available examples
Nurse Schedule Problem: Replication of Ikeda *et al.* DOI: 10.1038/s41598-019-49172-3
 - `Nurse Shift.py` performs forward (standard) annealing to find solutions and saves on pickle files.
 - `Reverse Nurse Shift.py` uses solutions and embeddings from `Nurse Shift.py` to perform reverse annealing, i.e., starting from a given initial state increase quantum fluctuations and anneal the state again locally around the given initial state. 
 - Using D-Wave IDE, which can be done by prefixing the URL to this repository with `https://ide.dwavesys.io/#` (requires D-Wave Leap account), first run `Nurse Shift.py` and then `Reverse Nurse Shift.py`.

### Notice
We do **not** have any direct association with D-Wave and their products. Questions about D-Wave services and products should be directed directly to them.
