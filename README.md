#Field-driven Ion Pairing Dynamics in Concentrated Electrolytes by Seokjin Moon and David T. Limmer
Analysis code accompanying the manuscript submitted to arXiv (2026).

---

## Overview
This repository contains the trajectory analysis framework used to compute:
- Cluster-based nearest-counterion distances
- Transition Path Theory (TPT) observables
- Nonequilibrium rate constants
- Conductivity
- Reactive trajectory mechanism analysis

and plotting tools for figures.

All analysis were performed on molecular dynamics trajectories generated using LAMMPS dcd format for concentrated electrolytes under applied electric fields.


---

## Repository Structure
.
├── simulations/ # LAMMPS input files
├── src/ # C++ source files
├── includes/ # C++ header files
├── library/ # C++ source files for dcd readers
├── scripts/ # Shell scripts
├── data/
│ ├── dumps/ # Input LAMMPS DCD files 
│ ├── binary/ # Preprocessed binary trajectories
│ ├── cnnDist/ # Cluster-based nearest-counterion distances
│ ├── cnnAngle/ # Cluster-based nearest-counterion angle
│ ├── cnnId/ # Cluster-based nearest-counterion id of central ion
│ └── cnnCid/ # Cluster-based nearest-counterion id of counter ion
├── results # directory for outputs
├── processedData # all processed data computed from raw trajectories
├── plots # Python plot code using processedData
├── Makefile
├── REPRODUCIBILITY.md
└── README.md

---

## Computational Workflow

The analysis pipeline consists of five stages:

### 1. Preprocessing

Convert LAMMPS '.dcd' trajectories into compact binary format: dcd -> ./data/traj/LiPF6inACN/traj.binary 
The dcd files only saves the position of Li+'s firstly, and subsequently P of PF6-.

Example: 
cd ./script
sh analysis.sh
select keywords, 1) preprocess 2) LiPF6inACN 3) field

---

## 2. Cluster-based Trajectory Analysis

Extract cluster-based nearest-counterion (CBNC) distances from binary trajectories.

Outputs:
- './data/cnnDist/'    > CBNC distance
- './data/cnnAngle/'   > CNBC vector angle to +z direction
- './data/cnnId/'      > the central atom index
- './data/cnnCid/'     > the counterion atom index

Example:
cd ./script
sh analysis.sh
select keywords, 1) processTraj 2) LiPF6inACN 3) field

---

## 3. TPT and Rate Analysis

From 'cnnDist' data, compute:

- Mean committors and populations of states
- Mean first-passage times
- Reactive flux
- Rate constants

Outputs:
- './results/rate/LiPF6inACND05E00.dat' 

Example:
cd ./script
sh analysis.sh
select keywords, 1) rate 2) LiPF6inACN 3) field

---

## 4. Conductivity Calculation

- zero-field conductivity via Einstein-Helfand relation
- Finite-field conductivity via finite-difference method

Outputs:
- './results/conductivity/LiPF6inACND05E00.dat'

Example:
cd ./script
sh analysis.sh
select keywords, 1) conductivity 2) LiPF6inACN 3) field

---

## 5. CBNC Probability Distribution 

From 'cnnDist' data, compute:

- probability distribution function of the CBNC distance

Outputs:
- './results/ionionprob/LiPF6inACND05E00.dat'

Example:
cd ./script
sh analysis.sh
select keywords, 1) ionionprob 2) LiPF6inACN 3) field



## 6. TPT Mechanism Analysis

From 'cnnDist', 'cnnAngle', 'cnnId', and 'cnnCid', extract reactive trajectory traces and compute:

- Hitting probability distributions at inner and outer boundaries

Outputs:
- './results/mechanism/LiPF6inACND05E00AB.r'   >  each row = 1 reactive trajectory, save distance
- './results/mechanism/LiPF6inACND05E00AB.a'   >  each row = 1 reactive trajectory, save angle
- './results/mechanism/LiPF6inACND05E00BA.r'   >  for reverse reaction, B to A
- './results/mechanism/LiPF6inACND05E00BA.a'

## Compilation

- GCC > 9
- C++17 standard

Compile using:
make


## Plots

need Seaborn


