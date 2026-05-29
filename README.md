# multi-canyon-pdem
PDEM-based stochastic seismic response analysis of valley groups under SH waves
![Computational framework](Flowchart.svg)
This repository contains the MATLAB code for the paper  
**"[Paper Title]"** submitted to *Computers & Geosciences*.  
It implements the **Probability Density Evolution Method (PDEM)** to quantify the uncertainty of surface ground motion caused by random canyon radii and spacings in a multi‑canyon topography under SH‑wave incidence.

![Computational framework](Fig02.png)
*Fig. 2: Overall analysis framework coupling stochastic simulation with a deterministic scattering solver.*

## Authors
- [Your Name] ([email] / ORCID)  
- [Co‑author names if any]

## Requirements
- **MATLAB** R2019b or later (tested on R2022a)
- **Statistics and Machine Learning Toolbox** (for `sobolset`)

## Installation
1. Clone this repository:
   ```bash
   git clone https://github.com/your-username/valley-group-pdem.git
