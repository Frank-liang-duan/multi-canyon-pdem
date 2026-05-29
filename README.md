# Stochastic Seismic Response Analysis of Valley Groups Based on Point Estimate and Kernel Density Estimation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains MATLAB code for the paper  

The code performs stochastic analysis of surface ground motion for a group of semi‑circular canyons under SH‑wave incidence, accounting for the uncertainty in canyon geometry—namely, the random variability of individual canyon radii and of the spacings between adjacent canyons. By coupling **number‑theoretic point selection**, a **deterministic wave‑function expansion solver**, and **adaptive kernel density estimation**, it efficiently yields the probability density function (PDF), cumulative distribution function (CDF), and statistical descriptors (mean, 95th percentile, coefficient of variation, failure probability) of the surface displacement amplitude.

![Computational framework](Flowchart.svg)  
*Fig. 2: Overall analysis framework coupling stochastic simulation with a deterministic scattering solver.*

## Authors
- [Xueliang Duan] ([email] / [ORCID])  
- [Co‑author names if any]

## Requirements
- MATLAB R2019b or later (tested on R2022a)
- Statistics and Machine Learning Toolbox (required for `sobolset`)

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/valley-group-pdem.git
