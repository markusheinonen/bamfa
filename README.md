# BAMFA

Bayesian metabolic flux analysis MATLAB package by Aalto University. Constructs a Bayesian metabolic model and samples the flux posterior, that is flux vectors that are compatible with steady-state, flux upper and lower bounds and flux measurements. The flux posterior is represented by a set of flux samples, or flux pair samples.

Drop-in replacement to standard flux balance analysis (FBA), metabolic flux analysis (MFA) and flux variability analysis (FVA).

Requirements:
- MATLAB
- OpenCobra toolbox installed

### How to run:

Run `demo.m` in MATLAB or run

```
% run FBA on the ecoli core metabolism (77 reactions)

addpath code;
initCobraToolbox;

% read model
model = readCbModel('models/Ec_core_flux1','fileType','SBML');

% run baseline FBA
fba = optimizeCbModel(model);

% sample, uses by default 10 chains of 200 samples with thinning = 50
sol = bfba(model, fba);

%% visualise

% plot first 15 flux distributions with fba
figure; plotfluxes(model, sol, 1:15, fba);

% plot 7x7 flux pair grid
figure; plotfluxpair(model, sol, [10 2 40 42 45 53 55], fba);

% plot 8 example flux pairs
figure; plotfluxpair2(model, sol, [30 40; 30 41; 40 47; 15 17; 4 5; 1 6; 10 18; 42 45], fba);
```


### Flux distributions

Below is visualised 9 flux pair distributions in three conditions (red/blue/green). The black dot is the classic FBA solution, while the scatter plots indicate all possible flux states.

<p align="center">
  <img src="figures/core_9.png" width="650"/>
</p>
