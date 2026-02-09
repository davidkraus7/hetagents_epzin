# Consumption and Portfolio Choice with Recursive Preferences and Uninsurable Labour Income Risk

Final year undergraduate thesis under the supervision of a faculty advisor. Investigates optimal consumption, savings, and portfolio choice in a continuous-time heterogeneous agents model with Epstein-Zin preferences and uninsurable labour income risk.

## Overview

Standard heterogeneous agents models typically assume power utility, in which the coefficient of relative risk aversion (RRA) is the reciprocal of the intertemporal elasticity of substitution (IES). This paper incorporates recursive preferences (Duffie & Epstein, 1992) to decouple these parameters, and investigates their distinct effects on household decisions and the wealth distribution. The model is solved numerically using a finite-difference scheme for the HJB equation and the adjoint relationship to obtain the stationary wealth distribution from the Kolmogorov-Forward equation.

## Repository structure

```
├── code/
│   ├── twoassetEZ.m                    # Baseline two-asset model with EZ preferences
│   ├── twoassetEZ_ep.m                    # Partial equilibrium extension with time-varying equity premium
│   └── plotting.m                       # Generates all figures
└── docs/
    ├── boundary_conditions.pdf             # Derivation of boundary conditions
    ├── boundary_conditions.tex
    ├── slides.pdf                       # Presentation slides
    └── slides.tex
```

## Notes on boundary conditions

The document `docs/boundary_conditions.pdf` derives the boundary conditions needed to solve the HJB equation. The first boundary condition simply enforces the borrowing constraint. The second boundary condition is imposed using an asymptotic result. Proposition 1 shows that the risky asset share with recursive preferences and uninsurable labour income risk converges to the Merton share as wealth tends to infinity, extending the result in Achdou et al. (2021) to recursive preferences. The proof proceeds via two auxiliary lemmas: Lemma 1 solves the household problem analytically without income risk, and Lemma 2 establishes a homogeneity property of the value function. This boundary condition is critical for the numerical solution.

## Numerical method

1. Discretise the state space
2. Approximate derivatives using an upwind scheme and central differences
3. Impose boundary conditions
4. Iteratively update the value function by solving a sparse linear system
5. Obtain the stationary wealth distribution by solving eigenvalue problem of the transposed generator

## Replication

Requires MATLAB (no additional toolboxes). Run `code/plotting.m` to generate all figures in the paper. The script calls `twoassetEZ(gamma, psi, lambda)` and `twoassetEZ_ep(mu, lambda_mu)` with the appropriate parameters for each experiment. Note that the time-varying equity premium extension (`twoassetEZ_ep.m`) is a partial equilibrium analysis of household decisions. In particular, asset prices are exogenous and do not reflect aggregate risk.