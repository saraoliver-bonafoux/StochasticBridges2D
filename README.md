# Generation of stochastic bridges for 2D continuous stochastic processes

This repository provides a Fortran implementation for generating stochastic bridges of two-dimensional continuous stochastic processes. We consider processes in continuous space and time described by a system of stochastic differential equations of the form:

$\dot{x} = F_x(x,y) + \sqrt{D} G_x(x,y) \xi_x(t)$

$\dot{y} = F_y(x,y) + \sqrt{D} G_y(x,y) \xi_y(t)$

where $\xi_x(t)$ and $\xi_y(t)$ are independent Gaussian white noise variables with zero mean and correlations $\langle \xi_i(t)\xi_j(t')\rangle = \delta_{i,j}\delta(t-t')$. The parameter $D > 0$ sets the noise intensity. The equations are interpreted in the Itô sense.

A stochastic bridge is a realization of the process conditioned to start at $(x(t=0), y(t=0)) = (x_0, y_0)$ and end at $(x(t=T), y(t=T)) = (x_T, y_T)$.

## Method

The code `Bridges2D_Continuous.f90` generates stochastic bridges using the backtracking method introduced in&nbsp;[1]. The theoretical derivation was originally developed for one-dimensional systems, while its extension to two dimensions is outlined in the Appendices of [2].

The random number generator used by the simulations is implemented in `dranxor.f90`.

### Example: genetic toggle switch

The repository includes an example based on the genetic toggle switch model studied in [2]. The current parameter choice is designed to generate differentiation bridges, connecting the pluripotent state at $t = 0$ to a differentiated state at $t = T$, under demographic noise with intensity $D = 0.025$.

The corresponding quasi-stationary distribution (QSD) of the metastable pluripotent state is also provided. 

The implementation is intended to facilitate reproducibility, but the code can also be readily adapted to other two-dimensional stochastic processes. Extension to higher dimensions is also possible.

## Repository Structure

The repository contains the following files:

```text
Bridges2D_Continuous.f90     	  # Stochastic bridge generator for continuous-state processes 
QSD_geneticswitch_Continuous.txt  # QSD of the genetic toggle switch model under demographic noise
dranxor.f90                       # Random number generator
```

The file `PQS_geneticswitch_Continuous.txt` contains the QSD of the metastable pluripotent state, centered around $(x_0, y_0) = (1, 1)$, for demographic noise with intensity $D = 0.025$.

## Compilation

All codes were compiled and tested with Intel Fortran Compiler:

   ifort (IFORT) 2021.10.0 20230609

The following compilation command was used: 

  ifort Bridges2D_Continuous.f90 dranxor.f90 -O3 -no-prec-div -fp-model fast=2 -march=sandybridge -mtune=core-avx2 -o Bridges2D_Continuous.x
  
## How to cite

If you use the stochastic bridge generator in your research, please cite the following references:

> [1] *Sampling rare trajectories using stochastic bridges*
> 
> Javier Aguilar, Joseph W. Baron, Tobias Galla, and Raúl Toral
>
> Phys. Rev. E 105, 064138 (2022)
>
> https://doi.org/10.1103/PhysRevE.105.064138

> [2] *The nature of stochastic fluctuations shapes transition dynamics in cell-type switching*
> 
> Sara Oliver-Bonafoux, Javier Aguilar, Tobias Galla, and Raúl Toral



