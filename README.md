# Generation of stochastic bridges for 2D stochastic processes

This repository provides a Fortran implementation for generating stochastic bridges of a two-dimensional stochastic process of the form:

$\dot{x} = F_x(x,y) + \sqrt{D} G_x(x,y) \xi_x(t),$

$\dot{y} = F_y(x,y) + \sqrt{D} G_y(x,y) \xi_y(t),$

where $\xi_x(t)$ and $\xi_y(t)$ are independent Gaussian white noise variables with zero mean and correlations $\langle \xi_i(t)\xi_j(t')\rangle = \delta_{i,j}\delta(t-t')$. The equations are interpreted in the Itô sense.

Stochastic bridges are realizations of the process conditioned to start at $(x(t=0), y(t=0)) = (x_0, y_0)$ and end at $(x(t=T), y(t=T)) = (x_T, y_T)$.

The code `Bridges2D.f90` allows the generation of stochastic bridges using the backtracking method proposed in *Sampling rare trajectories using stochastic bridges*, where the theoretical derivation is presented in detail only for one-dimensional systems. Its extension to two dimensions is outlined in the Appendices of *The nature of stochastic fluctuations shapes transition dynamics in cell-type switching*. The code `dranxor.f90` implements the random number generator used in the simulations.

The implementation is designed for reproducibility and can be readily adapted to a broad class of two-dimensional stochastic processes.

## Repository Structure

The repository contains the following files:

```text
Bridges2D.f90     	  # Stochastic bridge generator     
dranxor.f90			  # Random number generator
PQS_geneticswitch.txt # Quasi-stationary distribution for the genetic togle switch model under additive noise with intensity D = 0.0025.
```

## Compilation

All codes were compiled and tested with Intel Fortran Compiler:

   ifort (IFORT) 2021.10.0 20230609

The following compilation command was used: 

  ifort Bridges2D.f90 dranxor.f90 -O3 -no-prec-div -fp-model fast=2 -march=sandybridge -mtune=core-avx2 -o program_name.x
  
## How to cite

If you use the stochastic bridge generator in your research, please cite the following references:

> *Sampling rare trajectories using stochastic bridges*
> 
> Javier Aguilar, Joseph W. Baron, Tobias Galla, and Raúl Toral
>
> Phys. Rev. E 105, 064138 (2022)
>
> https://doi.org/10.1103/PhysRevE.105.064138

> *The nature of stochastic fluctuations shapes transition dynamics in cell-type switching*
> 
> Sara Oliver-Bonafoux, Javier Aguilar, Tobias Galla, and Raúl Toral



