
<a id='KeplerFlow-1'></a>

# KeplerFlow


DISCLAIMER: This package was made in my final degree project, and its final goal is to share the code that was written during that work. It is probable that there won't be any other updates in this package. Although the project was mainly written in basque, the package is fully documented in english. 


The code is mainly based on the following paper: [WHFast: A fast and unbiased implementation of a symplectic Wisdom-Holman integrator for long term gravitational simulations](https://arxiv.org/abs/1506.01084) by Hanno Rein and Daniel Tamayo. 

## Adding this package

In order to add this Julia package (Version 1.0 or newer), you need to enter in Julia, type "]" and add the package via URL:

```julia
julia> ]
(v1.0) pkg> add https://github.com/salanueva/KeplerFlow.jl
```

The package offers two main algorithms which solve two different problems: the Kepler problem, and the N-body problem. Some minor changes have been made in order to improve execution-times and the accuracy of the solvers.


<a id='Kepler-Problem-1'></a>

## Kepler Problem


With the initial position and velocity vectors of a body, Kepler's problem consists on calculating the new position and velocity vectors after `dt` time-step.


This is achieved solving the universal Kepler's equation, when the value of the *universal anomaly* is known. The new position **r** and velocity **v** vectors can be calculated, using Gauss's `f` and `g` coefficients, which are defined with the *universal anomaly*. 


These `f` and `g` coefficients are really small (values near 0); so, in order to avoid loss of information, a variation of this method that includes Kahan's summation algorithm is defined.

### Functions

```@docs
KeplerFlow.keplerFlow
KeplerFlow.keplerFlowWithKahan
```

### Usage


<a id='N-body-Problem-1'></a>

## N-body Problem

This problem's main goal is to evolve a planetary system under a time-step `dt`. This is a problem that has to be solved using a numerical approach. So, a Wisdom-Holfman map, defined in the above mentioned article, has been coded.

```@docs
KeplerFlow.simulation
```

### Usage