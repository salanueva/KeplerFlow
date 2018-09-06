# KeplerFlow


DISCLAIMER: This package was made in my final degree project, and its final goal is to share the code that was written during that work. It is probable that there won't be any other updates in this package. Although the project was mainly written in basque, the package is fully documented in english. 


The code is mainly based on the following paper: [WHFast: A fast and unbiased implementation of a symplectic Wisdom-Holman integrator for long term gravitational simulations](https://arxiv.org/abs/1506.01084) by Hanno Rein and Daniel Tamayo. 


## Adding this package

In order to add this Julia package (Version 1.0 or newer), you need to enter in Julia, type "]" and add the package via URL:

```julia
julia> ]
(v1.0) pkg> add https://github.com/salanueva/KeplerFlow.jl
```

Whenever it's wanted to use, don't forget to initialize the package executing the following:

```julia
julia> using KeplerFlow
```

If there is a doubt about how a function works (*calculateEnergy*, for example), check: what it does, what arguments it has and what it returns, typing ? plus the function in the console.

```julia
julia> ?
?help> calculateEnergy
```

The package offers two main algorithms which solve two different problems: the Kepler problem, and the N-body problem. Some minor changes have been made in order to improve execution-times and the accuracy of the solvers. Moreover, these algorithms are prepared to be executed with different floating point arithmetics (*Float64* and *BigFloat*, for example).


## Kepler Problem

With the initial position and velocity vectors of a body, Kepler's problem consists on calculating the new position and velocity vectors after `dt` time-step.

This is achieved solving the universal Kepler's equation, when the value of the *universal anomaly* is known. The new position **r** and velocity **v** vectors can be calculated, using Gauss's `f` and `g` coefficients, which are defined with the *universal anomaly*. 

These `f` and `g` coefficients are really small (values near 0); so, in order to avoid loss of information, a variation of this method that includes Kahan's summation algorithm is defined.

### Usage

Considering that `r0` and `v0` are the initial position and velocity 3-element arrays, and that `mu` equals to the sum of the standard gravitational parameter of the bodies involved in the Kepler Problem, *keplerFlow* function evolves the system under a `dt` time-step. This function returns the new `r` position and `v` velocity vectors after that time-step.  

```julia
julia> Gz = zeros(Float64, 4)
julia> r, v = keplerFlow(r0, v0, dt, mu, Gz, false)
```

`Gz` is an auxiliary array in which Stiefel's G-functions are saved (used to calculate `f` and `g` coefficients), and the last argument of the function specifies if we want to debug what values are calculated during its execution.


## N-body Problem

This problem's main goal is to evolve a planetary system under a time-step `dt` over and over until `t_max` is reached. This is a problem that has to be solved using a numerical approach.

So, in order to evolve N-body systems, a Wisdom-Holfman map is coded. Shortly, it evolves the system in two steps: it calculates keplerian orbits of each body as a base (*Drift* step), and then perturbates those orbits with forces that are created between the bodies (*Kick* step).

### Usage

Let's think we are working with a 4-body system. The following variables will be our input data: `pos` initial positions, `vel` initial velocities, `mass` masses of each body and `g` the gravitational constant which is defined depending on which units are we using (time is measured in days, mass in kilograms and length in Astronomical units). 

```julia
julia> pos
4×3 Array{Float64,2}:
  0.0        0.0         0.0       
 -0.140728  -0.443901   -0.0233456 
 -0.71863   -0.0225038   0.0411718 
 -0.168525   0.968783   -4.12097e-6  
julia> vel
4×3 Array{Float64,2}:
  0.0           0.0          0.0        
  0.0211689    -0.00709798  -0.00252283 
  0.000513533  -0.0203061   -0.000307175
 -0.0172339    -0.00300766   3.56293e-8 
julia> mass
4-element Array{Float64,1}:
 1.988544e30
 3.302e23   
 4.8685e24  
 5.97219e24   
julia> g
1.4881361162906673e-34
```

To determine some aspects of the simulation of this system, the following dictionary is defined, in which 4 different parameters (key/value pair) can be defined.

```julia
sldict = Dict() # Simulation Logic Dictionary
sldict["bool_energy"] = true # the energy of the whole system will be calculated on each output step
sldict["bool_print"] = false # true if we want to debug what values are calculated during the *Drift* step
sldict["ignore_H0"] = true # true if we want to ignore the linear path of the whole system
sldict["out_step"] = 10 # specifies how many iterations are needed to calculate each output
```

With the above data, we can evolve/simulate this 4-body system during 100 years (36524 days). Each output will be calculated after 10 iterations (see `sldict["out_step"]`), and each step/iteration evolves the system 0.1 days (so, outputs will be calculated each day).

```julia
pos_out, vel_out, energy = simulation(pos, vel, mass, sldict, g, 0.1, 36524)
```

If `energy` isn't calculated, there will only be two output variables: `pos_out` and `vel_out`. For more info, run *?simulation*.

### Saving output

In order to show or save outputs of this algorithm (both positions and velocities, and also energies), two different functions are coded. The first one shows the calculated values in a proper order (`names` being an array of strings or the names of the bodies, and `dt` being time-step; both are optional):

```julia
julia> showOutput(pos_out, vel_out, energy, names, dt)
--- Iteration 0 - Time 0.0 ---
Total energy: -2.005494155508561e21
NAME: X, Y, Z; Vx, Vy, Vz
Sun: 0.0, 0.0, 0.0; 0.0, 0.0, 0.0
Mercury: -0.1407280799640445, -0.443900957766333, -0.02334555971312334; 0.02116887137167173, -0.007097975438870807, -0.002522830951443754
Venus: -0.7186302169039649, -0.02250380105571625, 0.04117184137682463; 0.0005135327579269579, -0.02030614162239802, -0.0003071745100210852
Earth: -0.1685246489174995, 0.9687833048228511, -4.120973411130758e-6; -0.01723394583068879, -0.003007660259271771, 3.562931614781975e-8

--- Iteration 1 - Time 1.0 ---
Total energy: -2.005494155508552e21
NAME: X, Y, Z; Vx, Vy, Vz
Sun: 4.6175542887551935e-8, 6.024204220197825e-8, 1.2048192830016274e-9; 4.5364487630784114e-8, 6.055041989424461e-8, 1.2384793209705342e-9
Mercury: -0.11936426213119027, -0.4503482053289031, -0.025833125282562472; 0.021548190261186102, -0.005793942654687121, -0.0024511197135139097
Venus: -0.717832079354024, -0.04279828873880813, 0.04084840178334206; 0.0010825482167922886, -0.020280216222954638, -0.0003396607610764125
Earth: -0.18573142969073006, 0.9656251070531556, -4.083523899056571e-6; -0.017178769950213717, -0.0033086324500976043, 3.809688954591147e-8

...
``` 

The last one, *saveOutput*, works for saving output data in a .txt file. `delimeter` specifies what characters are written between two values (these characters aren't added at the end of the line). `digits` defines how many digits of each value are going to be written in the file, as maximum. 

```julia
julia> file_name = "output.txt"
julia> delimeter = " "
julia> digits = 16
julia> saveOutput(file_name, pos_out, v_out, energy, dt, delimeter, digits)
``` 