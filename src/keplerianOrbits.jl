"""
    stumpff!(beta, x, Gz)

Calculates the values of the c-functions and returns the first four g-functions, which will be used to solve Kepler's equation.

```math
c_{0} = \\cos(\\sqrt{z}) \\hspace{3cm} c_{1} = \\frac{\\sin(\\sqrt{z})}{\\sqrt{z}} \\hspace{3cm} c_{n} = \\frac{1}{n!} - z c_{n-2} \\hspace{3cm} G_{n} = x^n c_{n}
```

# Args

* `beta`: a coefficient calculated in function `KeplerFlow`, which equals to `µ`/`a` (where `a` is the semi-major axis and `mu` the standard gravitational parameters of the bodies involved).
* `x`: current estimation of the universal anomaly to solve the universal Kepler equation.
* `Gz`: array of four elements, in which the first four g-functions will be saved.
"""
function stumpff!(beta, x, Gz)
    z = beta * (x^2)
    
    # C-funtzioak
    if (z < 0)
        Gz[1] = cosh(sqrt(abs(z)))
        Gz[2] = sinh(sqrt(abs(z)))/sqrt(abs(z))
    else
        Gz[1] = cos(sqrt(z))
        Gz[2] = sin(sqrt(z))/sqrt(z)
    end
    
    for i = 3:4
        Gz[i] = (1.0 - Gz[i-2])/z
    end
    
    # G-funtzioak
    for i = 2:4
        Gz[i] *= x^(i-1)
    end
end




"""
    keplerSolve!(x0, beta, eta0, zeta0, r0, Gz, dt, trace)

Solves the universal equation of kepler. In this equation, G_{2} and G_{3} are the third and fourth

```math
r_{0}X + \\eta_{0}G_{2} + \\zeta_{0}G_{3} - dt = 0
```

# Args

* `x0`: initial estimation of the universal anomaly.
* `beta`: a coefficient calculated in function `KeplerFlow`, which equals to `µ`/`a` (where `a` is the semi-major axis and `mu` the standard gravitational parameters of the bodies involved).
* `eta0`: a coefficient that is needed to solve, calculated in function `KeplerFlow`.
* `zeta0`: another quantity of the equation it is needed to solve 
* `r0`: initial distance between the two bodies.
* `Gz`: array of the first four g-functions. Instead of defining it on each execution of the function, the same array is used in the whole simulation
* `dt`: time-step of the Kepler-flow.
* `trace`: boolean value that shows information of the execution if true.

# Returns

* `r`: distance between the two bodies after time-step dt.
"""
function keplerSolve!(x0, beta, eta0, zeta0, r0, Gz, dt, trace)
    
    xi = x0 # initial estimation
    r = r0
    imax = 50
    error = Inf
    tol = sqrt(eps(typeof(xi))) / 100
    
    for i in 1:imax
        
        stumpff!(beta, xi, Gz)
        fx = r0*xi + eta0*Gz[3] + zeta0*Gz[4] - dt
        r = r0 + eta0*Gz[2] + zeta0*Gz[3]
        eta = eta0*Gz[1] + zeta0*Gz[2]
        eps = fx / (r - (fx*eta)/(2*r))
        xi -= eps
        error = abs(eps)
        
        if trace
            println("### $i. Iteration ###")
            println("GZ: $Gz")
            println("R: $r")
            println("X$i: $xi")
            println("Epsilon: $eps\n")
        end
        
        if error < tol
            if trace
                println("Stop-condition achieved.")
                println("Universal anomaly X = $xi\n")
            end
            return r
        end
        
    end
    
    if trace
        println("Newton's method didn't converge.")
        println("Last calculated universal anomaly X = $xi\n")
    end
    
    return r
end




### KeplerFlow ###
# DESCRIPTION: defines the position and speeds of a particle -after a time-step- orbiting 
#              another body given its initial position and speeds
# INPUT: r0_bek, v0_bek, dt, mu, trace // AUXILIAR: Gz
# @param r0_bek: position vector of the secondary body orbiting the primary body
# @param v0_bek: speed vector of the secondary body orbiting the primary body
# @param dt: time-step, the function returns the new position and speed vectors after dt time has passed
# @param mu: standard gravitational parameter of the two bodies
# @param Gz: array of the 4 g-functions. Instead of defining it on each execution of 
#            the function, the same array is used in the whole simulation
# @param trace: boolean value that shows information of the execution if true
# @return r_bek: the position vector of the body after dt time from the initial values
# @return v_bek: the velocity vector of the body after dt time from the initial values

function keplerFlow(r0_bek, v0_bek, dt, mu, Gz, trace = false)
    
    r0 = norm(r0_bek)
    v02 = dot(v0_bek,v0_bek) 
    eta = dot(r0_bek,v0_bek)
    alpha = mu/r0
    beta = 2.0*alpha - v02
    zeta = mu - beta*r0
     
    # with beta, we can know the shape of the orbit
    if trace
        if beta > 0.0
            println("Orbit is elliptic.\n")
        elseif beta < 0.0
            println("Orbit is hyperbolic.\n")
        else beta == 0.0
            println("Orbit is parabolic.\n")
        end
    end
    
    ### Initial estimation ###
    x0 = (dt/r0) * (1.0 - eta/2.0) 
    if trace
        println("X0 = $x0\n")
    end
    
    ### Solve equation ###
    r = keplerSolve!(x0, beta, eta, zeta, r0, Gz, dt, trace)
        
    rinv = 1.0/r
    ### New position vector ###
    f = -alpha*Gz[3]
    g = r0*Gz[2] + eta*Gz[3]
   
    r_bek = f * r0_bek + g * v0_bek + r0_bek
    
    ### New speed vector ###
    df = -alpha*Gz[2]*rinv
    dg = -mu*Gz[3]*rinv
    
    v_bek = df * r0_bek + dg * v0_bek + v0_bek
    
    if trace
        println("New position vector: R = $r0_bek")
        println("New speed vector: V = $v0_bek")
    end
    
    return r_bek, v_bek
    
end




### KeplerFlowWithKahan ###
# DESCRIPTION: defines the position and speeds of a particle -after a time-step- orbiting 
#              another body given its initial position and speeds
# INPUT: r0_bek, v0_bek, dt, mu, trace, r_k, v_k // AUXILIAR: Gz // OUTPUT: r_k, v_k
# @param r0_bek: position vector of the secondary body orbiting the primary body
# @param v0_bek: speed vector of the secondary body orbiting the primary body
# @param dt: time-step, the function returns the new position and speed vectors after dt time has passed
# @param mu: standard gravitational parameter of the two bodies
# @param Gz: array of the 4 g-functions. Instead of defining it on each execution of 
#            the function, the same array is used in the whole simulation
# @param r_k: lost information in the sum of the f and g coefficients is saved here and 
#             used in the next iteration (for positions)
# @param v_k: lost information in the sum of the f and g coefficients is saved here and 
#             used in the next iteration (for velocities)
# @param trace: boolean value that shows information of the execution if true
# @return r_bek: the position vector of the body after dt time from the initial values
# @return v_bek: the velocity vector of the body after dt time from the initial values
function keplerFlowWithKahan(r0_bek, v0_bek, dt, mu, Gz, trace, r_k, v_k)
        
    r0 = norm(r0_bek)
    v02 = dot(v0_bek,v0_bek) 
    eta = dot(r0_bek,v0_bek)
    alpha = mu/r0
    beta = 2.0*alpha - v02
    zeta = mu - beta*r0
     
    # with beta, we can know the shape of the orbit
    if trace
        if beta > 0.0
            println("Orbit is elliptic.\n")
        elseif beta < 0.0
            println("Orbit is hyperbolic.\n")
        else beta == 0.0
            println("Orbit is parabolic.\n")
        end
    end
    
    ### Initial estimation ###
    x0 = (dt/r0) * (1.0 - eta/2.0) 
    if trace
        println("X0 = $x0\n")
    end
    
    ### Solve equation ###
    r = keplerSolve!(x0, beta, eta, zeta, r0, Gz, dt, trace)
        
    rinv = 1.0/r
    ### New position vector ###
    f = -alpha*Gz[3]
    g = r0*Gz[2] + eta*Gz[3]
   
    r0_bek_aux = copy(r0_bek)
    kahanSum!(r0_bek, f * r0_bek + g * v0_bek, r_k)
    
    ### New speed vector ###
    df = -alpha*Gz[2]*rinv
    dg = -mu*Gz[3]*rinv
    
    kahanSum!(v0_bek, df * r0_bek_aux + dg * v0_bek, v_k)
    
    if trace
        println("New position vector: R = $r0_bek")
        println("New speed vector: V = $v0_bek")
    end
    
    return r0_bek, v0_bek
end




### keplerStep! ###
# DESCRIPTION: drifts all particles/bodies under H_Kepler hamiltonian
# INPUT: r, v, dt, mu, trace // AUXILIAR: Gz // OUTPUT: r, v
# @param r: initially must contain the initial positions of the bodies, but the output
#           positions will be saved here as well (NxD, where N is the number of bodies 
#           and D is the number of dimensions)
# @param v: initially must contain the initial speeds of the bodies, but the output
#           speeds will be saved here as well (NxD, where N is the number of bodies 
#           and D is the number of dimensions)
# @param dt: time-step between the initial and final positions and speeds
# @param mu: standard gravitational parameters of the two bodies
# @param Gz: array of the 4 g-functions. Instead of defining it on each execution of 
#            the function, the same array is used in the whole simulation
# @param trace: boolean value that shows information of the execution if true

function keplerStep!(r, v, dt, mu, Gz, trace)
    # drift all particles under H_Kepler
    for i in 2:size(r)[1]
        r[i,:], v[i,:] = keplerFlow(r[i,:], v[i,:], dt, mu[i], Gz, trace)
    end
end




### keplerStepWithKahan! ###
# DESCRIPTION: drifts all particles/bodies under H_Kepler hamiltonian
# INPUT: r, v, dt, mu, trace, r_k, v_k // AUXILIAR: Gz // OUTPUT: r, v, r_k, v_k
# @param r: initially must contain the initial positions of the bodies, but the output
#           positions will be saved here as well (NxD, where N is the number of bodies 
#           and D is the number of dimensions)
# @param v: initially must contain the initial speeds of the bodies, but the output
#           speeds will be saved here as well (NxD, where N is the number of bodies 
#           and D is the number of dimensions)
# @param dt: time-step between the initial and final positions and speeds
# @param mu: standard gravitational parameters of the two bodies
# @param Gz: array of the 4 g-functions. Instead of defining it on each execution of 
#            the function, the same array is used in the whole simulation
# @param trace: boolean value that shows information of the execution if true
# @param r_k: auxiliary variable to save extra information of Kahan's algorithm (positions)
# @param v_k: auxiliary variable to save extra information of Kahan's algorithm (velocities)

function keplerStepWithKahan!(r, v, dt, mu, Gz, trace, r_k, v_k)
    ### drift all particles under H_Kepler  
    for i in 2:size(r)[1]
        r[i,:], v[i,:] =keplerFlowWithKahan(r[i,:], v[i,:], dt, mu[i], Gz, trace, r_k[i,:], v_k[i,:])
    end
end


