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

Solves the universal equation of kepler. In this equation, G_{2} and G_{3} are the third and fourth g-functions, respectively.

```math
r_{0}X + \\eta_{0}G_{2} + \\zeta_{0}G_{3} - dt = 0
```

# Args

* `x0`: initial estimation of the universal anomaly.
* `beta`: a coefficient calculated in function `KeplerFlow`, which equals to `µ`/`a` (where `a` is the semi-major axis and `mu` the standard gravitational parameters of the bodies involved).
* `eta0`: a coefficient that is needed to solve, calculated in function `KeplerFlow`.
* `zeta0`: another quantity of the equation it is needed to solve 
* `r0`: initial distance between the two bodies.
* `Gz`: array of the first four g-functions. Instead of defining it on each execution of the function, the same array is used in the whole simulation. That is why it's passed as an argument.
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




"""
    keplerFlow(r0_bek, v0_bek, dt, mu, Gz, trace = false)

Defines the position and speeds of a particle -after a time-step- orbiting another body given its initial position and speeds.

# Args

* `r0_bek`: position vector of the secondary body orbiting the primary body.
* `v0_bek`: velocity vector of the secondary body orbiting the primary body.
* `dt`: time-step, the function returns the new position and speed vectors after dt time has passed.
* `mu`: standard gravitational parameter of the two bodies.
* `r0`: initial distance between the two bodies.
* `Gz`: array of the first four g-functions. Instead of defining it on each execution of the function, the same array is used in the whole simulation. That is why it's passed as an argument.
* `trace`: boolean value that shows information of the execution if true.

# Returns

* `r_bek`: the position vector of the body after dt time from the initial values.
* `v_bek`: the position vector of the body after dt time from the initial values.
"""
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




"""
    keplerFlowWithKahan(r0_bek, v0_bek, dt, mu, Gz, trace, r_k, v_k)

Defines the position and speeds of a particle -after a time-step- orbiting another body given its initial position and speeds. When the new vectors are calculated, Kahan's algorithm is applied. 

# Args

* `r0_bek`: position vector of the secondary body orbiting the primary body.
* `v0_bek`: speed vector of the secondary body orbiting the primary body.
* `dt`: time-step, the function returns the new position and speed vectors after dt time has passed.
* `mu`: standard gravitational parameter of the two bodies.
* `r0`: initial distance between the two bodies.
* `Gz`: array of the first four g-functions. Instead of defining it on each execution of the function, the same array is used in the whole simulation. That is why it's passed as an argument.
* `trace`: boolean value that shows information of the execution if true.
* `r_k`: lost information in the sum of the ``f`` and ``g`` coefficients is saved here and used in the next iteration (for positions).
* `v_k`: lost information in the sum of the ``\\dot{f}`` and ``\\dot{g}`` coefficients is saved here and used in the next iteration (for velocities).

# Returns

* `r_bek`: the position vector of the body after dt time-step from the initial values.
* `v_bek`: the velocity vector of the body after dt time-step from the initial values.
"""
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




"""
    keplerStep!(r, v, dt, mu, Gz, trace)

Drifts all particles/bodies under ``\\mathcal{H}_{\\mathrm{Kepler}}`` hamiltonian.

```math
\\mathcal{H}_{\\mathrm{Kepler}}(\\textbf{q}, \\textbf{p}) = \\sum_{i=2}^{N} \\frac{\\textbf{p'}^{2}_{i}}{2 m_{i}'} - \\sum_{i=2}^{N} \\frac{Gm_{i}^{'}M_{i}}{\\mid \\textbf{q}_{i}^{'} \\mid} 
```

# Args

* `r`: NxD matrix, where N is the number of bodies and D the dimensions. Initially must contain the initial positions of the bodies, but the output positions will be saved here as well (NxD, where N is the number of bodies  and D is the number of dimensions).
* `v`: NxD matrix. Initially must contain the initial speeds of the bodies, but the output speeds will be saved here as well (NxD, where N is the number of bodies and D is the number of dimensions).
* `dt`: time-step between the initial and final positions and speeds.
* `mu`: N array. Standard gravitational parameter of the N bodies.
* `r0`: initial distance between the two bodies.
* `Gz`: array of the first four g-functions. Instead of defining it on each execution of the function, the same array is used in the whole simulation. That is why it's passed as an argument.
* `trace`: boolean value that shows information of the execution if true.
"""
function keplerStep!(r, v, dt, mu, Gz, trace)
    # drift all particles under H_Kepler
    for i in 2:size(r)[1]
        r[i,:], v[i,:] = keplerFlow(r[i,:], v[i,:], dt, mu[i], Gz, trace)
    end
end




"""
    keplerStepWithKahan!(r, v, dt, mu, Gz, trace, r_k, v_k)

Drifts all particles/bodies under ``\\mathcal{H}_{\\mathrm{Kepler}}`` hamiltonian.

```math
\\mathcal{H}_{\\mathrm{Kepler}}(\\textbf{q}, \\textbf{p}) = \\sum_{i=2}^{N} \\frac{\\textbf{p'}^{2}_{i}}{2 m_{i}'} - \\sum_{i=2}^{N} \\frac{Gm_{i}^{'}M_{i}}{\\mid \\textbf{q}_{i}^{'} \\mid} 
```

# Args

* `r`: NxD matrix, where N is the number of bodies and D the dimensions. Initially must contain the initial positions of the bodies, but the output positions will be saved here as well.
* `v`: NxD matrix. Initially must contain the initial speeds of the bodies, but the output speeds will be saved here as well.
* `dt`: time-step between the initial and final positions and speeds.
* `mu`: array of N elements, containing standard gravitational parameters of the N bodies.
* `r0`: initial distance between the two bodies.
* `Gz`: array of the first four g-functions. Instead of defining it on each execution of the function, the same array is used in the whole simulation. That is why it's passed as an argument.
* `trace`: boolean value that shows information of the execution if true.
* `r_k`: auxiliary variable to save extra information of Kahan's algorithm (positions). Same dimension as `r`.
* `v_k`: auxiliary variable to save extra information of Kahan's algorithm (velocities). Same dimension as `v`.
"""
function keplerStepWithKahan!(r, v, dt, mu, Gz, trace, r_k, v_k)
    ### drift all particles under H_Kepler  
    for i in 2:size(r)[1]
        r[i,:], v[i,:] =keplerFlowWithKahan(r[i,:], v[i,:], dt, mu[i], Gz, trace, r_k[i,:], v_k[i,:])
    end
end