"""
    simulation(r, v, m, sld::Dict{Any, Any}, g::Number = NaN, t_step::Number = NaN, t_max::Number = NaN)

Runs the Wisdom-Holfman map of the given N-body problem that defines the following Hamiltonian system.

```math
\\mathcal{H}(\\textbf{q}, \\textbf{p}) = \\frac{\\textbf{p'}^{2}_{1}}{2 m'_{1}} + \\sum_{i=2}^{N} \\frac{\\textbf{p'}^{2}_{i}}{2 m_{i}'} - \\sum_{i=2}^{N} \\frac{Gm_{i}^{'}M_{i}}{\\mid \\textbf{q}_{i}^{'} \\mid} + \\sum_{i=2}^{N} \\frac{Gm_{i}^{'}M_{i}}{\\mid \\textbf{q}_{i}^{'} \\mid} - \\sum_{i=1}^{N} \\sum_{j=i+1}^{N} \\frac{Gm_{i}m_{j}}{\\mid \\textbf{q}_{i} - \\textbf{q}_{j} \\mid}
```

# Args

* `r`: NxD matrix, where N is the number of bodies and D the number of dimensions, it contains the initial positions of each body.
* `v`: NxD matrix. It contains the initial velocities of each body.
* `m`: array of N elements. It contains the mass of each body.
* `sld`: simulation logic dictionary, a dictionary where all posible options of the simulation can be specified.
  * `bool_energy`: if true, energy of the hamiltonian will be calculated every output step.
  * `bool_print_kep`: information about the kepler-flow will be displayed in the console.
  * `ignore_H0`: boolean that if true, the linear movement of the system will be taken into account.
  * `out_step`: the output will be calculated on every 'out_step'-th step.
* `g`: gravitational constant.
* `t_step`: time that will be simulated in each step.
* `t_max`: specifies when the simulation will be stopped, after how much time.
"""
function simulation(r, v, m, sld::Dict{Any, Any}, g::Number = NaN, t_step::Number = NaN, t_max::Number = NaN)
    
    ### Analyze given input and calculate some variables
    exit = analyzeInput!(r, v, m, g, t_step, t_max, sld)
    if exit
        return
    end
    
    ### Allocate vectors we'll use during the execution of the program
    r_c = copy(r) # cartesian positions
    v_c = copy(v) # cartesian velocities
    
    r_j = zeros(r) # jacobi positions
    v_j = zeros(r) # jacobi velocities
    
    m_j = zeros(m) # jacobi masses
    mu = zeros(m) # standard gravitational parameters
    mu_sum = zeros(m) # standard gravitational parameters of the sum of inner bodies
    m_sum = zeros(m) # sum of masses of inner bodies
    
    Gz = zeros(typeof(r[1,1]),4) # g-functions will be saved here
    
    # if used values haven't got enough precision, define auxiliary varaibles for Kahan
    if eps(typeof(g)) > 1e-30
        r_kahan = zeros(r)
        v_kahan = zeros(r)
    end
    
    # if there are only 2 bodies, don't define the following
    if size(r)[1] > 2
        a_1 = zeros(r)
        a_2 = zeros(r)

        dp = zeros(r[1,:])
        dp_sum = zeros(r)
        
        if eps(typeof(g)) > 1e-30
            kick_kahan = zeros(r)
        end
    end
    
    # auxiliary variables to use when calculating energy
    if sld["bool_energy"]
        p_j = zeros(r)
        r_dist = zeros(m)
    end
    
    # iteration count and step variables
    dt = t_step # full time step
    half_dt = dt/2.0 # half time step
    iter = 0 # current iteration
    max_iter = Int64(round(t_max/dt)) # number of iterations that the simulation will last
    steps_until_out = Int64(sld["out_step"]) # how many steps are needed to get an output
    index = 1 # how many outputs have been saved
    n_out = Int64(ceil(max_iter / steps_until_out)) + 1 # number of total outputs
    
    # output will be saved here
    r_out = zeros(typeof(r[1,1]), n_out, size(r)[1], size(r)[2])
    v_out = zeros(typeof(r[1,1]), n_out, size(r)[1], size(r)[2])
    if sld["bool_energy"]
        e_out = zeros(typeof(r[1,1]), n_out)
    end
    
    # calculate different mass types and standard gravitational parameter (mu)
    calculateMassMu!(m, g, m_j, mu, mu_sum, m_sum)

    # Calculate jacobi coordinates 
    cartesian2jacobi!(r_c, r_j, m, m_sum)
    cartesian2jacobi!(v_c, v_j, m, m_sum)
    
    # In order to avoid calculating H_0 (or the movement of the whole system)
    if sld["ignore_H0"]
        v_j[1,:] -= v_j[1,:]
    end
    
    r_out[1,:,:] = r_c
    v_out[1,:,:] = v_c
    
    if sld["bool_energy"]
        for i in 1:size(r)[1]
            r_dist[i] = norm(r_j[i,:])
            p_j[i,:] = v_j[i,:]*m_j[i]
        end
        e_out[1] = calculateEnergy(r_c, r_j, r_dist, p_j, m, m_j, mu, mu_sum, sld["ignore_H0"])
    end
    
    ### drift the system under H0 for half a timestep dt
    if !sld["ignore_H0"]
        r_j[1,:] += v_j[1,:]*half_dt
    end
    
    ### drift all particles under H Kepler for half a timestep dt
    if eps(typeof(g)) <= 1e-30
        keplerStep!(r_j, v_j, half_dt, mu_sum, Gz, sld["bool_print"]) 
    else 
        keplerStepWithKahan!(r_j, v_j, half_dt, mu_sum, Gz, sld["bool_print"], r_kahan, v_kahan) 
    end
    total_dt = half_dt
    
    while iter < max_iter
        iter += 1
        steps_until_out -= 1
        
        # calculate 1st part of H Interaction in Jacobi coordinates
        if size(r)[1] > 2
            firstInteractionStep!(r_j, mu_sum, a_1)
            
            # update positions in the inertial frame
            jacobi2cartesian!(r_c, r_j, m, m_sum)

            # calculate 2nd part of H Interaction in inertial frame and convert them into Jacobi coordinates 
            secondInteractionStep!(r_c, m, m_sum, m_j, mu, dp, dp_sum, a_2)
            
            # apply kick from Jacobi accelerations to Jacobi velocities
            if eps(typeof(g)) <= 1e-30
                v_j += dt * (a_1 + a_2)
            else 
                kahanSum!(v_j, dt * (a_1 + a_2), kick_kahan)
            end
        end
    
        if iter != max_iter
            if steps_until_out == 0
                index += 1
                
                ### drift the system under H0 for half a timestep dt
                if !sld["ignore_H0"]
                    r_j[1,:] += v_j[1,:]*half_dt
                end

                ### drift all particles under H Kepler for half a timestep (half_dt)
                if eps(typeof(g)) <= 1e-30
                    keplerStep!(r_j, v_j, half_dt, mu_sum, Gz, sld["bool_print"]) 
                else 
                    keplerStepWithKahan!(r_j, v_j, half_dt, mu_sum, Gz, sld["bool_print"], r_kahan, v_kahan) 
                end
                total_dt += half_dt

                jacobi2cartesian!(r_c, r_j, m, m_sum)
                jacobi2cartesian!(v_c, v_j, m, m_sum)
                
                r_out[index,:,:] = r_c
                v_out[index,:,:] = v_c

                if sld["bool_energy"]
                    for i in 1:size(r)[1]
                        r_dist[i] = norm(r_j[i,:])
                        p_j[i,:] = v_j[i,:]*m_j[i]
                    end
                    e_out[index] = calculateEnergy(r_c, r_j, r_dist, p_j, m, m_j, mu, mu_sum, sld["ignore_H0"])
                end
                
                if !sld["ignore_H0"]
                    r_j[1,:] += v_j[1,:]*half_dt
                end

                if eps(typeof(g)) <= 1e-30
                    keplerStep!(r_j, v_j, half_dt, mu_sum, Gz, sld["bool_print"]) 
                else 
                    keplerStepWithKahan!(r_j, v_j, half_dt, mu_sum, Gz, sld["bool_print"], r_kahan, v_kahan) 
                end
                total_dt += half_dt

                steps_until_out = Int64(sld["out_step"])

            else
                ### drift the system under H0 for a full timestep dt
                if !sld["ignore_H0"]
                    r_j[1,:] += v_j[1,:]*dt
                end
                
                ### drift all particles under H Kepler for a full timestep (dt)
                if eps(typeof(g)) <= 1e-30
                    keplerStep!(r_j, v_j, dt, mu_sum, Gz, sld["bool_print"]) 
                else 
                    keplerStepWithKahan!(r_j, v_j, dt, mu_sum, Gz, sld["bool_print"], r_kahan, v_kahan) 
                end
                total_dt += dt
            end
        end
    end
    
    index += 1
    
    ### drift the system under H0 for half a timestep dt
    if !sld["ignore_H0"]
        r_j[1,:] += v_j[1,:]*half_dt
    end

    if eps(typeof(g)) <= 1e-30
        keplerStep!(r_j, v_j, half_dt, mu_sum, Gz, sld["bool_print"]) 
    else 
        keplerStepWithKahan!(r_j, v_j, half_dt, mu_sum, Gz, sld["bool_print"], r_kahan, v_kahan) 
    end
    total_dt += half_dt

    jacobi2cartesian!(r_c, r_j, m, m_sum)
    jacobi2cartesian!(v_c, v_j, m, m_sum)
    
    r_out[index,:,:] = r_c
    v_out[index,:,:] = v_c
    
    if sld["bool_energy"]
        for i in 1:size(r)[1]
            r_dist[i] = norm(r_j[i,:])
            p_j[i,:] = v_j[i,:]*m_j[i]
        end
        e_out[index] = calculateEnergy(r_c, r_j, r_dist, p_j, m, m_j, mu, mu_sum, sld["ignore_H0"])
    end
    
    if sld["bool_energy"]
        (r_out, v_out, e_out)
    else
        (r_out, v_out)
    end
    
end