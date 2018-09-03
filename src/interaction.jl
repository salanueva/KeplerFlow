"""
    firstInteractionStep!(r_j, mu_sum, a)

Calculates accelerations that are defined in the first step of the ``H_{\\mathrm{Interaction}}`` Hamiltonian.

```math
\\mathcal{H}_{\\mathrm{Interaction}_{1}}(\\textbf{q}, \\textbf{p}) = \\sum_{i=1}^{N} \\frac{\\sum_{i=2}^{N} \\frac{Gm'_{i}m_{sum_{i}}}{\\mid \\textbf{r'}_{i}\\mid} \\\\
```

# Args

* `r_j`: NxD matrix where N is the number of bodies and D is the number of dimensions, it contains the positions of the bodies in jacobi coordinates
* `mu_sum`: array of N elements, where the standard gravitational parameter of the bodies from 1 to i is saved in the i-th element (µ_sum)
* `a`: NxD matrix in which the accelerations defined by the first step of the interaction between bodies are saved.
"""
function firstInteractionStep!(r_j, mu_sum, a)
    
    for i in 2:size(r_j)[1]
        a[i,:] = (mu_sum[i] * r_j[i,:]) / norm(r_j[i,:])^3
    end
    
end




"""
    secondInteractionStep!(r, m, m_sum, m_j, mu, dp, dp_sum, a)

Calculates accelerations that are defined in the second step of the ``\\mathcal{H}_{Interaction_{2}}`` Hamiltonian.

```math
\\mathcal{H}_{Interaction_{2}}(\\textbf{q}, \\textbf{p}) = - \\sum_{i=1}^{N} \\sum_{j=i+1}^{N} \\frac{Gm_{i}m_{j}}{\\mid \\textbf{q}_{i} - \\textbf{q}_{j} \\mid} \\\\
```

# Args

* `r_j`: NxD matrix where N is the number of bodies and D is the number of dimensions, it contains the positions of the bodies in jacobi coordinates
* `m`: array of N elements, where the sum of the masses of the bodies from 1 to i is saved in the i-th element.
* `m_sum`: array of N elements, where the sum of the masses of the bodies from 1 to i is saved in the i-th element.
```math
m_{sum_{i}} = \\sum_{j = 1}^{i} m_{j} \\\\
```
* `m_j`: array of N elements, each containing the jacobi mass (m').
```math
m'_{1} = m_{sum_{N}} \\\\
m'_{i} = m_{i} \\frac{m_{sum_{i-1}}}{m_{sum_{i}}} \\\\
\\mathrm{,where } \\hspace{0.5cm} m_{sum_{i}} = \\sum_{j = 1}^{i} m_{j}
```
* `mu`: array of N elements, each containing the standard gravitational parameter of each body.
```math
mu_{i} = G \\cdot m_{j} \\\\
```
* `dp`: auxiliary array of D elements that is defined previously in order to decrease execution time and unneccesary allocations.
* `dp_sum`: NxD matrix that is defined previously in order to decrease execution time and unneccesary allocations.
* `a`: NxD matrix in which the accelerations defined by the first step of the interaction between bodies are saved.
"""
function secondInteractionStep!(r, m, m_sum, m_j, mu, dp, dp_sum, a)
    
    for i in 1:size(r)[1]
        for j in 1:length(dp)
            dp[j] = 0.0
        end
        for j in 1:size(r)[1]
            if i != j
                r_coor = r[i,:] - r[j,:]
                dp -= (mu[j]*m[i]*r_coor / (norm(r_coor)^3)) # dp[i] / dt 
            end
        end
        
        if i == 1
            dp_sum[1,:] = dp
        else
            dp_sum[i,:] = dp + dp_sum[i-1,:]
            a[i,:] = ((m_sum[i-1]*dp - m[i]*dp_sum[i-1,:]) / m_sum[i]) / m_j[i] # azelerazioak jacobira pasa
        end
    end
    
    a[1,:] = dp_sum[size(r)[1],:] / m_j[1]
    
end