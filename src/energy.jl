"""
    calculateEnergy(r, r_j, r_dist, p_j, m, m_j, mu, mu_sum, ignore_H0)

Calculates the sum of kynetic and potential energy of the N-body hamiltonian system:

```math
\\mathcal{H}_{\\mathrm{N-body}}(\\textbf{q}, \\textbf{p}) = \\sum_{i=1}^{N} \\frac{\\textbf{p'}^{2}_{i}}{2 m_{i}'} - \\sum_{i=1}^{N} \\sum_{j=i+1}^{N} \\frac{Gm_{i}m_{j}}{\\mid \\textbf{q}_{i} - \\textbf{q}_{j} \\mid} \\\\
```

# Args

* `r`: NxD matrix where N is the number of bodies and D is the number of dimensions, each row contains the positions of the bodies in cartesian coordinates.
* `r_j`:  NxD matrix, it contains the positions of the bodies in jacobi coordinates.
* `r_dist`: array of N elements, in which distances from the bodies to the barycenter are saved.
* `p_j`: NxD matrix, it contains the momentum of the bodies in jacobi coordinates.
* `m`: array of N elements, where element number i contains the mass of i-th body.
* `m_j`: array of N elements, each containing the i-th jacobi mass (``m'_{i}``)
```math
m'_{1} = m_{sum_{N}} \\\\
m'_{i} = m_{i} \\frac{m_{sum_{i-1}}}{m_{sum_{i}}} \\\\
\\mathrm{,where } \\hspace{0.5cm} m_{sum_{i}} = \\sum_{j = 1}^{i} m_{j}
```
* `mu`: array of N elements, where element with index i contains the standard gravitational parameter of the i-th body
```math
mu_{i} = G \\cdot m_{j} \\\\
```
* `mu_{sum}`: array of N elements, where the standard gravitational parameter of the bodies from 1 to i is saved in the i-th element
```math
mu_{sum_{i}} = \\sum_{j = 1}^{i} G \\cdot m_{j} \\\\
```
* `ignore_H0`: boolean that determines if Hamiltonian ``\\mathcal{H_{0}}`` will be ignored (when calculating the energy).

# Returns 

The energy of the ``\\mathcal{H}_{N-body}`` system.
"""
function calculateEnergy(r, r_j, r_dist, p_j, m, m_j, mu, mu_sum, ignore_H0)
    energy = 0.0
    
    # H_0
    if !ignore_H0
        energy += (dot(p_j[1,:],p_j[1,:]) / (2.0*m_j[1]))
    end
    
    # H_Kepler
    for i in 2:size(r)[1]
        energy += (dot(p_j[i,:],p_j[i,:]) / (2.0*m_j[i])) #- (mu_sum[i] * m_j[i] / r_dist[i])
    end
    #println(energy)
    
    # H_Interaction step 1
    #for i in 2:size(r)[1]
    #    energy += (mu_sum[i] * m_j[i] / r_dist[i])
    #end
    #println(energy)
    
    # H_Interaction step 2
    for i in 1:size(r)[1]
        for j in (i+1):size(r)[1]
            energy -= (mu[i] * m[j] / norm(r[i,:] - r[j,:]))
        end
    end
    
    return energy
end