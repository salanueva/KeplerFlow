#="""
    calculateEnergy(r, r_j, r_dist, p_j, m, m_j, mu, mu_sum, ignore_H0)

Calculates the sum of kynetic and potential energy of the N-body hamiltonian system:

```math
H_{\mathrm{N-body}}(\textbf{q}, \textbf{p}) = \frac{\textbf{p'}^{2}_{1}}{2 m'_{1}} + \sum_{i=2}^{N} \frac{\textbf{p'}^{2}_{i}}{2 m_{i}'} - \sum_{i=1}^{N} \sum_{j=i+1}^{N} \frac{Gm_{i}m_{j}}{\mid \textbf{q}_{i} - \textbf{q}_{j} \mid}
```

# INPUT: r, r_j, r_dist, p_j, m, m_j, mu, mu_sum, ignore_H0

# Args

* `r`: NxD matrix where N is the number of bodies and D is the number of dimensions, each row contains the positions of the bodies in cartesian coordinates.
* `r_j`:  NxD matrix, it contains the positions of the bodies in jacobi coordinates.
* `r_dist`: array of N elements, in which distances from the bodies to the barycenter are saved.
* `p_j`: NxD matrix, it contains the momentum of the bodies in jacobi coordinates.
* `m`: array of N elements, where element number i contains the mass of i-th body.
* `m_j`: array of N elements, each containing the i-th jacobi mass (m'[i])
```math
m'[1] = M[N] \\
m'[i] = m_{i} \frac{M_{i-1}}{M_{i}}
```
* `mu`: array of N elements, where element with index i contains the standard gravitational parameter of the i-th body
```math
mu[i] = G \cdot m[j]
```
* `$mu_{sum}$`: array of N elements, where the standard gravitational parameter of the bodies from 1 to i is saved in the i-th element
```math
mu_{sum}[i] = \sum_{j = 1}^{i} G \cdot m[j]
```
* `ignore_H0`: boolean that determines if Hamiltonian $H_{0}$ will be ignored (when calculating the energy).

# Returns 

The energy of the $H_{N-body}$ system.
"""=#
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