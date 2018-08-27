### calculateEnergy ###
# DESCRIPTION: calculates the sum of kynetic and potential energy of the hamiltonian 
#              system
# INPUT: r, r_j, r_dist, p_j, m, m_j, mu, mu_sum
# @param r: NxD matrix where N is the number of bodies and D the number of dimensions,
#           it contains the positions of the bodies in cartesian coordinates
# @param r_j: NxD matrix where N is the number of bodies and D the number of dimensions,
#             it contains the positions of the bodies in jacobi coordinates
# @param r_dist: array of N elements, in which distances from the bodies to the 
#                barycenter are saved
# @param p_j: NxD matrix where N is the number of bodies and D the number of dimensions,
#             it contains the momentum of the bodies in jacobi coordinates
# @param m: array of N elements, each containing the mass of each body
# @param m_j: array of N elements, each containing the jacobi mass (m')
# @param mu: array of N elements, each containing the standard gravitational parameter
#            of the body (µ)
# @param mu_sum: array of N elements, where the standard gravitational parameter of 
#                the bodies from 1 to i is saved in the i-th element (µ_sum)
# @param ignore_H0: a boolean that says if H0 will be ignored or not
# @return energy: returns the energy of the hamiltonian system

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