### firstInteractionStep! ###
# DESCRIPTION: calculates accelerations that are defined in the first step of the 
#              H_Interaction Hamiltonian.
# INPUT: r_j, mu_sum // OUTPUT: a
# @param r_j: NxD matrix where N is the number of bodies and D is the number of dimensions,
#             it contains the positions of the bodies in jacobi coordinates
# @param mu_sum: array of N elements, where the standard gravitational parameter of 
#                the bodies from 1 to i is saved in the i-th element (µ_sum)
# @param a: NxD matrix in which the accelerations defined by the first step of the 
#           interaction between bodies are saved
function firstInteractionStep!(r_j, mu_sum, a)
    
    for i in 2:size(r_j)[1]
        a[i,:] = (mu_sum[i] * r_j[i,:]) / norm(r_j[i,:])^3
    end
    
end




### secondInteractionStep! ###
# DESCRIPTION: calculates accelerations that are defined in the second step of the 
#              H_Interaction Hamiltonian.
# INPUT: r_j, m, m_sum, m_j, mu, dp, dp_sum // AUXILIAR: dp, dp_sum // OUTPUT: a
# @param r_j: NxD matrix where N is the number of bodies and D is the number of dimensions,
#             it contains the positions of the bodies in jacobi coordinates
# @param m: array of N elements, where the i-th element is the mass of the i-th body
# @param m_sum: array of N elements, where the sum of the masses of the bodies from 1 
#               to i is saved in the i-th element (µ_sum)
# @param m_j: array of N elements, each containing the jacobi mass (m')
# @param mu: array of N elements, each containing the standard gravitational parameter
#            of each body
# @param dp: auxiliary array of D elements that is defined previously in order to decrease
#            execution time and unneccesary allocations
# @param dp_sum: NxD matrix that is defined previously in order to decrease execution 
#                time and unneccesary allocations
# @param a: NxD matrix in which the accelerations defined by the first step of the 
#           interaction between bodies are saved
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