### calculateMassMu! ###
# DESCRIPTION: given the gravitational constant and the mass of each body, calculates 
#              sums of masses, standard gravitational parameters of the bodies and 
#              jacobian masses
# INPUT: m, g // OUTPUT: m_j, mu, mu_sum, m_sum
# @param m: array of N elements, each containing the mass of each body
# @param g: gravitational constant
# @param m_j: array of N elements, each containing the jacobi mass (m')
# @param mu: array of N elements, each containing the standard gravitational parameter
#            of the body (µ)
# @param mu_sum: array of N elements, where the standard gravitational parameter of 
#                the bodies from 1 to i is saved in the i-th element (µ_sum)
# @param m_sum: array of N elements, where the element i contains the mass of the bodies 
#               from 1 to i (M)

function calculateMassMu!(m, g, m_j, mu, mu_sum, m_sum)
    
    m_sum[1] = m[1]
    mu_sum[1] = m[1] * g
    mu[1] = m[1] * g
    
    for i in 2:length(m)
        mu[i] = m[i] * g
        m_sum[i] = m_sum[i-1] + m[i]
        mu_sum[i] = m_sum[i] * g
        m_j[i] = m[i] * m_sum[i-1] / m_sum[i]
    end
    m_j[1] = m_sum[length(m)]

end