### cartesian2jacobi! ###
# DESCRIPTION: transforms cartesian coordinates into jacobi coordinates
# INPUT: r, m, m_sum // OUTPUT: r_j
# @param r: NxD matrix where N is the number of bodies and D is the number of dimensions,
#           it contains the positions of the bodies in cartesian coordinates
# @param r_j: NxD matrix in which the results will be saved
# @param m: array of N elements, each containing the mass of each body
# @param m_sum: array of N elements, where the element i contains the mass of the bodies 
#               from 1 to i

function cartesian2jacobi!(r, r_j, m, m_sum) 
    R = m[1]*r[1,:]
    for i in 2:size(r)[1]
        r_j[i,:] = r[i,:] - (R/m_sum[i-1])
        R = R*(1.0+m[i]/m_sum[i-1]) + m[i]*r_j[i,:] 
    end
    r_j[1,:] = R/m_sum[size(r)[1]]
end


### jacobi2cartesian! ###
# DESCRIPTION: transforms jacobi coordinates into cartesian coordinates
# INPUT: r_j, m, m_sum // OUTPUT: r
# @param r: NxD matrix in which the results will be saved
# @param r_j: NxD matrix where N is the number of bodies and D is the number of dimensions,
#             it contains the positions of the bodies in jacobi coordinates
# @param m: array of N elements, each containing the mass of each body
# @param m_sum: array of N elements, where the element i contains the mass of the bodies 
#               from 1 to i
function jacobi2cartesian!(r, r_j, m, m_sum) 
    R = r_j[1,:]*m_sum[size(r)[1]]
    for i in size(r)[1]:-1:2
        R = (R - m[i]*r_j[i,:])/m_sum[i]
        r[i,:] = r_j[i,:] + R
        R *= m_sum[i-1]
    end
    r[1,:] = R/m[1]
end