"""
### function cartesian2jacobi!(r, r_j, m, m_sum) 

Transforms cartesian coordinates into jacobi coordinates.

# INPUT: r, m, m_sum // OUTPUT: r_j

# Args

* `r`: NxD matrix where N is the number of bodies and D is the number of dimensions, each row contains the positions of the bodies in cartesian coordinates.
* `r_j`: NxD matrix in which the results will be saved.
* `m`: array of N elements, where element number i contains the mass of i-th body.
* `m_sum`: array of N elements, where element with index i contains the mass of the bodies from 1 to i
```math
m_{sum}[i] = \sum_{j = 1}^{i} m[j]
```

"""
function cartesian2jacobi!(r, r_j, m, m_sum) 
    R = m[1]*r[1,:]
    for i in 2:size(r)[1]
        r_j[i,:] = r[i,:] - (R/m_sum[i-1])
        R = R*(1.0+m[i]/m_sum[i-1]) + m[i]*r_j[i,:] 
    end
    r_j[1,:] = R/m_sum[size(r)[1]]
end


"""
### function jacobi2cartesian!(r, r_j, m, m_sum) 

Transforms jacobi coordinates into cartesian coordinates.

# INPUT: r_j, m, m_sum // OUTPUT: r

# Args

* `r`: NxD matrix in which the results will be saved.
* `r_j`: NxD matrix where N is the number of bodies and D is the number of dimensions, each row contains the positions of the bodies in jacobi coordinates.
* `m`: array of N elements, where element number i contains the mass of i-th body.
* `m_sum`: array of N elements, where element with index i contains the mass of the bodies from 1 to i
```math
m_{sum}[i] = \sum_{j = 1}^{i} m[j]
```

"""
function jacobi2cartesian!(r, r_j, m, m_sum) 
    R = r_j[1,:]*m_sum[size(r)[1]]
    for i in size(r)[1]:-1:2
        R = (R - m[i]*r_j[i,:])/m_sum[i]
        r[i,:] = r_j[i,:] + R
        R *= m_sum[i-1]
    end
    r[1,:] = R/m[1]
end