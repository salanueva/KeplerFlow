"""
    calculateMassMu!(m, g, m_j, mu, mu_sum, m_sum)

Given the gravitational constant and the mass of each body, calculates sums of masses, standard gravitational parameters of the bodies and jacobian masses

```math
mu_{sum_{i}} = \\sum_{j = 1}^{i} m_{i} \\hspace{0.5cm} mu = G \\cdot m_{i} \\hspace{0.5cm} mu_{sum} = G \\cdot m_{sum_{i}} \\hspace{0.5cm} m' = m_{i} \\frac{m_{sum_{i-1}}}{m_{sum_{i}}} 
```

# Args

* `m`: array of N elements, each containing the mass of each body.
* `g`: gravitational constant.
* `m_j`: array of N elements. Each element will contain the jacobi mass ``m'`` of a body, as output.
* `mu`: array of N elements. Each element will contain the standard gravitational parameter of a body, as output.
* `mu_sum`: array of N elements, where the standard gravitational parameter of the bodies from 1 to i will be saved saved in the i-th element, as output.
* `m_sum`: array of N elements, where the element i contains the mass of the bodies from 1 to i, as output.
"""
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