"""
    kahanSumOneValue(a, b, c)

`a`+`b`=`sum` is done, where the lost information in that `sum` is saved in `c` (according to Kahan Summation Algorithm). Each input only represents a number.


# Args

* `a`: it is a number, where `a`>`b` should be true for the appropiate work of the algorithm.
* `b`: the value that will be added to `a`.
* `c`: the first time this function is called in an iteration, the value of `c` should be 0. After that, this value should contain the lost information of the last sum, for its proper function. 

# Returns

* `sum`: `a` + `b`.
* `c`: information lost in the sum `a`+`b`.
"""
function kahanSumOneValue(a, b, c)
    sum = a
    y = b - c          # So far, so good: c is zero.
    t = sum + y        # Alas, sum is big, y small, so low-order digits of y are lost.
    c = (t - sum) - y  # (t - sum) cancels the high-order part of y; subtracting y recovers negative (low part of y)
    sum = t            # Algebraically, c should always be zero. Beware overly-aggressive optimizing compilers!
    return sum, c 
end



"""
    kahanSumOneValue(a, b, c)

`a`+`b`=`sum` is done, where the lost information in that `sum` is saved in `c` (according to Kahan Summation Algorithm). Each input can represent an array or a matrix.


# Args

* `a`: it is a number, where `a`>`b` should be true for the appropiate work of the algorithm. `sum` values are overwritten here.
* `b`: the value that will be added to `a`.
* `c`: the first time this function is called in an iteration, the value of `c` should be 0. After that, this value should contain the lost information of the last sum, for its proper function. `c` values are overwritten here.
"""
function kahanSum!(a, b, c)
    
    if length(size(a)) == 1 # Vector
        for i in 1:length(a)
            y = b[i] - c[i]        # So far, so good: c is zero.
            t = a[i] + y           # Alas, a is big, y small, so low-order digits of y are lost.
            c[i] = (t - a[i]) - y  # (t - a) cancels the high-order part of y; subtracting y recovers negative (low part of y)
            a[i] = t               # Algebraically, c should always be zero. Beware overly-aggressive optimizing compilers!
        end
    elseif length(size(a)) == 2 # Matrix
        for i in 1:size(a)[1]
            for j in 1:size(a)[2]
                y = b[i,j] - c[i,j]        # So far, so good: c is zero.
                t = a[i,j] + y           # Alas, a is big, y small, so low-order digits of y are lost.
                c[i,j] = (t - a[i,j]) - y  # (t - a) cancels the high-order part of y; subtracting y recovers negative (low part of y)
                a[i,j] = t               # Algebraically, c should always be zero. Beware overly-aggressive optimizing compilers!
            end
        end 
    else
        printf("Warning: kahanSum! works only for arrays and matrices.")
    end
    
end