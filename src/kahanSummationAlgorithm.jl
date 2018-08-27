### kahanSumOneValue ###
# DESCRIPTION: a+b is done, where the lost information in that sum is saved in c 
# INPUT: a, b, c
# @param a: a > b should be true for the appropiate work of the algorithm
# @param b: the value that will be added to a
# @param c: inititally 0, the information that it is lost in a+b should be here
# @return sum: a+b 
# @return c: information lost in the sum a+b
function kahanSumOneValue(a, b, c)
    sum = a
    y = b - c          # So far, so good: c is zero.
    t = sum + y        # Alas, sum is big, y small, so low-order digits of y are lost.
    c = (t - sum) - y  # (t - sum) cancels the high-order part of y; subtracting y recovers negative (low part of y)
    sum = t            # Algebraically, c should always be zero. Beware overly-aggressive optimizing compilers!
    return sum, c 
end




### kahanSum! ###
# DESCRIPTION: a+b is done, where the lost information in that sum is saved in c 
# INPUT: a, b, c
# @param a: a > b should be true for the appropiate work of the algorithm
# @param b: the value that will be added to a
# @param c: inititally 0, the information that it is lost in a+b should be here
# @return sum: a+b 
# @return c: information lost in the sum a+b
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