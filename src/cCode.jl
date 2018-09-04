"""
    cKeplerFlow!(x_bek::Array{Float64,1}, v_bek::Array{Float64,1}, x0_bek::Array{Float64,1}, v0_bek::Array{Float64,1}, dt::Any, mu::Float64, traza::Int32)

Defines the position and speeds of a particle -after a time-step- orbiting another body given its initial position and speeds. It does the same as `keplerFlow` function, but it is faster (because it is written in C).

# Args

* `x_bek`: the position vector of the body after dt time from the initial values will be written in this array of D elements, where D is the number of used dimensions.
* `v_bek`: the velocity vector of the body after dt time from the initial values will be written in this array of D elements.
* `x0_bek`: position vector of the secondary body orbiting the primary body.
* `v0_bek`: speed vector of the secondary body orbiting the primary body.
* `dt`: time-step, the function returns the new position and speed vectors after dt time has passed.
* `mu`: standard gravitational parameter of the two bodies.
* `r0`: initial distance between the two bodies.
* `Gz`: array of the first four g-functions. Instead of defining it on each execution of the function, the same array is used in the whole simulation. That is why it's passed as an argument.
* `trace`: boolean value that shows information of the execution if true.
"""
cKeplerFlow!(x_bek::Array{Float64,1}, v_bek::Array{Float64,1}, x0_bek::Array{Float64,1}, v0_bek::Array{Float64,1}, dt::Any, mu::Float64, traza::Int32) = open(string(Base.source_path()[1:(end-8)],"kepler.c")) do f

    code = read(f, String)
    
    # liburutegiaren PATH-a hemen gordeko da
    Clib = mktemp()[1]

    # konpilatzeko, C-ko math liburutegia behar da
    # liburutegia konpilatuko dugu, -lm
    open(`gcc -fPIC -O3 -xc -shared -lm -o $(Clib * "." * Libdl.dlext) -`, "w") do f
         print(f, code)
    end

    # definitu Julia funtzio bat, C funtzioa deitzen duena
    cKeplerFlow!(x_bek::Array{Float64,1}, v_bek::Array{Float64,1}, x0_bek::Array{Float64,1}, v0_bek::Array{Float64,1}, dt::Any, mu::Float64, traza::Int32) = 
        ccall(("kepler", Clib), Nothing,(Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Float64,Float64,Int32), 
        x_bek, v_bek, x0_bek, v0_bek, convert(Float64,dt), mu, traza)
end