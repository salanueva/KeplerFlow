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
        ccall(("kepler", Clib), Void,(Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Float64,Float64,Int32), 
        x_bek, v_bek, x0_bek, v0_bek, convert(Float64,dt), mu, traza)
end