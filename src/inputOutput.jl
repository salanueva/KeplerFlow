### analyzeInput! ###
# DESCRIPTION: analyzes input values of the simulation to assure there is no mistake 
# INPUT: r, v, m, g, t_step, t_max, sld
# @param r: NxD matrix where N is the number of bodies and D the number of dimensions,
#           it contains the initial positions of the bodies
# @param v: NxD matrix where N is the number of bodies and D the number of dimensions,
#           it contains the initial velocities of the bodies
# @param m: array of N elements, each containing the mass of the particle
# @param g: gravitational constant
# @param t_step: time-step of the simulation
# @param t_max: specifies when stops the simulation
# @param sld: simulation logic dictionary, a dictionary where all posible options of the
#             simulation can be specified
#    @key bool_energy: if true, energy of the hamiltonian will be calculated every output step
#    @key bool_print_kep: information about the kepler-flow will be displayed in the console
#    @key ignore_H0: the lineal movement of the whole system will be ignored if this boolean is true.
#    @key out_step: the output will be calculated on every 'out_step'-th step  
# @return exit: if there are some errors, it's true; otherwise, false

function analyzeInput!(r, v, m, g, t_step, t_max, sld)
    
    exit = false
    
    n = size(r)[1]
    if n != size(v)[1] || n != size(m)[1]
        println("Error: Different amount of particles in input data.")
        exit = true
    end
      
    if size(r)[2] != 3 || size(v)[2] != 3
        println("Error: Initial positions and velocities must have 3 coordenates.")
        exit = true
    elseif length(size(m)) != 1
        println("Error: Each particle must have one value as its mass.")
        exit = true
    end
    
    if isnan(g) || g == 0.0
        g = 1.0
        println("Warning: 1 is assigned to the gravitational constant by default.") 
    end
    
    if t_step == 0.0 || isnan(t_step)
        println("Error: t_step can't be 0.")
        exit = true
    elseif t_step < 0.0 && t_max > t_step
        t_max = min(t_step, -t_max)
        println("Warning: t_step and t_max must be both positive or negative.")
    elseif t_step > 0.0 && t_max < t_step
        t_max = min(t_step, -t_max)
        println("Warning: t_step and t_max must be both positive or negative.")
    end
    
    # SIMULATION LOGIC DICTIONARY
    if isa(sld, Dict)
        
        if haskey(sld, "bool_energy")
            if !isa(sld["bool_energy"], Bool)
                sld["bool_energy"] = false
                println("Warning: bool_energy must be boolean. Default: false.")
            end
        else
            sld["bool_energy"] = false # DEFAULT false
        end

        
        if haskey(sld, "bool_print_kep")
            if !isa(sld["bool_print_kep"], Bool)
                sld["bool_print_kep"] = false
                println("Warning: bool_print_kep must be boolean. Default: false.")
            end
        else
            sld["bool_print_kep"] = false # DEFAULT false
        end
        
        if haskey(sld, "ignore_H0")
            if !isa(sld["ignore_H0"], Bool)
                sld["ignore_H0"] = false
                println("Warning: ignore_H0 wasn't well defined. Default: false.")
            end
        else 
            sld["ignore_H0"] = false
            println("Warning: ignore_H0 wasn't defined. Default: false.")
        end
        
        if haskey(sld, "out_step")
            if abs(t_max / sld["out_step"]) < abs(t_step)
                sld["out_step"] = 1
                println("Warning: out_step was out of bounds. All steps will have an output.")
            end
        else 
            sld["out_step"] = 1
            println("Warning: out_step was not defined. All steps will have an output.")
        end

    else
        println("Error: sld must be a dictionary.")
        exit = true
    end
    
    return exit
end




"""
    saveOutput(file::String, r = NaN, v = NaN, e = NaN, time_step=1.0, delimeter::String =" ", digits::Int64 =64)

Saves calculated positions, velocities and/or energies in a .txt file.

# Args

* `file`: the name of the output file (.txt format).
* `r`: OxNxD matrix, where O is the number of outputs, N the number of bodies and D the number of dimensions. It contains the positions of each body on each output time.
* `v`: OxNxD matrix,it contains the velocities of each body on each output time.
* `e`: array of O elements, it contains the energies of the system on each output time.
* `time_step`: specifies the time between two output steps, by default 1.
* `delimeter`: specifies what will be written between values.
* `digits`: specifies, at most, how many decimal digits will be written
"""
function saveOutput(file::String, r = NaN, v = NaN, e = NaN, time_step=1.0, delimeter::String =" ", digits::Int64 =64)

    x = size(r)[1]
    y = size(r)[2]
    z = size(r)[3]
    
    if !((size(v) == () || x == size(v)[1]) && (x == length(e)) || size(e) == ())
        println("Error: first dimensions of r, v and e don't match.")
        return
    elseif !(size(v) == () || (y == size(v)[2] && z == size(v)[3] && z == 3))
        println("Error: r and v's dimensions don't match.")
        return
    end
    
    time = 0.0
    
    f = 0
    if file != ""
        if file[(end-3):end] == ".txt"
            f = open(file, "w")
        else
            println("Warning: filename isn't a txt file. 'output.txt' will be its name.")
            f = open("output.txt")
        end
    else
        return
    end
        
    for m in 1:x
        if size(r) != ()
            for n in 1:y
                write(f, string("$(round(r[m,n,1], digits))", delimeter,"$(round(r[m,n,2], digits))", delimeter,"$(round(r[m,n,3], digits))",delimeter))
            end
            write(f, "\n")
        end
        
        if size(v) != ()
            for n in 1:y
                write(f, string("$(round(v[m,n,1], digits))", delimeter,"$(round(v[m,n,2], digits))", delimeter,"$(round(v[m,n,3], digits))",delimeter))  
            end
            write(f, "\n")
        end
        
        if size(e) != ()
            time += time_step
            write(f, "$time $(e[m])\n")
        end
    end
    
    close(f)
    
end




"""
    showOutput(r, v, e = NaN, names = "", time_step = 1.0)

Shows on the console the outputs that have been calculated.

# Args

* `file`: the name of the output file (.txt format).
* `r`: OxNxD matrix, where O is the number of outputs, N the number of bodies and D the number of dimensions, it contains the positions of each body on each output time.
* `v`: OxNxD matrix,it contains the velocities of each body on each output time.
* `e`: array of O elements, it contains the energies of the system on each output time.
* `names`: array of N strings, it contains the names of each body.
* `time_step`: specifies the time between two output steps, by default 1.0.
"""
function showOutput(r, v, e = NaN, names = "", time_step = 1.0)
    if !(size(r)[1] == size(v)[1] && (size(e) == () || x == length(e)))
        println("Error: first dimensions of r, v and e don't match.")
        return
    elseif size(r)[2] != size(v)[2] || size(r)[3] != size(v)[3] || size(r)[3] != 3
        println("Error: r and v's dimensions don't match.")
        return
    end
    
    total_dt = 0
    if names == ""
        names = string.(collect(1:n))
    end
    for i in 1:size(r)[1]
        total_dt += time_step
        println(" --- Iteration $i - Time $total_dt ---")
        if size(e) != ()
            println(" Total energy: $(e[i])")
        end
        println(" NAME: X, Y, Z; Vx, Vy, Vz")
        for j in 1:size(r)[2]
            println("""$(names[j]): $(r[i,j,1]), $(r[i,j,2]), $(r[i,j,3]); $(v[i,j,1]), $(v[i,j,2]), $(v[i,j,3])\n""")
        end
    end
end