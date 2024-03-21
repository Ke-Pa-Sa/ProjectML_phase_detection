using SpinMonteCarlo
using Makie
using CairoMakie
using Printf

const model = Ising 
const lat = "square lattice"
const L = 2
const update = SW_update!
const samples = 2

const Tc = 2.0/log1p(sqrt(2))
const Ts = Tc*range(0.85, stop=1.15, length=31)
const MCS = 8192
const Therm = MCS >> 3

#=
for T in Ts
    param = Parameter("Model"=>model, "Lattice"=>lat,
                      "L"=>L, "T"=>T, "J"=>1.0,
                      "Update Method"=>update,
                      "MCS"=>MCS, "Thermalization"=>Therm,
                     )
    result = runMC(param)
    @printf("%f %.15f %.15f\n",
            T, mean(result["Specific Heat"]), stderror(result["Specific Heat"]))
end
=#

param = Parameter("Model"=>model, "Lattice"=>lat,
                      "L"=>L, "T"=>Ts, "J"=>1.0,
                      "Update Method"=>update,
                      "MCS"=>MCS, "Thermalization"=>Therm,
                      "samples"=>samples
                     )

# generates the lattice based on params, for the moment its just a suqare lattice but this can change                     
l = generatelattice(param)
gen_snapshot!(Ising(l), Tc, 3)

const nT = 101
const minT = 0.2
const maxT = 7
temperatures = range(minT, stop=maxT, length=nT) # Temerature range might not me right for all models 
spins_by_temperature = Dict{Float64, Matrix}()
for T in temperatures
    print(T, " ")
    # the 3 is how many images are taken at each temperature, this can be changed
    spins = gen_snapshot!(model(l), T, samples) # Change out Ising to other model if needed, it might need more params like for Potts
    spins_by_temperature[T] = spins
end



for (T, spins) in spins_by_temperature
    print(T, " ")
    for n in 1:size(spins, 2)
        # Reshape the current row into a 16x16 grid
        spin_config = reshape(spins[:, n], L, L)
        
        # Create a heatmap for the current spin configuration
        # colorrange must be changed for different model, for ex (1,q) for Potts and (0,1) for XY
        fig, ax, hm = heatmap(spin_config, colorrange = (-1, 1), colormap = :coolwarm)
        ax.title = "Spin Configuration $n, T = $T"
        display(fig)
    end
end

# first time: uncomment 
#using Pkg
#Pkg.add("CSV")
using CSV

sorted_spins_by_temperature = sort(spins_by_temperature)

open("spins_data_L$(L)_model$(model)_N$(samples)_nT$(nT)_minT$(minT)_maxT$(maxT)_1.csv", "w") do file
    # Write the header row
    write(file, "Temperature,Spins\n")
    
    # Loop through each temperature and its spin configurations
    for (T, spins) in sorted_spins_by_temperature
        println("$T, ", size(spins))
        for n in 1:samples
            println(n)
            println("$T, ", size(spins[:,n]))
            
        # Flatten each spin configuration matrix into a 1D array
        # and join the elements into a comma-separated string
            spins_str = join(spins[:,n], " ")
        
        # Write the temperature and the flattened spin configuration to the file
            write(file, "$T,$spins_str\n")
        end
    end
end