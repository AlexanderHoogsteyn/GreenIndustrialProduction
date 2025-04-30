using JuMP, Gurobi # Optimization packages
using DataFrames, CSV, YAML, DataStructures # dataprocessing
using ProgressBars, Printf # progress bar
using JLD2
using Base.Threads: @spawn 
using Statistics, Random, Distributions
using ArgParse # Parsing arguments from the command line

include("src/agents.jl")
include("src/loadData.jl")
include("src/Backbone.jl")
include("src/ADMM.jl")
include("src/producer.jl")
include("src/fringe.jl")
include("src/traders.jl")


# Gurobi environment to suppress output
println("Define Gurobi environment...")
println("        ")
const GUROBI_ENV = Gurobi.Env()
# set parameters:
GRBsetparam(GUROBI_ENV, "OutputFlag", "0")   
GRBsetparam(GUROBI_ENV, "TimeLimit", "300")  # will only affect solutions if you're selecting representative days  
println("        ")

# Read the new CSV with scenarios 
scenarios_df = CSV.read(joinpath(@__DIR__, "data/scenarios.csv"), DataFrame)
sens_df = CSV.read(joinpath(@__DIR__, "data/sensetivities.csv"), DataFrame)

sens = 1


data = YAML.load_file(joinpath(@__DIR__, "data/assumptions.yaml"))
for nb in range(71,74)
    # Load Data
    local dataScen = merge(copy(data), 
    # Convert first row into a dictionary with String keys:
    Dict(string(k) => scenarios_df[nb, k] for k in names(scenarios_df)),
    Dict( "commodityPrices" => Dict{Any, Any}(string(k) => sens_df[sens, k] for k in names(sens_df)))
    )
    define_ETS_parameters!(dataScen)
    define_sector_parameters!(dataScen)
    define_stoch_parameters!(dataScen,2)
    
    # Solve agents
    agents, results = ADMM_rolling_horizon!(dataScen)

    # Write solution
    local sol = get_solution_summarized(agents,results)
        mkpath("results")
        CSV.write("results/scenario_deterministic_"* string(nb) * ".csv",sol)
        mkpath("results/detailed")
    local sol_detailed = get_solution(agents,results)
        CSV.write("results/detailed/scenario_deterministic_"* string(nb) * ".csv",sol_detailed)
    if haskey(results, "PriceConvergence")
        CSV.write("results/detailed/ets_prices_"* string(nb) * ".csv",results["PriceConvergence"])
    end

    println(" Scenario ", nb, " solved successfully")
end