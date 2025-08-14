## Topic: imperfect behaviour in EU ETS
# Paper: Hoogsteyn A., Meus J., Bruninx K., Delarue E. Barriers to efficient carbon pricing: 
# policy risk, myopic behaviour, and financial constraints. 2025. Working paper.
# Author: alexander Hoogsteyn
# Last update: August 2025
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

data = YAML.load_file(joinpath(@__DIR__, "data/assumptions.yaml"))

# Parse command line arguments
sim_number = parse_commandline()
scen = 11 #sim_number["scen"]
sens = 1# sim_number["sens"]

data["printoutlevel"] = sim_number["printoutlevel"]
data["commodityPrices"] =  Dict{Any, Any}(string(k) => sens_df[sens, k] for k in names(sens_df))

# Load Data
dataScen = merge(copy(data), 
# Convert first row into a dictionary with String keys:
Dict(string(k) => scenarios_df[scen, k] for k in names(scenarios_df)),
)
define_ETS_parameters!(dataScen)
define_sector_parameters!(dataScen)
define_stoch_parameters!(dataScen,2)

# Solve agents
agents, results = ADMM_rolling_horizon!(dataScen)


# Write solution
sol = get_solution_summarized(agents,results)
    mkpath("results")
    CSV.write("results/scenario_"* string(sens)* "_" * string(scen) * ".csv",sol)
    mkpath("results/detailed")
sol_detailed = get_solution(agents,results)
    CSV.write("results/detailed/scenario_"* string(sens)* "_" * string(scen) * ".csv",sol_detailed)
if haskey(results, "PriceConvergence")
    CSV.write("results/detailed/ets_prices_"* string(sens) * "_"* string(scen) * ".csv",results["PriceConvergence"])
end
println(" Scenario ", scen, " solved successfully")
