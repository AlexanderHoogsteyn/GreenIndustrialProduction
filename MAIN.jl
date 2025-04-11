using JuMP, Gurobi # Optimization packages
using DataFrames, CSV, YAML, DataStructures # dataprocessing
using ProgressBars, Printf # progress bar
using JLD2
using Base.Threads: @spawn 
using Statistics, Random, Distributions

include("src/agents.jl")
include("src/loadData.jl")
include("src/backbone.jl")
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

data = YAML.load_file(joinpath(@__DIR__, "data/assumptions.yaml"))
#for nb in range(43,52)
for nb in range(53,61)
    # Load Data
    local dataScen = merge(copy(data), 
    # Convert first row into a dictionary with String keys:
    Dict(string(k) => scenarios_df[nb, k] for k in names(scenarios_df))
    )
    define_ETS_parameters!(dataScen)
    define_sector_parameters!(dataScen)
    define_stoch_parameters!(dataScen,2)
    # Define agents
    local agents = Dict()
    agents["trader"] = build_stochastic_liquidity_constraint_trader!( Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))), dataScen)
    agents["fringe"] = build_stochastic_competitive_fringe!( Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))), dataScen)
    for (route, dict) in dataScen["technologies"]
        agents[route] = build_stochastic_producer!( Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))), dataScen, route)
    end

    local results, ADMM = define_results_stochastic(dataScen,agents) 

    # Solve agents
    ADMM_rolling_horizon!(results,ADMM,dataScen,agents)

    # Write solution
    local sol = get_solution_summarized(agents,results)
        mkpath("results")
        CSV.write("results/scenario_deterministic_"* string(nb) * ".csv",sol)
        mkpath("results/detailed")
    local sol_detailed = get_solution(agents,results)
        CSV.write("results/detailed/scenario_deterministic_"* string(nb) * ".csv",sol_detailed)
    println(" Scenario ", nb, " solved successfully")
end