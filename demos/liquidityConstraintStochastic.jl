using JuMP, Gurobi # Optimization packages
using DataFrames, CSV, YAML, DataStructures # dataprocessing
using ProgressBars, Printf # progress bar
using JLD2
using Base.Threads: @spawn 
using Statistics, Random, Distributions

include("../src/agents.jl")
include("../src/loadData.jl")
include("../src/backbone.jl")
include("../src/ADMM.jl")
include("../src/producer.jl")
include("../src/fringe.jl")
include("../src/traders.jl")

# Gurobi environment to suppress output
println("Define Gurobi environment...")
println("        ")
const GUROBI_ENV = Gurobi.Env()
# set parameters:
GRBsetparam(GUROBI_ENV, "OutputFlag", "0")   
GRBsetparam(GUROBI_ENV, "TimeLimit", "300")  # will only affect solutions if you're selecting representative days  
println("        ")

scenarios = YAML.load_file(joinpath(@__DIR__, "../data/scenarios.yaml"));

sector = "steelmaking"

data = YAML.load_file(joinpath(@__DIR__, "../data/assumptions_agents.yaml"));

 
#nb = 1
#scenario = scenarios[nb]
for nb in range(24,24)
    scenario = scenarios[nb]
    # Load Data
    dataScen = merge(copy(data),scenario)
    define_ETS_parameters!(dataScen)
    define_sector_parameters!(dataScen,sector)
    define_stoch_parameters!(dataScen,2)

    # Define agents
    agents = Dict()
    agents["trader"] = build_stochastic_liquidity_constraint_trader!( Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))), dataScen)
    agents["fringe"] = build_stochastic_competitive_fringe!( Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))), dataScen)
    for (route, dict) in dataScen["sectors"][sector]
        agents[route] = build_stochastic_producer!( Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))), dataScen, sector, route)
    end

    results_scen, ADMM = define_results_stochastic(dataScen,agents) 

    # Solve agents
    ADMM_rolling_horizon!(results_scen,ADMM,dataScen,sector,agents)

    # Write solution
    sol = get_solution_summarized(agents,results_scen)
        mkpath("results")
        CSV.write("results/liquidity_constraint_stochastic_"* string(nb) * ".csv",sol)
        mkpath("results/detailed")
    sol = get_solution(agents,results_scen)
        CSV.write("results/detailed/liquidity_constraint_stochastic_"* string(nb) * ".csv",sol)
    #print(sol)
end