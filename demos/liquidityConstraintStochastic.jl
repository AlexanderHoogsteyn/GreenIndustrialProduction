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

scenarios = YAML.load_file(joinpath(@__DIR__, "../data/scenarios_liquidity.yaml"));

sector = "steelmaking"

data = YAML.load_file(joinpath(@__DIR__, "../data/assumptions_agents.yaml"));

 
nb = 6
scenario = scenarios[nb]
#for (nb, scenario) in scenarios 
    # Load Data
    dataScen = merge(copy(data),scenario)
    define_ETS_parameters!(dataScen)
    define_sector_parameters!(dataScen,sector)
    define_stoch_parameters!(dataScen,2)

    # Define agents
    agents = Dict()
    agents["fringe"] = build_stochastic_liquidity_constraint_fringe!( Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))), dataScen)
    for (route, dict) in dataScen["sectors"][sector]
        agents[route] = build_stochastic_producer!( Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))), dataScen, sector, route)
    end

    results, ADMM = define_results_stochastic(dataScen,agents) 

    # Solve agents
    ADMM!(results,ADMM,dataScen,sector,agents)

    # Write solution
    sol = get_solution_summarized(agents,results)
        mkpath("results")
        CSV.write("results/liquidity_constraint_stochastic_"* string(nb) * ".csv",sol)
        mkpath("results/detailed")
    sol = get_solution(agents,results)
        CSV.write("results/detailed/liquidity_constraint_stochastic_"* string(nb) * ".csv",sol)
    #print(sol)
#end