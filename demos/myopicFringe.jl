using JuMP, Gurobi # Optimization packages
using DataFrames, CSV, YAML, DataStructures # dataprocessing
using ProgressBars, Printf # progress bar
using JLD2
using Base.Threads: @spawn 
using Base: split
using ArgParse # Parsing arguments from the command line

include("../src/agents.jl")
include("../src/loadData.jl")
include("../src/Backbone.jl")
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
#GRBsetparam(GUROBI_ENV, "Threads", "8")
#GRBsetparam(GUROBI_ENV, "Method", "2")  
GRBsetparam(GUROBI_ENV, "TimeLimit", "300")  # will only affect solutions if you're selecting representative days  
println("        ")

scenarios = YAML.load_file(joinpath(@__DIR__, "../data/scenarios.yaml"));

sector = "steelmaking"

data = YAML.load_file(joinpath(@__DIR__, "../data/assumptions_agents.yaml"));
define_ETS_parameters!(data)
define_sector_parameters!(data,sector)

#nb = 1
#scenario = scenarios[1]
for (nb, scenario) in scenarios
    # Load Data
 
    dataScen = merge(copy(data),scenario)
    #stoch = Dict()
    #define_stoch_parameters!(stoch,dataScen)
    # Define agents
    agents = Dict()
    agents["fringe"] = build_myopic_competitive_fringe!( Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))), dataScen)
    for (route, dict) in dataScen["sectors"][sector]
        agents[route] = build_myopic_banking_producer!( Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))), dataScen, sector, route)
    end

    results, ADMM = define_results(dataScen,agents) 

    # Solve agents
    ADMM!(results,ADMM,dataScen,sector,agents)

    # Write solution
    sol = get_solution(agents,results)
    CSV.write("results/myopicfringe_"* string(nb) * ".csv",sol)
    #print(sol)
end
