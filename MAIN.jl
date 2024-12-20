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

scenarios = YAML.load_file(joinpath(@__DIR__, "data/scenarios.yaml"));

sector = "steelmaking"

data = YAML.load_file(joinpath(@__DIR__, "data/assumptions_agents.yaml"));
define_ETS_parameters!(data)
define_sector_parameters!(data,sector)
dataScen = merge(copy(data),scenarios[1])
define_stoch_parameters!(dataScen,2)


#################################################################
#           PERFECT FORESIGHT
#################################################################

nb = 1
# Define agents
agents = Dict()
agents["fringe"] = build_stochastic_competitive_fringe!( Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))), dataScen)
for (route, dict) in dataScen["sectors"][sector]
    agents[route] = build_stochastic_producer!( Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))), dataScen, sector, route)
end

results, ADMM = define_results_stochastic(dataScen,agents) 

# Solve agents
ADMM!(results,ADMM,dataScen,sector,agents)

# Write solution
sol = get_solution_summarized(agents,results)
mkpath("results")
CSV.write("results/perfect_foresight_stochastic_"* string(nb) * ".csv",sol)


#################################################################
#           LIQUIDITY CONSTRAINED ETS TRADING
#################################################################

# Define agents
agents = Dict()
agents["fringe"] = build_stochastic_liquidity_constraint_fringe!( Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))), dataScen)
for (route, dict) in dataScen["sectors"][sector]
    agents[route] = build_stochastic_liquidity_constraint_producer!( Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))), dataScen, sector, route)
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

#################################################################
#           MYOPIC ETS TRADING
#################################################################

for (nb, scenario) in scenarios 
    # Load Data
    data_myopic = merge(YAML.load_file(joinpath(@__DIR__, "data/assumptions_agents.yaml")),scenario)
    define_ETS_parameters!(data_myopic)
    define_sector_parameters!(data_myopic,sector)
    define_stoch_parameters!(data_myopic,2)

    # Define agents
    agents_myopic = Dict()
    agents_myopic["fringe"] = build_stochastic_competitive_fringe!( Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))), data_myopic)
    for (route, dict) in data_myopic["sectors"][sector]
        agents[route] = build_stochastic_producer!( Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))), data_myopic, sector, route)
    end

    results_myopic, ADMM_myopic = define_results_stochastic(data_myopic,agents_myopic) 

    # Solve agents
    ADMM_rolling_horizon!(results_myopic,ADMM_myopic,data_myopic,sector,agents_myopic)

    # Write solution
    sol = get_solution_summarized(agents_myopic,results_myopic)
    mkpath("results")
    CSV.write("results/rolling_horizon_stochastic_"* string(nb) * ".csv",sol)
    mkpath("results/detailed")
    sol = get_solution(agents,results_myopic)
    CSV.write("results/detailed/rolling_horizon_stochastic_"* string(nb) * ".csv",sol)
    #print(sol)
end


#################################################################
#           LIQUIDITY CONSTRAINT & MYOPIC ETS TRADING
#################################################################

for (nb, scenario) in scenarios 
    # Load Data
    data_combined = merge(YAML.load_file(joinpath(@__DIR__, "data/assumptions_agents.yaml")),scenario)
    define_ETS_parameters!(data_combined)
    define_sector_parameters!(data_combined,sector)
    define_stoch_parameters!(data_combined,2)

    # Define agents
    agents = Dict()
    agents["fringe"] = build_stochastic_liquidity_constraint_fringe!( Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))), data_combined)
    for (route, dict) in data_combined["sectors"][sector]
        agents[route] = build_stochastic_liquidity_constraint_producer!( Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))), data_combined, sector, route)
    end

    results, ADMM = define_results_stochastic(data_combined,agents) 

    # Solve agents
    ADMM_rolling_horizon!(results,ADMM,data_combined,sector,agents)

    # Write solution
    sol = get_solution_summarized(agents,results)
        mkpath("results")
        CSV.write("results/combination_stochastic_"* string(nb) * ".csv",sol)
        mkpath("results/detailed")
    sol = get_solution(agents,results)
        CSV.write("results/detailed/combination_stochastic_"* string(nb) * ".csv",sol)
    #print(sol)
end