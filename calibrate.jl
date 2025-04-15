using JuMP, Gurobi # Optimization packages
using DataFrames, CSV, YAML, DataStructures # dataprocessing
using ProgressBars, Printf # progress bar
using JLD2
using Base.Threads: @spawn 
using Statistics, Random, Distributions
using ArgParse # Parsing arguments from the command line

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
sens_df = CSV.read(joinpath(@__DIR__, "data/sensetivities.csv"), DataFrame)

sens = 1

data = YAML.load_file(joinpath(@__DIR__, "data/assumptions.yaml"))
data["commodityPrices"] =  Dict{Any, Any}(string(k) => sens_df[sens, k] for k in names(sens_df))

nb = 32
# Load Data
dataScen = merge(copy(data), 
# Convert first row into a dictionary with String keys:
Dict(string(k) => scenarios_df[nb, k] for k in names(scenarios_df))
)
global MAC = 0.364230041879668
dataScen["MAC"] = MAC
global E_ref = dataScen["E_ref"]
τ = 1
while abs(τ) > 0.1
    define_ETS_parameters!(dataScen)
    define_sector_parameters!(dataScen)
    define_stoch_parameters!(dataScen,2)
    # Define agents
    agents = Dict()
    agents["trader"] = build_stochastic_liquidity_constraint_trader!( Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))), dataScen)
    agents["fringe"] = build_stochastic_competitive_fringe!( Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))), dataScen)
    for (route, dict) in dataScen["technologies"]
        agents[route] = build_stochastic_producer!( Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))), dataScen, route)
    end

    results, ADMM = define_results_stochastic(dataScen,agents) 

    # Solve agents once
    ADMM[:isRollingHorizon] = true
    ADMM[:start] = 1
    ADMM[:end]  = min(ADMM[:start] + dataScen["horizon_ets"]-1, dataScen["nyears"])
    mask = zeros(data["nyears"])
    mask[ADMM[:start]:ADMM[:end]] .= 1
    ADMM[:mask] = mask

    set_masks!(agents, ADMM, dataScen)

    set_lookahead_window!(agents)
    ADMM!(results, ADMM, dataScen, agents)
    move_lookahead_window!(agents, ADMM)
    ADMM!(results, ADMM, dataScen, agents)

    # Write solution
    sol = get_solution_summarized(agents,results)
    P_2025 = sol[2,:λ_ets_mean]
    global τ = P_2025 - dataScen["P_calibration"]
    global MAC = min(MAC/((1+τ/dataScen["P_calibration"])^(1/2)),100)
    dataScen["MAC"] = MAC
    #global E_ref
    dataScen["E_ref"] = E_ref
    print("Price Error:   "* string(τ)* "  New MAC:  "* string(MAC) * "\n")
end