import Pkg;
#Pkg.add("JuMP")
#Pkg.add("Gurobi")
#Pkg.add("YAML")
#Pkg.add("Plots")
#Pkg.add("DataFrames")
Pkg.instantiate()

using JuMP
using Gurobi
using YAML
using DataFrames
using CSV

include("EmissionMarket.jl")
include("ElectricityMarket.jl")
include("ADMM.jl")

frac_digit = 4

function run_ets_model()
    # Initialize model
    mod = Model(optimizer_with_attributes(Gurobi.Optimizer))
    define_sets_parameters!(mod,data)

    # Solve
    optimize!(mod)
    @show objective_value(mod);

    # Fetch solution
    buy_opt = convert(Array, JuMP.value.(mod.ext[:variables][:buy]))
    abate_opt = convert(Array, JuMP.value.(mod.ext[:variables][:a]))
    bank_opt = convert(Array, JuMP.value.(mod.ext[:expressions][:bank]))
    emission_opt = mod.ext[:parameters][:e][:,1] - abate_opt 
    supply = ETS["S"][:,1]
    prices =  -[shadow_price(mod.ext[:constraints][:buycons][i]) for i in 1:45]./mod.ext[:parameters][:A][:,1]
    
    # Print results
    sol = DataFrame(Y=2019:2063,Supply=supply,Buy=buy_opt, Abate=abate_opt, Emmit=emission_opt,Bank=bank_opt,Price=prices);
    print(sol)
    # Plot results
    plot(Matrix(sol)[:,2:end],labels=permutedims(names(sol[:,2:end])))
    return sol
end;

function run_stochastic_ets_model()
    # Initialize model
    mod = Model(optimizer_with_attributes(Gurobi.Optimizer))
    define_sets_parameters!(mod,data, stoch)

    # Solve
    optimize!(mod)
    @show objective_value(mod);

    # Fetch solution
    buy_opt = convert(Array, JuMP.value.(mod.ext[:variables][:buy]))
    abate_opt = convert(Array, JuMP.value.(mod.ext[:variables][:a]))
    bank_opt = convert(Array, JuMP.value.(mod.ext[:expressions][:bank]))
    emission_opt = mod.ext[:parameters][:e][:,1] - abate_opt 
    supply = ETS["S"][:,1]
    prices =  -[shadow_price(mod.ext[:constraints][:buycons][i]) for i in 1:45]./mod.ext[:parameters][:A][:,1]
    
    # Print results
    sol = DataFrame(Y=2019:2063,Supply=supply,Buy=buy_opt, Abate=abate_opt, Emmit=emission_opt,Bank=bank_opt,Price=prices);
    
    # Plot results
    plot(Matrix(sol)[:,2:end],labels=permutedims(names(sol[:,2:end])))
    return sol
end;

function run_joint_ets_model_iter()
    # Initial solve model
    mod = Model(optimizer_with_attributes(Gurobi.Optimizer))
    define_sets_parameters!(mod,data)
        # For one agent which can use H-DRI
    optimize!(mod);

    # Set up agent
    agent = Model(optimizer_with_attributes(Gurobi.Optimizer))
    route = routes["HDRI-EAF"]
    build_agent!(agent,data,route)

    prices = update_prices!(agent,mod);
    optimize!(agent);
    update_abatement!(mod,agent);
    optimize!(mod);
    updated_prices = update_prices!(agent,mod);
    prices = hcat(prices,updated_prices)

    δ = sum((prices[:,end-1] - prices[:,end]).^2)
    iter = 1

    #Keep solving iteratively
    while δ  > 1e-6
        optimize!(agent);
        update_abatement!(mod,agent);
        optimize!(mod);
        updated_prices = update_prices!(agent,mod);
        prices = hcat(prices,updated_prices)
        δ = sum((prices[:,end-1] - prices[:,end]).^2)
        iter += 1 
        if iter == 100
            print("Did not converge after " * string(iter) * " runs     δ= " * string(δ) * "\n")
            print("--------------------------------------------------\n")
            break
        end
        print("Run: " * string(iter) * "     δ= " * string(δ) * "\n")
        print("--------------------------------------------------\n")
    end;

    # final optimization to circonvent OptimizeNotCalled()
    optimize!(mod);
    optimize!(agent);

    frac_digit = 4
    buy_opt = round.(convert(Array, JuMP.value.(mod.ext[:variables][:buy])),digits=frac_digit)
    abate_opt = round.(convert(Array, JuMP.value.(mod.ext[:variables][:a])),digits=frac_digit)
    bank_opt = round.(convert(Array, JuMP.value.(mod.ext[:expressions][:bank])),digits=frac_digit)
    emission_opt = round.(mod.ext[:parameters][:e][:,1] - abate_opt ,digits=frac_digit)
    supply = round.(ETS["S"][:,1],digits=frac_digit)

    cap_opt = round.(convert(Array, JuMP.value.(agent.ext[:variables][:cap])),digits=frac_digit) # MW electricity
    g_opt = round.(convert(Array, JuMP.value.(agent.ext[:variables][:g])),digits=frac_digit) # MWh electricity
    abatement_Mton = round.(convert(Array, JuMP.value.(agent.ext[:expressions][:abatement_Mton])),digits=frac_digit)

    # Print results
    sol = DataFrame(Y=2019:2063,Buy=buy_opt, Abate=abate_opt, Emmit=emission_opt,Bank=bank_opt,InitPrice=prices[:,1], PriceUpdate=prices[:,end], Capacity=cap_opt,Agent=abatement_Mton);
    return sol
end;

function get_solution(mod::Model, agents::Dict)
    buy_opt = round.(convert(Array, JuMP.value.(mod.ext[:variables][:buy])),digits=frac_digit)
    abate_opt = round.(convert(Array, JuMP.value.(mod.ext[:variables][:a])),digits=frac_digit)
    bank_opt = round.(convert(Array, JuMP.value.(mod.ext[:expressions][:bank])),digits=frac_digit)
    emission_opt = round.(mod.ext[:parameters][:e][:,1] - abate_opt ,digits=frac_digit)
    supply = round.(ETS["S"][:,1],digits=frac_digit)
    sol = DataFrame(Y=2019:2063,buy=buy_opt, abate=abate_opt, emmit=emission_opt,bank=bank_opt);

    for (key,agent) in agents
        for variable in keys(agent.ext[:variables])
            values = round.(convert(Array, JuMP.value.(agent.ext[:variables][variable])),digits=frac_digit)
            variable_name = Symbol(string(key) * "_" * string(variable))
            sol[!,variable_name] = values 
        end
    end
    return sol
end

scenarios = YAML.load_file(joinpath(@__DIR__, "../data/scenarios.yaml"));
 

for (nb, scenario) in scenarios
    mod = Model(optimizer_with_attributes(Gurobi.Optimizer))
    define_sets_parameters!(mod,data)

    agents = Dict()
    agent_alk = Model(optimizer_with_attributes(Gurobi.Optimizer))
    agents["Alkaline"] = build_alkaline_agent!(agent_alk,data,"SMR")
    agent_PEM = Model(optimizer_with_attributes(Gurobi.Optimizer))
    agents["PEM"] = build_PEM_agent!(agent_PEM,data,"SMR")
    ADMM(mod,agents,0,1e-3,scenario)

    sol = get_solution(mod,agents)
    CSV.write("results\\scenario_"* string(nb) * ".csv",sol)
    print(sol)
end