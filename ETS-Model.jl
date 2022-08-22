import Pkg;
#Pkg.add("JuMP")
#Pkg.add("Gurobi")
#Pkg.add("YAML")
#Pkg.add("Plots")
#Pkg.add("DataFrames")


using JuMP
using Gurobi
using YAML
using Plots
using DataFrames
using Plots

using Statistics
using Distributions, Random
# If there exist different carbon reducing technologies with different CAPEX OPEX (efficiencies) what is an optimal policy

# Read ETS data
data = YAML.load_file(joinpath(@__DIR__, "data_assumptions_ets.yaml"));
data_agents = YAML.load_file(joinpath(@__DIR__, "data_assumptions.yaml"));
routes = data_agents["SteelRoutes"]
properties = YAML.load_file(joinpath(@__DIR__, "data_properties.yaml"));

# Parameters/variables ETS 
ETS = Dict()
function define_ETS_parameters!(ETS::Dict,data::Dict)
    # LRF 2017 - 2021, no Brexit, expansion of scope (aviation, maritime) or Green Deal
    ETS["LRF"] = zeros(data["nyears"],1);
    ETS["LRF"][1:4] = data["LRF_2017"]*ones(4,1); 
    ETS["LRF"][5:14] = ETS["LRF"][1]*data["LRF_2021"]/0.0174*ones(10,1);                            # 2021-2030
    ETS["LRF"][15:end] = ETS["LRF"][1]*data["LRF_2031"]/0.0174*ones(data["nyears"]-14,1);           # 2030 - end ETS
       
    # CAP
    ETS["CAP"] = zeros(data["nyears"],1);
    for y =1:data["nyears"]
        ETS["CAP"][y]= maximum([data["CAP_2016"]-sum(ETS["LRF"][1:y]) 0])
    end
    ETS["CAP"][1] = data["S_2017"]+data["EX_2016"] # real supply in 2017 and surplus in market at end 2016

    # Supply 
    ETS["S"] = copy(ETS["CAP"])

    return ETS
end
define_ETS_parameters!(ETS, data)

# Stochastic optimization parameters
stoch = Dict()
function define_stoch_parameters!(stoch::Dict,data::Dict)
    nsamples = data["nsamples"]
    avg = data["MAC"]
    std = data["std"]

    if data["use_cvar"]
        d = Normal(avg,std)
        nselect = floor(Int,data["cvar"]*nsamples)
        stoch["MAC"] = sort(rand(d,nsamples), rev=true)[1:nselect]
        stoch["nsamples"] = size(SO["MAC"])[1]
    else
        stoch["MAC"] = avg
        stoch["nsamples"] = 1
    end
    return stoch
end
define_stoch_parameters!(stoch,data)

function define_sets_parameters!(mod::Model, data::Dict, stoch::Dict) 
    # Define sets
    mod.ext[:sets] = Dict()
    mod.ext[:sets][:Y] = 1:data["nyears"] 
    mod.ext[:sets][:S] = 1:stoch["nsamples"]
    S = mod.ext[:sets][:S]   
    Y = mod.ext[:sets][:Y]

    # define parameters 
    mod.ext[:parameters] = Dict()

    # Emissions representative agents, bound to historical values in 2017-2019
    mod.ext[:parameters][:e] = ones(data["nyears"],1)*data["E_2019"]
    E = mod.ext[:parameters][:e]

    mod.ext[:parameters][:A] = ones(data["nyears"],1)
    for y in 4:data["nyears"]
        mod.ext[:parameters][:A][y] = 1/(1+data["discount_rate"])^(y-3);
    end
    A = mod.ext[:parameters][:A]
    
    # Price structure
    mod.ext[:parameters][:λ] = zeros(data["nyears"],1)

    # Define variables
    mod.ext[:variables] = Dict() 
    buy = mod.ext[:variables][:buy] =  @variable(mod, [y=Y], lower_bound=0, base_name="allowances bougth")
    a = mod.ext[:variables][:a] =  @variable(mod, [y=Y], lower_bound=0, base_name="abatement") # Mton CO2

    # Define expressions
    mod.ext[:expressions] = Dict()
    mod.ext[:expressions][:bank] = @expression(mod, [y=Y], sum(buy[1:y])+sum(a[1:y])-sum(E[1:y]))
    mod.ext[:expressions][:netto_emiss] = @expression(mod, [y=Y], E[y]-a[y])

    MAC = stoch["MAC"] # Marginal abatement cost
    # Objective
    mod.ext[:objective] = @objective(mod, Min, sum(sum(MAC[s]*a[y]^2*A[y] for y in Y) for s in S)/stoch["nsamples"])
 
    # Define constraint
    mod.ext[:constraints] = Dict()
    mod.ext[:constraints][:buycons] = @constraint(mod,[y=Y], buy[y] <= ETS["S"][y])
    mod.ext[:constraints][:auctioncons] = @constraint(mod, [y=Y], sum(buy[1:y]) >= sum(E[1:y]) - sum(a[1:y]))
    mod.ext[:constraints][:nonnegemiss] = @constraint(mod, [y=Y], E[y] - a[y] >= 0)

    return mod
end

function define_sets_parameters!(mod::Model, data::Dict)
    # Alias of define_sets_parameters without stochasticity
    stoch = Dict()
    stoch["nsamples"] = 1
    stoch["MAC"] = [data["MAC"]]
    define_sets_parameters!(mod, data, stoch)
end

# Model Representation of the agents
function build_agent!(agent::Model,data::Dict,route::Dict,stoch::Dict)
    agent.ext[:sets] = Dict()
    agent.ext[:sets][:Y] = 1:data["nyears"]
    agent.ext[:sets][:S] = 1:stoch["nsamples"]
    S = agent.ext[:sets][:S]
    Y = agent.ext[:sets][:Y]

    # define parameters 
    agent.ext[:parameters] = Dict()
    agent.ext[:parameters][:λ_ets] = zeros(data["nyears"],1)
    agent.ext[:parameters][:A] = ones(data["nyears"],1)
    agent.ext[:parameters][:I] = ones(data["nyears"],1)
    for y in 4:data["nyears"]
        agent.ext[:parameters][:A][y] = 1/(1+data["discount_rate"])^(y-3);
        agent.ext[:parameters][:I][y] = (1+data["inflation"])^(y-3);
    end
    A = agent.ext[:parameters][:A]
    I = agent.ext[:parameters][:I]
    λ_ets = agent.ext[:parameters][:λ_ets]
    CAPEX = agent.ext[:parameters][:CAPEX] = data["CAPEX"] # €/MW
    η = agent.ext[:parameters][:η] = data["efficiency_PEM"]
    f = agent.ext[:parameters][:f_steel] = data["f_steel"]

    # Define variables
    agent.ext[:variables] = Dict() 
    cap = agent.ext[:variables][:cap] = @variable(agent, [y=Y], lower_bound=0, base_name="capacity") # MW electricity
    g = agent.ext[:variables][:g] = @variable(agent, [y=Y], lower_bound=0, base_name="generation") # MWh electricity

    # Define expressions
    agent.ext[:expressions] = Dict()
    hydrogen = agent.ext[:expressions][:hydrogen] = @expression(agent, [y=Y], g[y]/η) # MWh hydrogen
    abatement_Mton = agent.ext[:expressions][:abatement_Mton] = @expression(agent, [y=Y], g[y]*η/route["Needs"]["Hydrogen"]*1.850*1e-6) # Mton CO2
    abatement_ton = agent.ext[:expressions][:abatement_ton] = @expression(agent, [y=Y], g[y]*η/route["Needs"]["Hydrogen"]*1.850) # ton CO2
    λ_elec = agent.ext[:expressions][:λ_elec] = @expression(agent, [y=Y], data["p_elek"]*I[y] + g[y]*1/10000*I[y] ) # EUR/MWh

    # Objective
    agent.ext[:objective] = @objective(agent, Min,sum(sum((CAPEX*cap[y] + λ_elec[y]*g[y]/η - λ_ets[y]*abatement_ton[y])*A[y] for y in Y) for s in S)/stoch["nsamples"])
  
    # Define constraint
    agent.ext[:constraints] = Dict()
    agent.ext[:constraints][:capacitycons] = @constraint(agent, [y=Y], sum(cap[1:y])*8760 >= g[y])
    #agent.ext[:constraints][:capcons] = @constraint(agent, [y=Y], cap[y]<= 100)
    agent.ext[:constraints][:maxabatecons] = @constraint(agent, [y=Y], abatement_Mton[y] <= ETS["S"][y]*f)
    return agent
end

function build_agent!(agent::Model,data::Dict,route::Dict)
    stoch = Dict()
    stoch["nsamples"] = 1
    stoch["MAC"] = [data["MAC"]]
    build_agent!(agent, data, route, stoch) 
end

# Linking of the agent models and competivive fringe
function update_prices!(agent::Model,mod::Model)
    # Fetch sets
    S = agent.ext[:sets][:S]
    nb_S = size(S)[1]
    Y = agent.ext[:sets][:Y]

    # Fetch variables
    cap = agent.ext[:variables][:cap] # MW electricity
    g = agent.ext[:variables][:g] # MWh electricity
    a = JuMP.value.(mod.ext[:variables][:a])  # Mton CO2

    # Fetch expressions
    hydrogen = agent.ext[:expressions][:hydrogen]
    abatement_Mton = agent.ext[:expressions][:abatement_Mton]
    abatement_ton = agent.ext[:expressions][:abatement_ton]
    netto_emiss = mod.ext[:expressions][:netto_emiss]
    λ_elec = agent.ext[:expressions][:λ_elec]

    # Fetch Parameters
    A = agent.ext[:parameters][:A]
    I = agent.ext[:parameters][:I]
    CAPEX = agent.ext[:parameters][:CAPEX]
    λ_ets = -[shadow_price(mod.ext[:constraints][:buycons][i]) for i in 1:45]./A[:,1]
    e = JuMP.value.(mod.ext[:expressions][:netto_emiss])
    η = agent.ext[:parameters][:η]
    f = agent.ext[:parameters][:f_steel]

    # Update prices
    mod.ext[:parameters][:λ] = λ_ets
    agent.ext[:parameters][:λ] = λ_ets
    
    # Objective
    agent.ext[:objective] = @objective(agent, Min,sum(sum((CAPEX*cap[y] + λ_elec[y]*g[y]/η- λ_ets[y]*abatement_ton[y])*A[y] for y in Y) for s in S)/nb_S)

    # Adjust max abatement constraint for change competitive fringe
    for y in Y
        set_normalized_rhs(agent.ext[:constraints][:maxabatecons][y], e[y]*f)
    end
    return λ_ets
end    

function update_abatement!(mod::Model, agent::Model)
    # Fetch sets
    S = mod.ext[:sets][:S]   
    Y = mod.ext[:sets][:Y]

    # Fetch parameters 
    E = mod.ext[:parameters][:e]
    A = mod.ext[:parameters][:A]

    # Fetch variables
    buy = mod.ext[:variables][:buy]
    a = mod.ext[:variables][:a]  # Mton CO2

    # Fetch exprsessions
    netto_emiss = mod.ext[:expressions][:netto_emiss]
    abatement_Mton = JuMP.value.(agent.ext[:expressions][:abatement_Mton])

    # Define constraint
    #mod.ext[:constraints][:buycons] = @constraint(mod,[y=Y], buy[y] <= ETS["S"][y] + abatement_Mton[y])
    for y in Y 
        set_normalized_rhs(mod.ext[:constraints][:buycons][y], ETS["S"][y] + abatement_Mton[y])
    end

    # Update bank
    mod.ext[:expressions][:bank] = @expression(mod, [y=Y], sum(buy[1:y])+sum(a[1:y])-sum(E[1:y]))
    return mod
end

##################################################################
#                                                                #
#                       Run Algorithms   	                     #
#                                                                #
##################################################################

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
end 

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
    end

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