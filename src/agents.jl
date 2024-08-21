import Pkg;
#Pkg.add("JuMP")
#Pkg.add("Gurobi")
#Pkg.add("YAML")
#Pkg.add("Plots")
#Pkg.add("DataFrames")


using JuMP
using YAML
using DataFrames

using Statistics
using Distributions, Random
# If there exist different carbon reducing technologies with different CAPEX OPEX (efficiencies) what is an optimal policy

# Parameters/variables ETS 
function define_ETS_parameters!(data::Dict)
    # LRF 2017 - 2021, no Brexit, expansion of scope (aviation, maritime) or Green Deal
    data["LRF"] = zeros(data["nyears"],1);
    data["LRF"][1:4] = data["LRF_2017"]*ones(4,1); 
    data["LRF"][5:14] = data["LRF"][1]*data["LRF_2021"]/0.0174*ones(10,1);                            # 2021-2030
    data["LRF"][15:end] = data["LRF"][1]*data["LRF_2031"]/0.0174*ones(data["nyears"]-14,1);           # 2030 - end ETS
       
    # CAP
    data["CAP"] = zeros(data["nyears"],1);
    for y =1:data["nyears"]
        data["CAP"][y]= maximum([data["CAP_2016"]-sum(data["LRF"][1:y]) 0])
    end
    data["CAP"][1] = data["S_2017"]+data["EX_2016"] # real supply in 2017 and surplus in market at end 2016

    # Supply 
    data["S"] = copy(data["CAP"])

    return data
end

function define_sector_parameters!(data::Dict,route::String)
    # Demand good
    data["D"] = ones(data["nyears"],1)*data["demand"][route];

    return data
end


# Stochastic optimization parameters
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

# Model Representation of the agents
function build_agent!(agent::Model, data::Dict)
    # Build structure of agents which is shared among all types
    agent.ext[:sets] = Dict()
    Y = agent.ext[:sets][:Y] = 1:data["nyears"]

    # define parameters 
    agent.ext[:parameters] = Dict()
    agent.ext[:parameters][:λ_ets] = zeros(data["nyears"],1)
    agent.ext[:parameters][:λ_elec] = zeros(data["nyears"],1)
    agent.ext[:parameters][:λ_good] = zeros(data["nyears"],1)
    agent.ext[:parameters][:A] = ones(data["nyears"],1)
    agent.ext[:parameters][:I] = ones(data["nyears"],1)
    for y in 4:data["nyears"]
        agent.ext[:parameters][:A][y] = 1/(1+data["WACC"])^(y-3);
        agent.ext[:parameters][:I][y] = (1+data["inflation"])^(y-3);
    end

    # define expressions
    agent.ext[:expressions] = Dict()

    # define variables
    agent.ext[:variables] = Dict()
    b = agent.ext[:variables][:b] =  @variable(agent, [y=Y], lower_bound=0, base_name="allowances bougth")

    # Define constraint
    agent.ext[:constraints] = Dict()

    # Objective
    agent.ext[:objective] = @objective(agent, Min, 0)
end

function build_agent!(agent::Model,data::Dict,stoch::Dict)
    build_agent!(agent,data)
    agent.ext[:sets][:S] = 1:stoch["nsamples"]
end

function build_competitive_fringe!(agent::Model, data::Dict, stoch::Dict) 
    # Define sets
    build_agent!(agent,data,stoch)
    A = agent.ext[:parameters][:A]
    S = agent.ext[:sets][:S]   
    Y = agent.ext[:sets][:Y]
    λ_ets = agent.ext[:parameters][:λ_ets]

    # Emissions representative agents, bound to historical values in 2017-2019
    E = agent.ext[:parameters][:e] = ones(data["nyears"],1)*data["E_2019"]

    # Define variables
    a = agent.ext[:variables][:a] =  @variable(agent, [y=Y], lower_bound=0, base_name="abatement") # Mton CO2
    b = agent.ext[:variables][:b]

    # Define expressions
    agent.ext[:expressions][:bank] = @expression(agent, [y=Y], sum(b[1:y])+sum(a[1:y])-sum(E[1:y]))
    agent.ext[:expressions][:netto_emiss] = @expression(agent, [y=Y], E[y]-a[y])
 
    # Define constraint
    agent.ext[:constraints][:buycons] = @constraint(agent,[y=Y], b[y] <= 1.1*data["S"][y])
    agent.ext[:constraints][:auctioncons] = @constraint(agent, [y=Y], sum(b[1:y]) >= sum(E[1:y]) - sum(a[1:y]))
    agent.ext[:constraints][:nonnegemiss] = @constraint(agent, [y=Y], E[y] - a[y] >= 0)

    # zero production
    g = agent.ext[:variables][:g] = @variable(agent, [y=Y], lower_bound=0, base_name="production") # ton product
    agent.ext[:constraints][:zerogen] = @constraint(agent, [y=Y], g[y] == 0)


    return agent
end

function build_competitive_fringe!(agent::Model, data::Dict)
    # Alias of define_sets_parameters without stochasticity
    stoch = Dict()
    stoch["nsamples"] = 1
    stoch["MAC"] = [data["MAC"]]
    build_competitive_fringe!(agent, data, stoch)
end

function solve_competitive_fringe!(agent::Model, data::Dict,stoch::Dict)
    # update Objective
    A = agent.ext[:parameters][:A]
    S = agent.ext[:sets][:S]
    Y = agent.ext[:sets][:Y]
    MAC = stoch["MAC"] # Marginal abatement cost
    λ_ets = agent.ext[:parameters][:λ_ets] 
    a = agent.ext[:variables][:a]
    b = agent.ext[:variables][:b]
    agent.ext[:objective] = @objective(agent, Min, 
                                        sum(A[y]*(MAC[s]*a[y]^2 + λ_ets*b[y]) for y in Y, s in S)
                                        + sum(ρ_EUA/2*(b[y]-b_bar[y])^2 for y in Y, s in S)
    )
    optimize!(agent::Model)
    return agent
end

function solve_competitive_fringe!(agent::Model, data::Dict)
    stoch = Dict()
    stoch["nsamples"] = 1
    stoch["MAC"] = [data["MAC"]]
    build_competitive_fringe!(agent, data, stoch)
    optimize!(agent::Model)
    return agent
end

function build_producer!(agent::Model,data::Dict,sector::String,route::String)
    build_agent!(agent,data,stoch)
    A = agent.ext[:parameters][:A]
    S = agent.ext[:sets][:S]   
    Y = agent.ext[:sets][:Y]
    OPEX = agent.ext[:parameters][:OPEX] = ones(data["nyears"])*route_costs(data["commodityPrices"], data["sectors"][sector][route])
    CAPEX = agent.ext[:parameters][:CAPEX] = ones(data["nyears"])*data["sectors"][sector][route]["CAPEX"]
    EF = agent.ext[:parameters][:EF] = data["sectors"][sector][route]["EUA"]

    λ_ets = agent.ext[:parameters][:λ_ets]
    λ_good = agent.ext[:parameters][:λ_good]

    b = agent.ext[:variables][:b]
    cap = agent.ext[:variables][:cap] = @variable(agent, [y=Y], lower_bound=0, base_name="capacity") # ton/y production capacity
    g = agent.ext[:variables][:g] = @variable(agent, [y=Y], lower_bound=0, base_name="production") # ton product

    agent.ext[:expressions][:netto_emiss] = @expression(agent, [y=Y], g[y]*EF)

    agent.ext[:constraints][:capacitycons] = @constraint(agent, [y=Y], sum(cap[1:y]) >= g[y])
    agent.ext[:constraints][:buycons] = @constraint(agent,[y=Y], sum(b[1:y]) >= sum(g[1:y]*EF))
    return agent
end

function solve_producer!(agent::Model,data::Dict,sector::String,route::String)
    A = agent.ext[:parameters][:A]
    S = agent.ext[:sets][:S]   
    Y = agent.ext[:sets][:Y]
    OPEX = agent.ext[:parameters][:OPEX]
    CAPEX = agent.ext[:parameters][:CAPEX]

    cap = agent.ext[:variables][:cap]
    #cap_bar = agent.ext[:parameters][:cap_bar]
    b = agent.ext[:variables][:b]
    b_bar = agent.ext[:parameters][:b_bar]
    g = agent.ext[:variables][:g]
    g_bar = agent.ext[:parameters][:g_bar]
    ρ_EUA = agent.ext[:parameters][:ρ_EUA]
    ρ_product = agent.ext[:parameters][:ρ_product]

    λ_ets = agent.ext[:parameters][:λ_ets]
    λ_product = agent.ext[:parameters][:λ_product]
    agent.ext[:objective] = @objective(agent, Min,
                            sum(A[y]*(CAPEX[y]*cap[y] + λ_ets[y]*b[y] + (OPEX[y]-λ_product[y])*g[y]) for y in Y, s in S)
                            + sum(ρ_EUA/2*(b[y]-b_bar[y])^2 for y in Y, s in S)
                            + sum(ρ_product/2*(g[y]-g_bar[y])^2 for y in Y, s in S)
                            #+ sum(ρ_cap/2*(cap[y]-cap_bar[y])^2 for y in Y, s in S)
                            )
    optimize!(agent::Model)
    return agent
end