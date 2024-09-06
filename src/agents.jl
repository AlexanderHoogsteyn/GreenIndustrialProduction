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
    data["LRF"][1:4] = data["LRF_2024"]/100*ones(4,1);                            # 2021-2030
    data["LRF"][5:end] = data["LRF_2028"]/100*ones(data["nyears"]-4,1);           # 2030 - end ETS
       
    # CAP
    data["CAP"] = zeros(data["nyears"],1);
    for y =1:data["nyears"]
        data["CAP"][y]= maximum([data["CAP_2024"]*(1-sum(data["LRF"][1:y])) 0])
    end
    data["CAP"][1] = data["CAP_2024"] #+data["EX_2016"] # real supply in 2017 and surplus in market at end 2016

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
    agent.ext[:sets][:S] = 1:1

    # define parameters 
    agent.ext[:parameters] = Dict()
    agent.ext[:parameters][:λ_ets] = zeros(data["nyears"],1)
    agent.ext[:parameters][:λ_elec] = zeros(data["nyears"],1)
    agent.ext[:parameters][:λ_product] = zeros(data["nyears"],1)
    agent.ext[:parameters][:A] = ones(data["nyears"],1)
    agent.ext[:parameters][:i] = ones(data["nyears"],1)
    agent.ext[:parameters][:r_equity] = ones(data["nyears"],1)
    agent.ext[:parameters][:r_debt] = ones(data["nyears"],1)
    for y in 1:data["nyears"]
        agent.ext[:parameters][:A][y] = 1/(1+data["discount_rate"])^(y);
        agent.ext[:parameters][:r_debt][y] = 1/(1+data["debt_rate"])^(y);
        agent.ext[:parameters][:r_equity][y] = 1/(1+data["equity_rate"])^(y);
        agent.ext[:parameters][:i][y] = (1+data["inflation"])^(y);
    end

    # define expressions
    agent.ext[:expressions] = Dict()

    # define variables
    agent.ext[:variables] = Dict()
    b = agent.ext[:variables][:b] =  @variable(agent, [y=Y], base_name="allowances bougth")

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
    Y = agent.ext[:sets][:Y]

    # Emissions representative agents, bound to historical values in 2017-2019
    E = agent.ext[:parameters][:e] = zeros(data["nyears"],1)
    E[1] = data["E_2019"]

    # Define variables
    b = agent.ext[:variables][:b]

    # Define expressions
    agent.ext[:expressions][:bank] = @expression(agent, [y=Y], sum(b[1:y])-sum(E[1:y]))
    agent.ext[:expressions][:netto_emiss] = @expression(agent, [y=Y], E[y])
 
    # Define constraint
    agent.ext[:constraints][:con1]  = @constraint(agent,[y=Y], sum(b[1:y]) >= sum(E[1:y]))

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

function build_myopic_competitive_fringe!(agent::Model,data::Dict)
    build_competitive_fringe!(agent,data)

    agent.ext[:parameters][:isMyopic] = true

    Y = agent.ext[:sets][:Y][1:end-data["horizon_ets"]]
    E = agent.ext[:parameters][:e]
    b = agent.ext[:variables][:b]


    # Add constraint on look-ahead of banking horizon
    agent.ext[:constraints][:myopic_banking] = @constraint(agent,[y=Y], 0 <= 0)
    return agent
end

function solve_competitive_fringe!(agent::Model)
    # update Objective
    A = agent.ext[:parameters][:A]
    S = agent.ext[:sets][:S]
    Y = agent.ext[:sets][:Y]
    λ_ets = agent.ext[:parameters][:λ_ets] 
    ρ_ets = agent.ext[:parameters][:ρ_ets]
    E = agent.ext[:parameters][:e]
    b = agent.ext[:variables][:b]
    b_bar = agent.ext[:parameters][:b_bar]

    agent.ext[:objective] = @objective(agent, Min, 
                                        sum(A[y]*λ_ets[y]*b[y] for y in Y, s in S)
                                        + sum(A[y]*ρ_ets/2*(b[y]-b_bar[y])^2 for y in Y, s in S)
    )
    # Update constraints 
    for y in Y
        delete(agent,agent.ext[:constraints][:con1][y])
    end 

    agent.ext[:constraints][:con1]  = @constraint(agent,[y=Y], sum(b[1:y]) >= sum(E[1:y]))

    agent.ext[:expressions][:bank] = @expression(agent, [y=Y], sum(b[1:y])-sum(E[1:y]))
    agent.ext[:expressions][:netto_emiss] = @expression(agent, [y=Y], E[y])

    optimize!(agent::Model)
    return agent
end

function solve_myopic_competitive_fringe!(agent::Model)
    @assert is_myopic(agent) "Agent is not myopic"

    # Update constraints 
    Y = agent.ext[:sets][:Y][1:end-data["horizon_ets"]]
    b = agent.ext[:variables][:b]
    E = agent.ext[:parameters][:e]

    for y in Y
        delete(agent,agent.ext[:constraints][:myopic_banking][y])
    end 
    agent.ext[:constraints][:myopic_banking] = @constraint(agent,[y=Y], sum(b[1:y])-sum(E[1:y]) <= sum(E[y+1:y+data["horizon_ets"]]))
    solve_competitive_fringe!(agent)
end

#function solve_competitive_fringe!(agent::Model)
#    stoch = Dict()
#    stoch["nsamples"] = 1
#    solve_competitive_fringe!(agent, stoch)
#    optimize!(agent::Model)
#    return agent
#end

function build_producer!(agent::Model,data::Dict,sector::String,route::String)
    build_agent!(agent,data)

    @assert haskey(data, "nyears") "Data must contain key 'nyears'"
    @assert haskey(data, "commodityPrices") "Data must contain key 'commodityPrices'"
    @assert haskey(data, "sectors") "Data must contain key 'sectors'"
    @assert haskey(data["sectors"], sector) "Data must contain the given sector"
    @assert haskey(data["sectors"][sector], route) "Sector data must contain the given route"

    Y = agent.ext[:sets][:Y]
    agent.ext[:parameters][:OPEX] = ones(data["nyears"])*route_costs(data["commodityPrices"], data["sectors"][sector][route])
    agent.ext[:parameters][:CAPEX] = ones(data["nyears"])*data["sectors"][sector][route]["CAPEX"]
    EF = agent.ext[:parameters][:EF] = data["sectors"][sector][route]["ETS"]

    legacy_cap = ones(data["nyears"]).*range(data["sectors"][sector][route]["legacy_capacity"],0,data["nyears"])

    b = agent.ext[:variables][:b]
    cap = agent.ext[:variables][:cap] = @variable(agent, [y=Y], lower_bound=0, base_name="capacity") # ton/y production capacity
    g = agent.ext[:variables][:g] = @variable(agent, [y=Y], lower_bound=0, base_name="production") # ton product

    agent.ext[:expressions][:netto_emiss] = @expression(agent, [y=Y], g[y]*EF)
    agent.ext[:expressions][:bank] = @expression(agent, [y=Y], sum(b[1:y])-sum(g[1:y]*EF))

    agent.ext[:constraints][:capacitycons] = @constraint(agent, [y=Y], legacy_cap[y] + sum(cap[1:y]) >= g[y])

    # Allow banking:
    agent.ext[:constraints][:buycons] = @constraint(agent,[y=Y], sum(b[1:y]) >= sum(g[1:y]*EF))
    
    return agent
end

function build_myopic_banking_producer!(agent::Model,data::Dict,sector::String,route::String)
    build_producer!(agent,data,sector,route)

    Y = agent.ext[:sets][:Y]
    b = agent.ext[:variables][:b]
    E = agent.ext[:expressions][:netto_emiss]

    agent.ext[:constraints][:myopic_banking] = @constraint(agent,[y=Y[1:end-data["horizon_ets"]]], sum(b[1:y])-sum(E[1:y]) <= sum(E[y+1:y+data["horizon_ets"]]))

    return agent
end

function build_myopic_producer!(agent::Model,data::Dict,sector::String,route::String)
    build_producer!(agent,data,sector,route)
    Y = agent.ext[:sets][:Y]
    cap = agent.ext[:variables][:cap] 
    g = agent.ext[:variables][:g] 

    # Define myopic specific variables
    agent.ext[:parameters][:isMyopic] = true
    agent.ext[:variables_myopic] = Dict()
    g_τ = agent.ext[:variables_myopic][:g_τ] = @variable(agent, [y=Y,τ=Y], lower_bound=0, base_name="production") # ton product

    agent.ext[:constraints][:gen_total] = @constraint(agent, [y=Y], sum(g_τ[y,:]) == g[y])
    agent.ext[:constraints][:capacitycons_myopic_a] = @constraint(agent, [y in Y, τ in Y; τ > y], g_τ[y,τ] == 0 )
    agent.ext[:constraints][:capacitycons_myopic_b] = @constraint(agent, [y in Y, τ in Y; τ <= y], cap[τ] >= g_τ[y,τ])

    for y in Y
        delete(agent,agent.ext[:constraints][:capacitycons][y])
    end 

    M = agent.ext[:parameters][:M] = ones(data["nyears"],data["nyears"])
    for y in Y
        for τ in Y
            M[y,τ] = (y > τ + data["pay_back_time"] ? 0 : 1)
        end
    end

    E = agent.ext[:expressions][:netto_emiss]
    b = agent.ext[:variables][:b]


    # Add constraint on look-ahead of banking horizon
    #agent.ext[:constraints][:myopic_banking] = @constraint(agent,[y=Y[1:end-data["horizon_ets"]]], sum(b[1:y])-sum(E[1:y]) <= sum(E[y+1:y+data["horizon_ets"]]))
    return agent
end

function solve_producer!(agent::Model)
    @assert !is_myopic(agent) "Agent is myopic"

    A = agent.ext[:parameters][:A]
    r_equity = agent.ext[:parameters][:r_equity]
    r_debt = agent.ext[:parameters][:r_debt]
    i = agent.ext[:parameters][:i]
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
    ρ_ets = agent.ext[:parameters][:ρ_ets]
    ρ_product = agent.ext[:parameters][:ρ_product]

    EF = agent.ext[:parameters][:EF]

    λ_ets = agent.ext[:parameters][:λ_ets]
    λ_product = agent.ext[:parameters][:λ_product]
    agent.ext[:objective] = @objective(agent, Min,
                            sum((r_debt[y]*i[y]*CAPEX[y]*cap[y] + r_equity[y]*λ_ets[y]*b[y] + r_equity[y]*(i[y]*OPEX[y]-λ_product[y])*g[y]) for y in Y, s in S)
                            + sum(r_equity[y]*ρ_ets/2*(b[y]-b_bar[y])^2 for y in Y, s in S)
                            + sum(r_equity[y]*ρ_product/2*(g[y]-g_bar[y])^2 for y in Y, s in S)
                            #+ sum(ρ_cap/2*(cap[y]-cap_bar[y])^2 for y in Y, s in S)
                            )
    optimize!(agent::Model)
    return agent
end

function solve_myopic_producer!(agent::Model)
    # Asserts
    @assert is_myopic(agent) "Agent is not myopic"

    A = agent.ext[:parameters][:A]
    r_equity = agent.ext[:parameters][:r_equity]
    r_debt = agent.ext[:parameters][:r_debt]
    i = agent.ext[:parameters][:i]
    M = agent.ext[:parameters][:M]
    S = agent.ext[:sets][:S]   
    Y = agent.ext[:sets][:Y]
    OPEX = agent.ext[:parameters][:OPEX]
    CAPEX = agent.ext[:parameters][:CAPEX]

    cap = agent.ext[:variables][:cap]
    b = agent.ext[:variables][:b]
    b_bar = agent.ext[:parameters][:b_bar]
    g = agent.ext[:variables][:g]
    g_bar = agent.ext[:parameters][:g_bar]
    g_τ = agent.ext[:variables_myopic][:g_τ]
    g_bar_τ = agent.ext[:parameters][:g_bar_τ]
    ρ_ets = agent.ext[:parameters][:ρ_ets]
    ρ_product = agent.ext[:parameters][:ρ_product]

    λ_ets = agent.ext[:parameters][:λ_ets]
    λ_product = agent.ext[:parameters][:λ_product]

    agent.ext[:objective] = @objective(agent, Min,
    sum((r_debt[y]*i[y]*CAPEX[y]*cap[y] + r_equity[y]*λ_ets[y]*b[y] + sum(r_equity[y]*M[y,τ]*(i[y]*OPEX[y]-λ_product[y])*g_τ[y,τ] for τ in Y)) for y in Y, s in S)
    #sum(A[y]*(CAPEX[y]*cap[y] + λ_ets[y]*b[y] + (OPEX[y]-λ_product[y])*g[y]) for y in Y, s in S)
    + sum(r_equity[y]*ρ_ets/2*(b[y]-b_bar[y])^2 for y in Y, s in S)
    + sum(r_equity[y]*ρ_product/2*(g[y]-g_bar[y])^2 for y in Y, s in S)
    + sum(r_equity[y]*ρ_product/20*(g_τ[y,τ]-g_bar_τ[y,τ])^2 for y in Y, τ in Y, s in S)
    )


    optimize!(agent::Model)
    return agent
end