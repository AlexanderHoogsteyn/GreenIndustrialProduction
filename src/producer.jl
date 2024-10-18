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

    CAP_SV = agent.ext[:parameters][:CAP_SV] = [maximum([0,1-(data["nyears"]-y+1)/data["sectors"][sector][route]["lifetime"]]) for y=1:data["nyears"]] 
    LEG_CAP = agent.ext[:parameters][:LEG_CAP] = zeros(data["nyears"],1)
    agent.ext[:parameters][:LEG_CAP][1] = data["sectors"][sector][route]["legacy_capacity"]  
    agent.ext[:parameters][:LEG_CAP][2:data["nyears"]] = [data["sectors"][sector][route]["legacy_capacity"]*maximum([0,(data["sectors"][sector][route]["legcap_out"]-y+1)/data["sectors"][sector][route]["legcap_out"]]) for y=1:data["nyears"]-1]  
    agent.ext[:parameters][:CAP_LT] = zeros(data["nyears"],data["nyears"]) 
    for y=1:data["nyears"]
        if y+data["sectors"][sector][route]["leadtime"] < data["nyears"]
            for yy = y+data["sectors"][sector][route]["leadtime"]:minimum([y+data["sectors"][sector][route]["leadtime"]+data["sectors"][sector][route]["lifetime"]-1,data["nyears"]])
                agent.ext[:parameters][:CAP_LT][y,yy] = 1
            end
        end
    end
    CAP_LT = agent.ext[:parameters][:CAP_LT] # lead time on new capacity

    b = agent.ext[:variables][:b]
    cap = agent.ext[:variables][:cap] = @variable(agent, [y=Y], lower_bound=0, base_name="capacity") # ton/y production capacity
    g = agent.ext[:variables][:g] = @variable(agent, [y=Y], lower_bound=0, base_name="production") # ton product

    agent.ext[:expressions][:netto_emiss] = @expression(agent, [y=Y], g[y]*EF)
    agent.ext[:expressions][:bank] = @expression(agent, [y=Y], sum(b[1:y])-sum(g[1:y]*EF))
    agent.ext[:expressions][:capacity] = @expression(agent, [y=Y], sum(CAP_LT[y2,y]*cap[y2] for y2=1:y) + LEG_CAP[y] )

    agent.ext[:constraints][:capacitycons] = @constraint(agent, [y=Y], sum(CAP_LT[y2,y]*cap[y2] for y2=1:y) + LEG_CAP[y] >= g[y])

    # Allow banking:
    agent.ext[:constraints][:buycons] = @constraint(agent,[y=Y], sum(b[1:y]) >= sum(g[1:y]*EF))

    return agent
end

function build_stochastic_producer!(agent::Model,data::Dict,sector::String,route::String)
    build_stochastic_agent!(agent,data)

    @assert is_stochastic(agent) "Agent is not stochastic"
    @assert haskey(data, "nyears") "Data must contain key 'nyears'"
    @assert haskey(data, "commodityPrices") "Data must contain key 'commodityPrices'"
    @assert haskey(data, "sectors") "Data must contain key 'sectors'"
    @assert haskey(data["sectors"], sector) "Data must contain the given sector"
    @assert haskey(data["sectors"][sector], route) "Sector data must contain the given route"

    Y = agent.ext[:sets][:Y]
    S = agent.ext[:sets][:S] 
    agent.ext[:parameters][:OPEX] = ones(data["nyears"])*route_costs(data["commodityPrices"], data["sectors"][sector][route])
    agent.ext[:parameters][:CAPEX] = ones(data["nyears"])*data["sectors"][sector][route]["CAPEX"]
    EF = agent.ext[:parameters][:EF] = data["sectors"][sector][route]["ETS"]

    legacy_cap = ones(data["nyears"]).*range(data["sectors"][sector][route]["legacy_capacity"],0,data["nyears"])

    b = agent.ext[:variables][:b]
    cap = agent.ext[:variables][:cap] = @variable(agent, [y=Y], lower_bound=0, base_name="capacity") # ton/y production capacity
    g = agent.ext[:variables][:g] = @variable(agent, [y=Y,s=S], lower_bound=0, base_name="production") # ton product

    agent.ext[:expressions][:netto_emiss] = @expression(agent, [y=Y,s=S], g[y,s]*EF)
    agent.ext[:expressions][:bank] = @expression(agent, [y=Y,s=S], sum(b[1:y,s])-sum(g[1:y,s]*EF))

    agent.ext[:constraints][:capacitycons] = @constraint(agent, [y=Y,s=S], legacy_cap[y] + sum(cap[1:y]) >= g[y,s])

    # Allow banking:
    agent.ext[:constraints][:buycons] = @constraint(agent,[y=Y,s=S], sum(b[1:y,s]) >= sum(g[1:y,s]*EF))

    return agent
end

function build_producer_trader!(agent::Model,data::Dict,sector::String,route::String)
    build_producer!(agent,data,sector,route)

    Y = agent.ext[:sets][:Y]
    b = agent.ext[:variables][:b]
    g = agent.ext[:variables][:g]
    EF = agent.ext[:parameters][:EF]

    if haskey(agent.ext[:constraints], :buycons)
        delete.(agent,agent.ext[:constraints][:buycons])
    end
    # Allow borrowing 
    agent.ext[:constraints][:buycons] = @constraint(agent, sum(b) >= sum(g)*EF)
    agent.ext[:constraints][:buycons_supply] = @constraint(agent,[y=Y], b[y] <= data["CAP"][y])
    
    return agent
end

function build_stochastic_producer_trader!(agent::Model,data::Dict,sector::String,route::String)
    build_stochastic_producer!(agent,data,sector,route)

    # TO DO  : Implement

    return agent
end

function build_myopic_banking_producer!(agent::Model,data::Dict,sector::String,route::String)
    build_producer!(agent,data,sector,route)

    Y = agent.ext[:sets][:Y]
    b = agent.ext[:variables][:b]
    E = agent.ext[:expressions][:netto_emiss]
    horizon_ets = agent.ext[:parameters][:horizon_ets] = data["horizon_ets"]


    agent.ext[:constraints][:myopic_banking] = @constraint(agent,[y=Y[1:end-data["horizon_ets"]]], sum(b[1:y])-sum(E[1:y]) <= sum(E[y+1:y+horizon_ets]))

    return agent
end

function build_stochastic_myopic_banking_producer!(agent::Model,data::Dict,sector::String,route::String)
    build_stochastic_producer!(agent,data,sector,route)

    Y = agent.ext[:sets][:Y]
    S = agent.ext[:sets][:S] 
    b = agent.ext[:variables][:b]
    E = agent.ext[:expressions][:netto_emiss]
    horizon_ets = agent.ext[:parameters][:horizon_ets] = data["horizon_ets"]

    agent.ext[:constraints][:myopic_banking] = @constraint(agent,[y=Y[1:end-data["horizon_ets"]], s=S], sum(b[1:y,s])-sum(E[1:y,s]) <= sum(E[y+1:y+horizon_ets,s]))
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
    CAP_SV = agent.ext[:parameters][:CAP_SV]

    λ_ets = agent.ext[:parameters][:λ_ets]
    λ_product = agent.ext[:parameters][:λ_product]
    agent.ext[:objective] = @objective(agent, Min,
                            sum((r_debt[y]*i[y]*(1-CAP_SV[y])CAPEX[y]*cap[y] + r_equity[y]*λ_ets[y]*b[y] + r_equity[y]*(i[y]*OPEX[y]-λ_product[y])*g[y]) for y in Y, s in S)
                            + sum(r_equity[y]*ρ_ets/2*(b[y]-b_bar[y])^2 for y in Y, s in S)
                            + sum(r_equity[y]*ρ_product/2*(g[y]-g_bar[y])^2 for y in Y, s in S)
                            )
    optimize!(agent)
    return agent
end


function solve_stochastic_producer!(agent::Model)
    @assert !is_myopic(agent) "Agent is myopic"
    @assert is_stochastic(agent) " Agent is not stochastic"

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
                            sum((r_debt[y]*i[y]*CAPEX[y]*cap[y] + r_equity[y]*λ_ets[y,s]*b[y,s] + r_equity[y]*(i[y]*OPEX[y]-λ_product[y,s])*g[y,s]) for y in Y, s in S)
                            + sum(r_equity[y]*ρ_ets/2*(b[y,s]-b_bar[y,s])^2 for y in Y, s in S)
                            + sum(r_equity[y]*ρ_product/2*(g[y,s]-g_bar[y,s])^2 for y in Y, s in S)
                            )
    optimize!(agent)
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
    + sum(r_equity[y]*ρ_ets/2*(b[y]-b_bar[y])^2 for y in Y, s in S)
    + sum(r_equity[y]*ρ_product/2*(g[y]-g_bar[y])^2 for y in Y, s in S)
    + sum(r_equity[y]*ρ_product/20*(g_τ[y,τ]-g_bar_τ[y,τ])^2 for y in Y, τ in Y, s in S)
    )

    optimize!(agent::Model)
    return agent
end