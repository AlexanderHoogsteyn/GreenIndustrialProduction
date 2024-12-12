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
    #agent.ext[:constraints][:buycons] = @constraint(agent,[y=Y], sum(b[1:y]) >= sum(g[1:y]*EF))

    # Prohibit banking:
    agent.ext[:constraints][:buycons] = @constraint(agent,[y=Y], b[y] >= g[y]*EF)


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

    # Prohibit banking:
    agent.ext[:constraints][:buycons] = @constraint(agent,[y=Y,s=S], b[y,s] >= g[y,s]*EF)

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

function build_liquidity_constraint_producer!(agent::Model,data::Dict,sector::String,route::String)
    build_producer!(agent,data,sector,route)

    agent.ext[:parameters][:isLiquidityConstraint] = true

    Y = agent.ext[:sets][:Y]
    EF = agent.ext[:parameters][:EF]
    g = agent.ext[:variables][:g] 
    b = agent.ext[:variables][:b] 

    agent.ext[:constraints][:nobanking] = @constraint(agent,[y=Y], b[y] >= g[y]*EF)

    return agent 
end

function build_stochastic_liquidity_constraint_producer!(agent::Model,data::Dict,sector::String,route::String)
    build_stochastic_producer!(agent,data,sector,route)

    agent.ext[:parameters][:isLiquidityConstraint] = true

    S = agent.ext[:sets][:S]   
    Y = agent.ext[:sets][:Y]
    EF = agent.ext[:parameters][:EF]
    g = agent.ext[:variables][:g] 
    b = agent.ext[:variables][:b] 

    agent.ext[:constraints][:nobanking] = @constraint(agent,[y=Y,s=S], b[y,s] >= g[y,s]*EF)

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
                            )


    if is_liquidity_constraint(agent)
        if haskey(agent.ext[:constraints], :liquidity_constraint)
            delete.(agent,agent.ext[:constraints][:liquidity_constraint])
        end 
        #agent.ext[:constraints][:liquidity_constraint] = @constraint(
        #    agent, [y=Y],
        #     sum(b[i] - g[i]*EF for i in 1:y) <= 0.0* data["TNAC_2023"] 
        #)
    end
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

    if is_liquidity_constraint(agent)
        if haskey(agent.ext[:constraints], :liquidity_constraint)
             delete.(agent,agent.ext[:constraints][:liquidity_constraint])
         end 
        # agent.ext[:constraints][:liquidity_constraint] = @constraint(
        #     agent, [y=Y,s=S],
        #     (data["TNAC_2023"] + sum(b[i,s] - g[i,s]*EF for i in 1:y)) * λ_ets[y,s] <= data["TNAC_2023"] * data["P_2023"]
        #     )
    end

    optimize!(agent)
    return agent
end
