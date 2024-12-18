function build_competitive_fringe!(agent::Model, data::Dict) 
    # Define sets
    build_agent!(agent,data) 

    @assert haskey(data, "nyears") "Data must contain key 'nyears'"

    Y = agent.ext[:sets][:Y]

    # Emissions representative agents, bound to historical values in 2017-2019
    E_ref = agent.ext[:parameters][:E_REF] = data["E_ref"]*ones(data["nyears"],1)
    MAC = agent.ext[:parameters][:MAC]  = data["MAC"]
    agent.ext[:parameters][:mask] = ones(data["nyears"])


    # Define variables
    b = agent.ext[:variables][:b]
    e = agent.ext[:variables][:e] =  @variable(agent, [y=Y], lower_bound=0, base_name="emissions")

    # Define expressions
    bank = agent.ext[:expressions][:bank] = @expression(agent, [y=Y], data["TNAC_2023"] + sum(b[1:y])-sum(e[1:y]))
    agent.ext[:expressions][:netto_emiss] = @expression(agent, [y=Y], e[y] )
 
    # Define constraint
    agent.ext[:constraints][:con1]  = @constraint(agent,[y=Y], bank >= 0)
    #agent.ext[:constraints][:con2] = @constraint(agent,[y=Y], b[y] <= 1.2* data["S"][y])
    #agent.ext[:constraints][:con3]  = @constraint(agent, [y=Y], E_ref[y] - a[y] >= 0)
    # zero production
    g = agent.ext[:variables][:g] = @variable(agent, [y=Y], lower_bound=0, base_name="production") # ton product
    agent.ext[:constraints][:zerogen] = @constraint(agent, [y=Y], g[y] == 0)

    return agent
end

function build_stochastic_competitive_fringe!(agent::Model, data::Dict)
    build_stochastic_agent!(agent,data)
    A = agent.ext[:parameters][:A]
    Y = agent.ext[:sets][:Y]
    S = agent.ext[:sets][:S]

    E_ref = agent.ext[:parameters][:E_REF] = repeat(data["E_ref"]', data["nyears"])
    MAC = agent.ext[:parameters][:MAC]  = data["MAC"]
    agent.ext[:parameters][:mask] = ones(data["nyears"])
    λ_ets = agent.ext[:parameters][:λ_ets] 
    mask = agent.ext[:parameters][:mask]


    # Define variables
    b = agent.ext[:variables][:b]
    e = agent.ext[:variables][:e] =  @variable(agent, [y=Y,s=S], lower_bound=0, base_name="emissions")

    # Define expressions
    bank = agent.ext[:expressions][:bank] = @expression(agent, [y=Y,s=S], data["TNAC_2023"] + sum(b[1:y,s])-sum(e[1:y,s]))
    agent.ext[:expressions][:netto_emiss] = @expression(agent, [y=Y,s=S], e[y,s])
 
    # Define constraint
    agent.ext[:constraints][:con1]  = @constraint(agent,[y=Y,s=S], bank[y,s] >= 0)
   #agent.ext[:constraints][:con3] = @constraint(agent,[y=Y,s=S],E_ref[y,s] - a[y,s] >= 0 )

    # zero production
    g = agent.ext[:variables][:g] = @variable(agent, [y=Y,s=S], lower_bound=0, base_name="production") # ton product
    agent.ext[:constraints][:zerogen] = @constraint(agent, [y=Y,s=S], g[y,s] == 0)

    π_ets = agent.ext[:variables][:π_ets] = @variable(agent, [y=Y,s=S], base_name = "carbon cost")
    agent.ext[:constraints][:π_ets] = @constraint(agent, [y=Y,s=S], π_ets[y,s] <= 0) # Correctly set in solve_risk_averse_fringe!()
    π_MAC = agent.ext[:expressions][:π_MAC] = @expression(agent,[y=Y,s=S], -mask[y]*A[y]*MAC[s]*(E_ref[y,s]-e[y,s])^2)
    
    return agent
end

function build_liquidity_constraint_fringe!(agent::Model,data::Dict)
    build_competitive_fringe!(agent,data)

    agent.ext[:parameters][:isLiquidityConstraint] = true

    bank = agent.ext[:expressions][:bank] 
    Y = agent.ext[:sets][:Y]
    λ_ets = agent.ext[:parameters][:λ_ets] 
 
    agent.ext[:constraints][:liquidity_constraint] = @constraint(agent, [y=Y], bank[y] <= 1000000)

    return agent
end

function build_stochastic_liquidity_constraint_fringe!(agent::Model,data::Dict)
    build_stochastic_competitive_fringe!(agent,data)

    agent.ext[:parameters][:isLiquidityConstraint] = true


    bank = agent.ext[:expressions][:bank] 
    Y = agent.ext[:sets][:Y]
    S = agent.ext[:sets][:S]
    λ_ets = agent.ext[:parameters][:λ_ets] 
 
    agent.ext[:constraints][:liquidity_constraint] = @constraint(agent, [y=Y, s=S], bank[y,s] <= 1000000)

    return agent
end

function build_risk_averse_fringe!(agent::Model, data::Dict)
    build_stochastic_competitive_fringe!(agent,data)

    agent.ext[:parameters][:isRiskAverse] = true

    S = agent.ext[:sets][:S]
    A = agent.ext[:parameters][:A]
    Y = agent.ext[:sets][:Y]
    e = agent.ext[:variables][:e]
    MAC = agent.ext[:parameters][:MAC]
    E_ref = agent.ext[:parameters][:E_REF]
    mask = agent.ext[:parameters][:mask]

    # If agents are risk-averse
    α = agent.ext[:variables][:α] = @variable(agent, base_name ="alpha")
    β = agent.ext[:parameters][:β] = data["cvar"]
    γ = agent.ext[:parameters][:γ] = data["risk_averseness"]
    π_ets = agent.ext[:variables][:π_ets]
    u = agent.ext[:variables_anonymous][:u] = @variable(agent, [s=S], lower_bound=0, base_name = "utility")
    agent.ext[:expressions][:CVAR] = @expression(agent, α - 1/β * sum(u[s] for s in S))
    #agent.ext[:constraints][:π_ets] = @constraint(agent, [y=Y,s=S], π_ets[y,s] <= 0) # Correctly set in solve_risk_averse_fringe!()
    π_MAC = agent.ext[:expressions][:π_MAC] 
    agent.ext[:constraints][:u] = @constraint(agent, [s=S], u[s] >= α - sum(π_ets[y,s]+π_MAC[y,s] for y in Y))
    return agent
end

function solve_competitive_fringe!(agent::Model)
    # update Objective
    A = agent.ext[:parameters][:A]
    Y = agent.ext[:sets][:Y]
    λ_ets = agent.ext[:parameters][:λ_ets] 
    ρ_ets = agent.ext[:parameters][:ρ_ets]
    e = agent.ext[:variables][:e]
    b = agent.ext[:variables][:b]
    b_bar = agent.ext[:parameters][:b_bar]
    MAC = agent.ext[:parameters][:MAC]
    E_ref = agent.ext[:parameters][:E_REF]
    mask = agent.ext[:parameters][:mask]


    agent.ext[:objective] = @objective(agent, Min, 
                                        sum(mask[y]*A[y]*λ_ets[y]*b[y] for y in Y)
                                        + sum(mask[y]*A[y]*MAC*(E_ref[y]-e[y])^2 for y in Y)
                                        + sum(mask[y]*A[y]*ρ_ets/2*(b[y]-b_bar[y])^2 for y in Y)
    )

    # Add liquidity constraint if applicable
    if is_liquidity_constraint(agent)
        for y in Y
            set_normalized_rhs(agent.ext[:constraints][:liquidity_constraint][y], data["TNAC_2023"] * data["P_2023"] / λ_ets[y] - data["TNAC_2023"])
        end
    end
    optimize!(agent)
    return agent
end

function solve_stochastic_competitive_fringe!(agent::Model)
    @assert is_stochastic(agent) " Agent is not stochastic"

    # update Objective
    A = agent.ext[:parameters][:A]
    S = agent.ext[:sets][:S]
    Y = agent.ext[:sets][:Y]
    λ_ets = agent.ext[:parameters][:λ_ets] 
    ρ_ets = agent.ext[:parameters][:ρ_ets]
    b = agent.ext[:variables][:b]
    b_bar = agent.ext[:parameters][:b_bar]
    MAC = agent.ext[:parameters][:MAC]
    e = agent.ext[:variables][:e]
    E_ref = agent.ext[:parameters][:E_REF]
    mask = agent.ext[:parameters][:mask]


    if is_risk_averse(agent)
        α = agent.ext[:variables][:α]
        β = agent.ext[:parameters][:β]
        γ = agent.ext[:parameters][:γ]
        u = agent.ext[:variables_anonymous][:u] 
        CVAR = agent.ext[:expressions][:CVAR]
        π_ets = agent.ext[:variables][:π_ets]
        π_MAC = agent.ext[:expressions][:π_MAC]
        #for y in Y, s in S
            delete.(agent,agent.ext[:constraints][:π_ets])
            #set_normalized_coefficient(agent.ext[:constraints][:π_ets][y,s], b[y,s], mask[y]*A[y]*λ_ets[y,s])
        #end
        agent.ext[:constraints][:π_ets] = @constraint(agent,[y=Y,s=S], π_ets[y,s] + mask[y]*A[y]*λ_ets[y,s]*b[y,s] <= 0 )

        # TO DO: for the rolling horizon approach the window changes and possibly I do need to update this expression
        #agent.ext[:constraints][:u] = @constraint(agent, [s=S], u[s] >= α + sum(π_ets[y,s]+π_MAC[y,s] for y in Y))

        agent.ext[:objective] = @objective(agent, Max,  γ * sum(π_ets[y,s] + π_MAC[y,s] for y in Y, s in S) # Profit
                                                     + (1-γ)*(α - (1/β) * sum(u[s] for s in S)) # CVaR
                                                      - sum(mask[y]*A[y]*ρ_ets/2*(b[y,s]-b_bar[y,s])^2 for y in Y, s in S) # Penalty term
        )
        # TO DO: for the rolling horizon approach the window changes and possibly I do need to update this constraint
        #agent.ext[:constraints][:u] = @constraint(agent, [s=S], u[s] >= α + sum(π_ets[y,s]+π_MAC[y,s] for y in Y))
    else
        π_ets = agent.ext[:variables][:π_ets]
        π_MAC = agent.ext[:expressions][:π_MAC]
        for y in Y, s in S
            delete.(agent,agent.ext[:constraints][:π_ets][y,s])
            #set_normalized_coefficient(agent.ext[:constraints][:π_ets][y,s], b[y,s], mask[y]*A[y]*λ_ets[y,s])
        end
        agent.ext[:constraints][:π_ets] = @constraint(agent,[y=Y,s=S], π_ets[y,s] + mask[y]*A[y]*λ_ets[y,s]*b[y,s] <= 0 )
        agent.ext[:objective] = @objective(agent, Max, 
                                                sum(π_ets[y,s] + π_MAC[y,s] for y in Y, s in S) 
                                        #    sum(mask[y]*A[y]*λ_ets[y,s]*b[y,s] for y in Y, s in S)
                                        #    + sum(mask[y]*A[y]*MAC[s]*(E_ref[y,s]-e[y,s])^2 for y in Y, s in S)
                                            - sum(mask[y]*A[y]*ρ_ets/2*(b[y,s]-b_bar[y,s])^2 for y in Y, s in S)
        )
    end

    if is_liquidity_constraint(agent)
        for y in Y, s in S
            set_normalized_rhs(agent.ext[:constraints][:liquidity_constraint][y,s], data["TNAC_2023"] * data["P_2023"] / λ_ets[y, s] - data["TNAC_2023"])
        end
    end

    optimize!(agent)
    return agent
end

