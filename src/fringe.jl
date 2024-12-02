function build_competitive_fringe!(agent::Model, data::Dict) 
    # Define sets
    build_agent!(agent,data) 

    @assert haskey(data, "nyears") "Data must contain key 'nyears'"

    Y = agent.ext[:sets][:Y]

    # Emissions representative agents, bound to historical values in 2017-2019
    E_ref = agent.ext[:parameters][:E_REF] = data["E_ref"]*ones(data["nyears"],1)
    MAC = agent.ext[:parameters][:MAC]  = data["MAC"]

    # Define variables
    b = agent.ext[:variables][:b]
    a = agent.ext[:variables][:a] =  @variable(agent, [y=Y], lower_bound=0, base_name="abatement")

    # Define expressions
    bank = agent.ext[:expressions][:bank] = @expression(agent, [y=Y], data["TNAC_2023"] + sum(b[1:y])-sum(E_ref[1:y]) + sum(a[1:y]))
    agent.ext[:expressions][:netto_emiss] = @expression(agent, [y=Y], E_ref[y] - a[y])
 
    # Define constraint
    agent.ext[:constraints][:con1]  = @constraint(agent,[y=Y], bank >= 0)
    #agent.ext[:constraints][:con2] = @constraint(agent,[y=Y], b[y] <= 1.2* data["S"][y])
    agent.ext[:constraints][:con3]  = @constraint(agent, [y=Y], E_ref[y] - a[y] >= 0)
    # zero production
    g = agent.ext[:variables][:g] = @variable(agent, [y=Y], lower_bound=0, base_name="production") # ton product
    agent.ext[:constraints][:zerogen] = @constraint(agent, [y=Y], g[y] == 0)
    return agent
end

function build_stochastic_competitive_fringe!(agent::Model, data::Dict)
    build_stochastic_agent!(agent,data)
    Y = agent.ext[:sets][:Y]
    S = agent.ext[:sets][:S]

    E_ref = agent.ext[:parameters][:E_REF] = repeat(data["E_ref"]', data["nyears"])
    MAC = agent.ext[:parameters][:MAC]  = data["MAC"]


    # Define variables
    b = agent.ext[:variables][:b]
    a = agent.ext[:variables][:a] =  @variable(agent, [y=Y,s=S], lower_bound=0, base_name="abatement")

    # Define expressions
    bank = agent.ext[:expressions][:bank] = @expression(agent, [y=Y,s=S], data["TNAC_2023"] + sum(b[1:y,s])-sum(E_ref[1:y,s]) + sum(a[1:y,s]))
    agent.ext[:expressions][:netto_emiss] = @expression(agent, [y=Y,s=S], E_ref[y,s] - a[y,s])
 
    # Define constraint
    agent.ext[:constraints][:con1]  = @constraint(agent,[y=Y,s=S], bank[y,s] >= 0)
    agent.ext[:constraints][:con3] = @constraint(agent,[y=Y,s=S],E_ref[y,s] - a[y,s] >= 0 )

    # zero production
    g = agent.ext[:variables][:g] = @variable(agent, [y=Y,s=S], lower_bound=0, base_name="production") # ton product
    agent.ext[:constraints][:zerogen] = @constraint(agent, [y=Y,s=S], g[y,s] == 0)
    return agent
end

function build_liquidity_constraint_fringe!(agent::Model,data::Dict)
    build_competitive_fringe!(agent,data)

    agent.ext[:parameters][:isLiquidityConstraint] = true

    return agent
end

function build_stochastic_liquidity_constraint_fringe!(agent::Model,data::Dict)
    build_stochastic_competitive_fringe!(agent,data)

    agent.ext[:parameters][:isLiquidityConstraint] = true

    return agent
end

function solve_competitive_fringe!(agent::Model)
    # update Objective
    A = agent.ext[:parameters][:A]
    Y = agent.ext[:sets][:Y]
    λ_ets = agent.ext[:parameters][:λ_ets] 
    ρ_ets = agent.ext[:parameters][:ρ_ets]
    a = agent.ext[:variables][:a]
    b = agent.ext[:variables][:b]
    b_bar = agent.ext[:parameters][:b_bar]
    bank = agent.ext[:expressions][:bank]
    MAC = agent.ext[:parameters][:MAC]

    agent.ext[:objective] = @objective(agent, Min, 
                                        sum(A[y]*λ_ets[y]*b[y] for y in Y)
                                        + sum(A[y]*MAC*(a[y])^2 for y in Y)
                                        + sum(A[y]*ρ_ets/2*(b[y]-b_bar[y])^2 for y in Y)
                                        #+ sum(A[y]*ρ_ets/2*(a[y]-a_bar[y])^2 for y in Y)
    )

    # Add liquidity constraint if applicable
    if is_liquidity_constraint(agent)
        if haskey(agent.ext[:constraints], :liquidity_constraint)
            delete.(agent,agent.ext[:constraints][:liquidity_constraint])
        end 
        agent.ext[:constraints][:liquidity_constraint] = @constraint(
            agent, [y=Y],
            bank[y] * λ_ets[y] <= data["TNAC_2023"] * data["P_2023"]
        )
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
    bank =  agent.ext[:expressions][:bank]
    MAC = agent.ext[:parameters][:MAC]
    a = agent.ext[:variables][:a]


    agent.ext[:objective] = @objective(agent, Min, 
                                        sum(A[y]*λ_ets[y,s]*b[y,s] for y in Y, s in S)
                                        + sum(A[y]*MAC[s]*(a[y,s])^2 for y in Y, s in S)
                                        + sum(A[y]*ρ_ets/2*(b[y,s]-b_bar[y,s])^2 for y in Y, s in S)
    )

    if is_liquidity_constraint(agent)
        if haskey(agent.ext[:constraints], :liquidity_constraint)
            delete.(agent,agent.ext[:constraints][:liquidity_constraint])
        end 
        agent.ext[:constraints][:liquidity_constraint] = @constraint(agent, [y=Y, s=S], bank[y,s] * λ_ets[y, s] <= data["TNAC_2023"] * data["P_2023"])
    end

    optimize!(agent)
    return agent
end

