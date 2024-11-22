function build_competitive_fringe!(agent::Model, data::Dict) 
    # Define sets
    build_agent!(agent,data) 

    @assert haskey(data, "nyears") "Data must contain key 'nyears'"

    Y = agent.ext[:sets][:Y]

    # Emissions representative agents, bound to historical values in 2017-2019
    E = agent.ext[:parameters][:e] = zeros(data["nyears"],1)
    E[1] = data["E_2019"]

    # Define variables
    b = agent.ext[:variables][:b]

    # Define expressions
    agent.ext[:expressions][:bank] = @expression(agent, [y=Y], data["TNAC_2023"] + sum(b[1:y])-sum(E[1:y]))
    agent.ext[:expressions][:netto_emiss] = @expression(agent, [y=Y], E[y])
 
    # Define constraint
    agent.ext[:constraints][:con1]  = @constraint(agent,[y=Y], data["TNAC_2023"] + sum(b[1:y]) >= sum(E[1:y]))
    #agent.ext[:constraints][:con2] = @constraint(agent,[y=Y], b[y] <= 1.2* data["S"][y])
    # zero production
    g = agent.ext[:variables][:g] = @variable(agent, [y=Y], lower_bound=0, base_name="production") # ton product
    agent.ext[:constraints][:zerogen] = @constraint(agent, [y=Y], g[y] == 0)
    return agent
end

function build_stochastic_competitive_fringe!(agent::Model, data::Dict)
    build_stochastic_agent!(agent,data)
    Y = agent.ext[:sets][:Y]
    S = agent.ext[:sets][:S]

    # Emissions representative agents, bound to historical values in 2017-2019
    E = agent.ext[:parameters][:e] = zeros(data["nyears"],data["nsamples"])
    E[1,:] .= data["E_2019"]

    # Define variables
    b = agent.ext[:variables][:b]

    # Define expressions
    agent.ext[:expressions][:bank] = @expression(agent, [y=Y,s=S], data["TNAC_2023"] + sum(b[1:y,s])-sum(E[1:y,s]))
    agent.ext[:expressions][:netto_emiss] = @expression(agent, [y=Y,s=S], E[y,s])
 
    # Define constraint
    agent.ext[:constraints][:con1]  = @constraint(agent,[y=Y,s=S], data["TNAC_2023"] + sum(b[1:y,s]) >= sum(E[1:y,s]))

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
    E = agent.ext[:parameters][:e]
    b = agent.ext[:variables][:b]
    b_bar = agent.ext[:parameters][:b_bar]

    agent.ext[:objective] = @objective(agent, Min, 
                                        sum(A[y]*λ_ets[y]*b[y] for y in Y)
                                        + sum(A[y]*ρ_ets/2*(b[y]-b_bar[y])^2 for y in Y)
    )

    # Check if constraints exist and delete them if they do
    if haskey(agent.ext[:constraints], :con1)
        delete.(agent, agent.ext[:constraints][:con1])
    end

    # Add the new constraints
    agent.ext[:constraints][:con1] = @constraint(agent, [y=Y], data["TNAC_2023"] + sum(b[1:y]) >= sum(E[1:y]))

    agent.ext[:expressions][:bank] = @expression(agent, [y=Y],data["TNAC_2023"] + sum(b[1:y])-sum(E[1:y]))
    agent.ext[:expressions][:netto_emiss] = @expression(agent, [y=Y], E[y])

    # Add liquidity constraint if applicable
    if is_liquidity_constraint(agent)
        if haskey(agent.ext[:constraints], :liquidity_constraint)
            delete.(agent,agent.ext[:constraints][:liquidity_constraint])
        end 
        agent.ext[:constraints][:liquidity_constraint] = @constraint(
            agent, [y=Y],
            (data["TNAC_2023"] + sum(b[i] - E[i] for i in 1:y)) * λ_ets[y] <= data["TNAC_2023"] * data["P_2023"]
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
    E = agent.ext[:parameters][:e]
    b = agent.ext[:variables][:b]
    b_bar = agent.ext[:parameters][:b_bar]
    bank =  agent.ext[:expressions][:bank]

    agent.ext[:objective] = @objective(agent, Min, 
                                        sum(A[y]*λ_ets[y,s]*b[y,s] for y in Y, s in S)
                                        + sum(A[y]*ρ_ets/2*(b[y,s]-b_bar[y,s])^2 for y in Y, s in S)
    )

    # Check if constraints exist and delete them if they do
    if haskey(agent.ext[:constraints], :con1)
        delete.(agent, agent.ext[:constraints][:con1])
    end

    # Add the new constraints
    agent.ext[:constraints][:con1] = @constraint(agent, [y=Y, s=S], data["TNAC_2023"] + sum(b[1:y,s]) >= sum(E[1:y,s]))

    bank = agent.ext[:expressions][:bank] = @expression(agent, [y=Y, s=S], data["TNAC_2023"] + sum(b[1:y,s])-sum(E[1:y,s]))
    agent.ext[:expressions][:netto_emiss] = @expression(agent, [y=Y,s=S], E[y,s])

    if is_liquidity_constraint(agent)
        if haskey(agent.ext[:constraints], :liquidity_constraint)
            delete.(agent,agent.ext[:constraints][:liquidity_constraint])
        end 
        agent.ext[:constraints][:liquidity_constraint] = @constraint(agent, [y=Y, s=S], bank[y,s] * λ_ets[y, s] <= data["TNAC_2023"] * data["P_2023"])
    end

    optimize!(agent)
    return agent
end

