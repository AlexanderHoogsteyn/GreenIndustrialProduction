function build_trader!(agent::Model, data::Dict) 
    # Define sets
    build_agent!(agent,data) 
    Y = agent.ext[:sets][:Y]

    # Define variables
    b = agent.ext[:variables][:b]

    # Define expressions
    agent.ext[:expressions][:bank] = @expression(agent, [y=Y], sum(b[1:y]))
    agent.ext[:expressions][:netto_emiss] = @expression(agent, [y=Y], 0)
 
    # Define constraint
    agent.ext[:constraints][:buycons] = @constraint(agent, sum(b) >= 0)
    agent.ext[:constraints][:buycons_supply] = @constraint(agent,[y=Y], b[y] <= data["CAP"][y])

    # zero production
    g = agent.ext[:variables][:g] = @variable(agent, [y=Y], lower_bound=0, base_name="production") # ton product
    agent.ext[:constraints][:zerogen] = @constraint(agent, [y=Y], g[y] == 0)
    return agent
end

function build_stochastic_trader!(agent::Model, data::Dict)
    build_stochastic_agent!(agent,data)
    Y = agent.ext[:sets][:Y]
    S = agent.ext[:sets][:S]

    # Define variables
    b = agent.ext[:variables][:b]

    # Define expressions
    agent.ext[:expressions][:bank] = @expression(agent, [y=Y,s=S], sum(b[1:y,s]))
    agent.ext[:expressions][:netto_emiss] = @expression(agent, [y=Y,s=S], 0)
 
    # Define constraint
    agent.ext[:constraints][:buycons] = @constraint(agent,[s=S], sum(b[:,s]) >= 0)
    agent.ext[:constraints][:buycons_supply] = @constraint(agent,[y=Y,s=S], b[y,s] <= data["CAP"][y])

    # zero production
    g = agent.ext[:variables][:g] = @variable(agent, [y=Y,s=S], lower_bound=0, base_name="production") # ton product
    agent.ext[:constraints][:zerogen] = @constraint(agent, [y=Y,s=S], g[y,s] == 0)
    return agent
end

function solve_trader!(agent::Model)
    # update Objective
    A = agent.ext[:parameters][:A]
    Y = agent.ext[:sets][:Y]
    λ_ets = agent.ext[:parameters][:λ_ets] 
    ρ_ets = agent.ext[:parameters][:ρ_ets]
    b = agent.ext[:variables][:b]
    b_bar = agent.ext[:parameters][:b_bar]
    r_equity = agent.ext[:parameters][:r_equity]


    agent.ext[:objective] = @objective(agent, Min, 
                                        sum(A[y]*λ_ets[y]*b[y] for y in Y)
                                        + sum(A[y]*ρ_ets/2*(b[y]-b_bar[y])^2 for y in Y)
    )

    optimize!(agent)
    return agent
end

function solve_stochastic_trader!(agent,data)
        # update Objective
        A = agent.ext[:parameters][:A]
        S = agent.ext[:sets][:S]
        Y = agent.ext[:sets][:Y]
        λ_ets = agent.ext[:parameters][:λ_ets] 
        ρ_ets = agent.ext[:parameters][:ρ_ets]
        b = agent.ext[:variables][:b]
        b_bar = agent.ext[:parameters][:b_bar]
        r_equity = agent.ext[:parameters][:r_equity]
    
        agent.ext[:objective] = @objective(agent, Min, 
                                            sum(A[y]*λ_ets[y,s]*b[y,s] for y in Y, s in S)
                                            + sum(A[y]*ρ_ets/2*(b[y,s]-b_bar[y,s])^2 for y in Y, s in S)
        )
        optimize!(agent)
    return agent
end