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