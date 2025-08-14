function build_competitive_compliance_agent!(agent::Model, data::Dict) 
    """
    Builds a competitive compliance_agent agent model with specific parameters and variables.
    Arguments:
    - `agent::Model`: Parse an empty JuMP model that will be mutatted to contain the agent.
    - `data::Dict`: Dictionary that stores agent data.
    """

    # Define sets
    build_agent!(agent,data) 

    @assert haskey(data, "nyears") "Data must contain key 'nyears'"

    Y = agent.ext[:sets][:Y]

    # Emissions representative agents, bound to historical values in 2017-2019
    E_ref = agent.ext[:parameters][:E_REF] = data["E_ref"]*ones(data["nyears"],1)
    MAC = agent.ext[:parameters][:MAC]  = data["MAC"]


    # Define variables
    b = agent.ext[:variables][:b]
    e = agent.ext[:variables][:e] = @variable(agent, [y=Y], lower_bound=0, base_name="emissions")
    g = agent.ext[:variables][:g] = @variable(agent, [y=Y], lower_bound=0, base_name="production") # Will be imposed zero


    # Define expressions
    bank = agent.ext[:expressions][:bank] = @expression(agent, [y=Y], sum(b[1:y])-sum(e[1:y]))
    agent.ext[:expressions][:netto_emiss] = @expression(agent, [y=Y], e[y] )
 
    # Define constraints
        # zero production
    agent.ext[:constraints][:zerogen] = @constraint(agent, [y=Y], g[y] == 0)

        # Zero banking
    agent.ext[:constraints][:zerobank] = @constraint(agent, [y=Y], b[y] >= e[y])

    return agent
end

function build_stochastic_competitive_compliance_agent!(agent::Model, data::Dict)
    """
    Builds a stochastic variant of the competitive compliance_agent agent model with specific parameters and variables.
    Arguments:
    - `agent::Model`: Parse an empty JuMP model that will be mutatted to contain the agent.
    - `data::Dict`: Dictionary that stores agent data.
    """
    # Define sets
    build_stochastic_agent!(agent,data)
    Y = agent.ext[:sets][:Y]
    S = agent.ext[:sets][:S]

    A = agent.ext[:parameters][:A]
    E_ref = agent.ext[:parameters][:E_REF] = repeat(data["E_ref"]', data["nyears"])
    MAC = agent.ext[:parameters][:MAC]  = data["MAC"]
    mask = agent.ext[:parameters][:mask] = ones(data["nyears"])

    # Define variables
    b = agent.ext[:variables][:b]
    e = agent.ext[:variables][:e] =  @variable(agent, [y=Y,s=S], lower_bound=0, base_name="emissions")
    g = agent.ext[:variables][:g] = @variable(agent, [y=Y,s=S], lower_bound=0, base_name="production") # Will be imposed zero


    # Define expressions
    bank = agent.ext[:expressions][:bank] = @expression(agent, [y=Y,s=S], sum(b[1:y,s])-sum(e[1:y,s]))
    agent.ext[:expressions][:netto_emiss] = @expression(agent, [y=Y,s=S], e[y,s])
    agent.ext[:expressions][:π_MAC] = @expression(agent,[y=Y,s=S], mask[y]*A[y]*MAC[s]*(E_ref[y,s]-e[y,s])^2)
    agent.ext[:expressions][:π] = @expression(agent,[y=Y,s=S], A[y]*MAC[s]*(E_ref[y,s]-e[y,s])^2)
    agent.ext[:expressions][:cost] = @expression(agent,[y=Y,s=S], A[y]*MAC[s]*(E_ref[y,s]-e[y,s])^2)

 
    # Define constraints
        # zero production
    agent.ext[:constraints][:zerogen] = @constraint(agent, [y=Y,s=S], g[y,s] == 0)
        # Zero banking
    agent.ext[:constraints][:zerobank] = @constraint(agent, [y=Y,s=S], b[y,s] >= e[y,s])

    return agent
end

function solve_competitive_compliance_agent!(agent::Model)
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

    agent.ext[:expressions][:cumprofit] = @expression(agent,[yy=Y],  
                                            sum(A[y]*λ_ets[y]*b[y]+ A[y]*MAC*(E_ref[y]-e[y])^2 for y in 1:yy)
    )

    # Add liquidity constraint if applicable
    if is_liquidity_constraint(agent)
        for y in Y
            set_normalized_rhs(agent.ext[:constraints][:liquidity_constraint][y], data["TNAC_2023"] * data["P_2024"] / λ_ets[y] - data["TNAC_2023"])
        end
    end
    optimize!(agent)
    return agent
end

function solve_stochastic_competitive_compliance_agent!(agent::Model,data::Dict)
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



    agent.ext[:objective] = @objective(agent, Min, 
                                        sum(mask[y]*A[y]*λ_ets[y,s]*b[y,s] for y in Y, s in S)
                                        + sum(mask[y]*A[y]*MAC[s]*(E_ref[y,s]-e[y,s])^2 for y in Y, s in S)
                                        + sum(mask[y]*A[y]*ρ_ets/2*(b[y,s]-b_bar[y,s])^2 for y in Y, s in S)
    )
    if is_liquidity_constraint(agent)
        for y in Y, s in S
            set_normalized_rhs(agent.ext[:constraints][:liquidity_constraint][y,s], data["TNAC_2023"] * data["P_2024"] * data["liquidity_factor"] / λ_ets[y,s] - data["TNAC_2023"])
        end
    end

    optimize!(agent)
    return agent
end

