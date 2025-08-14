function build_financial_agent!(agent::Model, data::Dict) 
    # Define sets
    build_agent!(agent,data) 
    Y = agent.ext[:sets][:Y]

    # Define variables
    b = agent.ext[:variables][:b]
    g = agent.ext[:variables][:g] = @variable(agent, [y=Y], lower_bound=0, base_name="production") # ton product

    # Define expressions
    agent.ext[:expressions][:netto_emiss] = @expression(agent, [y=Y], 0)
    bank = agent.ext[:expressions][:bank] = @expression(agent, [y=Y], data["TNAC_2023"] + sum(b[1:y]))

    # Define constraint
    agent.ext[:constraints][:con1]  = @constraint(agent,[y=Y], bank[y] >= 0)
    agent.ext[:constraints][:zerogen] = @constraint(agent, [y=Y], g[y] == 0)

    return agent
end

function build_stochastic_financial_agent!(agent::Model, data::Dict)
    build_stochastic_agent!(agent,data)
    Y = agent.ext[:sets][:Y]
    S = agent.ext[:sets][:S]

    # Define variables
    b = agent.ext[:variables][:b]
    g = agent.ext[:variables][:g] = @variable(agent, [y=Y,s=S], lower_bound=0, base_name="production") # ton product

    # Define expressions
    bank = agent.ext[:expressions][:bank] = @expression(agent, [y=Y,s=S], data["TNAC_2023"] + sum(b[1:y,s]))
    agent.ext[:expressions][:netto_emiss] = @expression(agent, [y=Y,s=S], 0)
 
    # Define constraint
    agent.ext[:constraints][:buycons_supply] = @constraint(agent,[y=Y,s=S], b[y,s] <= data["CAP"][y])
    agent.ext[:constraints][:con1]  = @constraint(agent,[y=Y,s=S], bank[y,s] >= 0)

    # zero production
    agent.ext[:constraints][:zerogen] = @constraint(agent, [y=Y,s=S], g[y,s] == 0)

    return agent
end

function build_liquidity_constraint_financial_agent!(agent::Model,data::Dict)
    build_competitive_financial_agent!(agent,data)

    agent.ext[:parameters][:isLiquidityConstraint] = true

    bank = agent.ext[:expressions][:bank]
    b = agent.ext[:variables][:b]
    Y = agent.ext[:sets][:Y]
    λ_ets = agent.ext[:parameters][:λ_ets]
 
    # Prevent out of bounds error in first iteration
    agent.ext[:constraints][:liquidity_constraint] = @constraint(agent, [y=Y],  sum(b[1:y]) <= 1000000)

    return agent
end

function build_stochastic_liquidity_constraint_financial_agent!(agent::Model,data::Dict)
    build_stochastic_financial_agent!(agent,data)

    agent.ext[:parameters][:isLiquidityConstraint] = true


    bank = agent.ext[:expressions][:bank] 
    b = agent.ext[:variables][:b]
    Y = agent.ext[:sets][:Y]
    S = agent.ext[:sets][:S]
    λ_ets = agent.ext[:parameters][:λ_ets] 
 
    # Prevent out of bounds error in first iteration
    agent.ext[:constraints][:liquidity_constraint] = @constraint(agent, [y=Y, s=S], sum(b[1:y,s]) <= 1000000)

    return agent
end

function solve_stochastic_financial_agent!(agent::Model,data::Dict)
        # update Objective
        A = agent.ext[:parameters][:A]
        S = agent.ext[:sets][:S]
        Y = agent.ext[:sets][:Y]
        λ_ets = agent.ext[:parameters][:λ_ets] 
        ρ_ets = agent.ext[:parameters][:ρ_ets]
        b = agent.ext[:variables][:b]
        b_bar = agent.ext[:parameters][:b_bar]
        r_equity = agent.ext[:parameters][:r_equity]
        mask = agent.ext[:parameters][:mask]
        start = agent.ext[:parameters][:start] 
        end_ = agent.ext[:parameters][:end]  

        agent.ext[:objective] = @objective(agent, Min, 
                                            sum(mask[y]*r_equity[y]*λ_ets[y,s]*b[y,s] for y in Y, s in S)
                                            + sum(mask[y]*r_equity[y]*ρ_ets/2*(b[y,s]-b_bar[y,s])^2 for y in Y, s in S)
        )

        agent.ext[:expressions][:π] = @expression(agent, [y=Y, s=S], 
                                                    r_equity[y]*λ_ets[y,s]*b[y,s]
        )
        agent.ext[:expressions][:cumprofit] = @expression(agent, [yy=Y, s=S], 
                                                    sum(r_equity[y]*λ_ets[y,s]*b[y,s] for y in 1:yy)
        )

        if is_liquidity_constraint(agent)
            for y in Y[start:end_], s in S
                set_normalized_rhs(agent.ext[:constraints][:liquidity_constraint][y,s], (data["TNAC_2023"] * data["P_2024"] * data["liquidity_factor"]) / λ_ets[y, s] - data["TNAC_2023"])
            end
        end
 

        optimize!(agent)
    return agent
end