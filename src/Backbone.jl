function needed_eua_price(hydrogen_price,benchmark_route::Dict,route::Dict)
    return 1000*hydrogen_price*route["Needs"]["Hydrogen"] / (benchmark_route["Needs"]["ETS"]-route["Needs"]["ETS"] )
end

function route_costs(commodityPrices::Dict{Any,Any}, route::Dict{Any,Any},policies::Dict{Any,Any})
    cost = 0;
    needs = route["OPEX"]

    # Calculate route costs
    for (need, amount) in needs
        cost += commodityPrices[need]*amount;
    end

    # Account for cost reducing policies
    if "floor" in keys(policies)
        cost = cost - floor_revenue(policies["floor"], commodityPrices, route)
    end
    if "CCfD" in keys(policies)
        cost = cost - CCfD_revenue(policies["CCfD"], commodityPrices, route)
    end
    if "Grandfathering" in keys(policies)
        cost = cost - policies["Grandfathering"]["Benchmark"]*commodityPrices["ETS"]
    end
    return cost
end

function route_costs(commodityPrices::Dict{Any,Any}, route::Dict{Any,Any})
    cost = route_costs(commodityPrices, route, Dict())
    return cost
end

function cost_component(commodityPrices::Dict{String,Float64}, route::Dict{Any,Any},component)
    if component in keys(route["Needs"])
        return commodityPrices[component]*route["Needs"][component]
    else 
        return 0.0
    end
end

function is_myopic(mod::Model)
    if haskey(mod.ext[:parameters], :isMyopic) && mod.ext[:parameters][:isMyopic] == true
        return true
    else 
        return false
    end
end

function is_stochastic(mod::Model)
    if haskey(mod.ext[:parameters], :isStochastic) && mod.ext[:parameters][:isStochastic] == true
        return true
    else 
        return false
    end
end

function is_rolling_horizon(ADMM::Dict)
    if haskey(ADMM, :isRollingHorizon) && ADMM[:isRollingHorizon] == true
        return true
    else 
        return false
    end
end

function set_lookahead_window!(agent::Model,ADMM::Dict)
    # Constraints an agents decision variables outside the lookahead to what they currently are
    agent.ext[:constraints_rolling_horizon] = Dict()
    Y = agent.ext[:sets][:Y][ADMM[:end]+1:end]

    for variable in agent.ext[:variables]
        agent.ext[:constraints_rolling_horizon][variable] = 
            @constraint(agent, [y=Y], agent.ext[:variables][variable][y] == 0) 
    end
    return agent
end 

function set_lookahead_window!(agents::Dict,ADMM::Dict)
    for agent in agents 
        set_lookahead_window!(agent,ADMM)
    end
    return agents
end

function move_lookahead_window!(agent::Model,ADMM::Dict)
    # This releases the constraints on a model's decision variables on the next year, e.g. moves the lookahead window one further.
    # It fixes the decsion variable of the start of the window to its current vlue
    ADMM[:end] += 1
    for variable in agent.ext[:variables]
        delete.(agent,agent.ext[:constraints_rolling_horizon][variable][ADMM[:end]])
        values = JuMP.value.(agent.ext[:variables][variable][ADMM[:start]])
        @constraint(agent, agent.ext[:constraints_rolling_horizon][variable][ADMM[:start]] == values)
    end
    ADMM[:start] += 1
    mask = zeros(size(ADMM["Imbalances"]["ETS"][end]))
    mask[ADMM[:start]:ADMM[:end]] .= 1
    ADMM[:mask] = mask
        return agent 
end

function move_lookahead_window!(agents::Dict,ADMM::Dict)
    for agent in agents 
        move_lookahead_window!(agent,ADMM)
    end
    return agents
end