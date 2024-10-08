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

function is_stochastic(variable::JuMP.Containers.DenseAxisArray)
    if ndims(variable) > 1 
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
    Y = agent.ext[:sets][:Y]
    Y_window = Y[ADMM[:end]+1:end]

    for (variable_name,variable) in agent.ext[:variables]
        agent.ext[:constraints_rolling_horizon][variable_name] = Dict()  # Initialize a Dict to hold constraints
        agent.ext[:constraints_rolling_horizon][variable_name] =  JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(undef, Y)
        for y in Y_window
            agent.ext[:constraints_rolling_horizon][variable_name][y] = 
            @constraint(agent, agent.ext[:variables][variable_name][y] == 0) 
        end
    end
    return agent
end 

function set_lookahead_window!(agents::Dict,ADMM::Dict)
        for (agent,model) in agents
            if is_stochastic(model)
                set_lookahead_window_stochastic!(model,ADMM)
            else
                set_lookahead_window!(model,ADMM)
            end
        end
    return agents
end

function set_lookahead_window_stochastic!(agent::Model,ADMM::Dict)
    agent.ext[:constraints_rolling_horizon] = Dict()
    Y = agent.ext[:sets][:Y]
    Y_window = Y[ADMM[:end]+1:end]
    S = agent.ext[:sets][:S]   

    for (variable_name,variable) in agent.ext[:variables]
        agent.ext[:constraints_rolling_horizon][variable_name] = Dict()  # Initialize a Dict to hold constraints
        if is_stochastic(variable)
            agent.ext[:constraints_rolling_horizon][variable_name] = JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(undef,Y,S)
            for y in Y_window, s in S
                agent.ext[:constraints_rolling_horizon][variable_name][y,s] = 
                @constraint(agent, agent.ext[:variables][variable_name][y,s] == 0) 
            end
        else
            agent.ext[:constraints_rolling_horizon][variable_name] =  JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(undef, Y)
            for y in Y_window
                agent.ext[:constraints_rolling_horizon][variable_name][y] = 
                @constraint(agent, agent.ext[:variables][variable_name][y] == 0) 
            end
        end
    end
end

function move_lookahead_window!(agent::Model,ADMM::Dict)
    # This releases the constraints on a model's decision variables on the next year, e.g. moves the lookahead window one further.
    # It fixes the decsion variable of the start of the window to its current vlue
    for (variable_name, variable) in agent.ext[:variables]
        optimize!(agent) # TO DO check the OptimizenotCalled error
        old_value = JuMP.value.(agent.ext[:variables][variable_name][ADMM[:start]])
        try
            delete.(agent,agent.ext[:constraints_rolling_horizon][variable_name][ADMM[:end]])
        catch error
            print("Reached end of window")
        end
        agent.ext[:constraints_rolling_horizon][variable_name][ADMM[:start]] = @constraint(agent, agent.ext[:variables][variable_name][ADMM[:start]] == old_value)
    end
    return agent 
end

function move_lookahead_window!(agents::Dict,ADMM::Dict)
    ADMM[:end] = min(ADMM[:end] + 1, data["nyears"])
    for (agent,model) in agents 
        move_lookahead_window!(model,ADMM)
    end
    ADMM[:start] += 1
    mask = zeros(size(ADMM["Imbalances"]["ETS"][end]))
    mask[ADMM[:start]:ADMM[:end]] .= 1
    ADMM[:mask] = mask
    return agents
end