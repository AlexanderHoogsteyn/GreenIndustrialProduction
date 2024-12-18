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

function is_liquidity_constraint(mod::Model)
    if haskey(mod.ext[:parameters], :isLiquidityConstraint) && mod.ext[:parameters][:isLiquidityConstraint] == true
        return true
    else 
        return false
    end
end

function is_risk_averse(mod::Model)
    if haskey(mod.ext[:parameters], :isRiskAverse) && mod.ext[:parameters][:isRiskAverse] == true
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
    old_values = Dict()
    for (variable_name,  variable)  in agent.ext[:variables]
        println(variable_name)
        old_values[variable_name] = JuMP.value.(agent.ext[:variables][variable_name][ADMM[:start]])
    end
    for (variable_name, variable) in agent.ext[:variables]
        try
            delete.(agent,agent.ext[:constraints_rolling_horizon][variable_name][ADMM[:end]])
        catch error
            print("Reached end of window")
        end
        agent.ext[:constraints_rolling_horizon][variable_name][ADMM[:start]] = @constraint(agent, 
                    agent.ext[:variables][variable_name][ADMM[:start]] == old_values[variable_name]
                    )
    end
    if agent == "fringe"
        agent.ext[:parameters][:E_REF][ADMM[:end]] = data["E_ref"]
    end
    return agent 
end

function move_lookahead_window!(agents::Dict,ADMM::Dict)
    ADMM[:end] = min(ADMM[:end] + 1, data["nyears"])
    for (agent,model) in agents 
        if is_stochastic(model)
            move_lookahead_window_stochastic!(model,ADMM)
        else
            move_lookahead_window!(model,ADMM)
        end
    end
    ADMM[:start] += 1
    mask = zeros(data["nyears"])
    mask[ADMM[:start]:ADMM[:end]] .= 1
    ADMM[:mask] = mask
    agents["fringe"].ext[:parameters][:mask] = mask
    return agents
end

function move_lookahead_window_stochastic!(agent::Model,ADMM::Dict)
    # This releases the constraints on a model's decision variables on the next year, e.g. moves the lookahead window one further.
    # It fixes the decsion variable of the start of the window to its current vlue
    @assert is_stochastic(agent) " Agent is not stochastic"

    old_values = Dict()
    S = agent.ext[:sets][:S]   

    for (variable_name, variable) in agent.ext[:variables]
        if ndims(agent.ext[:variables][variable_name]) == 2
            old_values[variable_name] = JuMP.value.(agent.ext[:variables][variable_name][ADMM[:start], :]) 
        else 
            old_values[variable_name] = JuMP.value.(agent.ext[:variables][variable_name][ADMM[:start]])
        end
    end

    for (variable_name, variable) in agent.ext[:variables]
        if is_stochastic(variable)
            for s in S
                try 
                    delete.(agent,agent.ext[:constraints_rolling_horizon][variable_name][ADMM[:end],s])
                catch error
                    print("Reached end of window")
                end
                agent.ext[:constraints_rolling_horizon][variable_name][ADMM[:start],s] = @constraint(agent, agent.ext[:variables][variable_name][ADMM[:start],s] == old_values[variable_name][s])
            end
        else
            # Routine for non-stochastic variables (such as cap[y]) 
            try
                delete.(agent,agent.ext[:constraints_rolling_horizon][variable_name][ADMM[:end]])
            catch error
                print("Reached end of window")
            end
            agent.ext[:constraints_rolling_horizon][variable_name][ADMM[:start]] = @constraint(agent, agent.ext[:variables][variable_name][ADMM[:start]] == old_values[variable_name])
        end
    end
    if agent == "fringe"
        agent.ext[:parameters][:E_REF][ADMM[:end],:] .= data["E_ref"]
    end
    return agent
end