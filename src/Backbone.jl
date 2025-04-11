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

function is_producer(mod::Model)
    """
    Check if the model is a producer. This is used to determine if the model should be 
    treated as a producer or not.
    """
    if haskey(mod.ext[:parameters], :isProducer) && mod.ext[:parameters][:isProducer] == true
        return true
    else 
        return false
    end
end

function set_masks!(agents::Dict, ADMM::Dict, data::Dict)
    """
    Set the mask for the rolling horizon window. This function sets the mask to true for the years
    within the rolling horizon window and false for the years outside of it.
    """
    # Set both masks up for the ADMM routine
    ADMM[:start] = 1
    ADMM[:end] = min(ADMM[:start] + data["horizon_ets"] - 1, data["nyears"])
    ADMM[:end_product] = min(ADMM[:start] + data["horizon_product"] - 1, data["nyears"])

    # Directly initialize and assign masks to the ADMM dictionary
    ADMM[:mask] = zeros(Bool, data["nyears"])
    ADMM[:mask][ADMM[:start]:ADMM[:end]] .= true

    ADMM[:mask_product] = zeros(Bool, data["nyears"])
    ADMM[:mask_product][ADMM[:start]:ADMM[:end_product]] .= true

    # Store the correct mask in each agent depending on which one applies
    for (agent, model) in agents
        model.ext[:parameters][:start] = 1
        if is_producer(model)
            model.ext[:parameters][:mask] = ADMM[:mask_product]
            model.ext[:parameters][:end] = min(ADMM[:start] + data["horizon_product"] - 1, data["nyears"])
        else
            model.ext[:parameters][:mask] = ADMM[:mask]
            model.ext[:parameters][:end] = min(ADMM[:start] + data["horizon_ets"] - 1, data["nyears"])
        end
    end
end

function set_lookahead_window!(agent::Model)
    """
    Constraints an agents decision variables outside the lookahead to zero.
    This initialises the rolling horizon window for the agent.
    """
    agent.ext[:constraints_rolling_horizon] = Dict()
    Y = agent.ext[:sets][:Y]
    Y_window = Y[agent.ext[:parameters][:end]+1:end]

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

function set_lookahead_window!(agents::Dict)
        for (agent,model) in agents
            if is_stochastic(model)
                set_lookahead_window_stochastic!(model)
            else
                set_lookahead_window!(model)
            end
        end
    return agents
end

function set_lookahead_window_stochastic!(agent::Model)
    agent.ext[:constraints_rolling_horizon] = Dict()
    Y = agent.ext[:sets][:Y]
    Y_window = Y[agent.ext[:parameters][:end]+1:end]
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

function move_lookahead_window!(agent::Model)
    # This releases the constraints on a model's decision variables on the next year, e.g. moves the lookahead window one further.
    # It fixes the decsion variable of the start of the window to its current vlue
    old_values = Dict()
    start = agent.ext[:parameters][:start]
    end_ = agent.ext[:parameters][:end]
    mask = agent.ext[:parameters][:mask]

    # Move the lookahead window one year forward
    end_ = min(end_ + 1, length(agent.ext[:sets][:Y]))

    for (variable_name,  variable)  in agent.ext[:variables]
        println(variable_name)
        old_values[variable_name] = JuMP.value.(agent.ext[:variables][variable_name][start])
    end
    for (variable_name, variable) in agent.ext[:variables]
        try
            delete.(agent,agent.ext[:constraints_rolling_horizon][variable_name][end_])
        catch error
            @assert end_ == length(agent.ext[:sets][:Y]) "Tried deleting constraint at end_, but end of window not reached"
        end
        agent.ext[:constraints_rolling_horizon][variable_name][start] = @constraint(agent, 
                    agent.ext[:variables][variable_name][start] == old_values[variable_name]
                    )
    end

    # Move the lookahead window one year forward
    start += 1
    mask[1:end] .= 0
    mask[start:end_] .= 1

    # Push updates to the agent
    agent.ext[:parameters][:start] = start
    agent.ext[:parameters][:end] = end_
    agent.ext[:parameters][:mask] = mask

    # Reference emissions are initialised to zero, so there is no imbalance outside the rolling horizon window
    # Once a year is added to the rolling horizon, the reference emissions are set to E_ref
    if agent == "fringe"
        agent.ext[:parameters][:E_REF][end_] = data["E_ref"]
    end
    return agent 
end

function move_lookahead_window!(agents::Dict,ADMM::Dict)
    # Set the new end of the window
    for (agent,model) in agents 
        if is_stochastic(model)
            move_lookahead_window_stochastic!(model)
        else
            move_lookahead_window!(model)
        end
    end

    #update the ADMM mask
    ADMM[:start] += 1
    ADMM[:end] = min(ADMM[:end] + 1, ADMM["nyears"])
    ADMM[:end_product] = min(ADMM[:end_product] + 1, ADMM["nyears"])
    # Set the new mask
    ADMM[:mask] = zeros(Bool, ADMM["nyears"])
    ADMM[:mask][ADMM[:start]:ADMM[:end]] .= true
    ADMM[:mask_product] = zeros(Bool, ADMM["nyears"])
    ADMM[:mask_product][ADMM[:start]:ADMM[:end_product]] .= true

    return agents
end

function move_lookahead_window_stochastic!(agent::Model)
    # This releases the constraints on a model's decision variables on the next year, e.g. moves the lookahead window one further.
    # It fixes the decsion variable of the start of the window to its current vlue
    @assert is_stochastic(agent) " Agent is not stochastic"

    old_values = Dict()
    S = agent.ext[:sets][:S]  
    
    start = agent.ext[:parameters][:start]
    end_ = agent.ext[:parameters][:end]
    mask = agent.ext[:parameters][:mask]

    # Move the lookahead window one year forward
    end_ = min(end_ + 1, length(agent.ext[:sets][:Y]))

    # Get the current values of the decision variables
    for (variable_name, variable) in agent.ext[:variables]
        if ndims(agent.ext[:variables][variable_name]) == 2
            old_values[variable_name] = JuMP.value.(agent.ext[:variables][variable_name][start, :]) 
        else 
            old_values[variable_name] = JuMP.value.(agent.ext[:variables][variable_name][start])
        end
    end

    # Delete the constraints on the end of the window
    for (variable_name, variable) in agent.ext[:variables]
        if is_stochastic(variable)
            for s in S
                try 
                    delete.(agent,agent.ext[:constraints_rolling_horizon][variable_name][end_,s])
                catch error
                    @assert end_ == length(agent.ext[:sets][:Y]) "Tried deleting constraint at end_, but end of window not reached"
                end
                agent.ext[:constraints_rolling_horizon][variable_name][start,s] = @constraint(agent, agent.ext[:variables][variable_name][start,s] == old_values[variable_name][s])
            end
        else
            # Routine for non-stochastic variables (such as cap[y]) 
            try
                delete.(agent,agent.ext[:constraints_rolling_horizon][variable_name][end_])
            catch error
                @assert end_ == length(agent.ext[:sets][:Y]) "Tried deleting constraint at end_, but end of window not reached"
            end
            agent.ext[:constraints_rolling_horizon][variable_name][start] = @constraint(agent, agent.ext[:variables][variable_name][start] == old_values[variable_name])
        end
    end

    # Move the lookahead window one year forward
    start += 1
    mask[1:end] .= 0
    mask[start:end_] .= 1

    # Push updates to the agent
    agent.ext[:parameters][:start] = start
    agent.ext[:parameters][:end] = end_
    agent.ext[:parameters][:mask] = mask

    if agent == "fringe"
        agent.ext[:parameters][:E_REF][end_,:] .= data["E_ref"]
    end

    return agent
end