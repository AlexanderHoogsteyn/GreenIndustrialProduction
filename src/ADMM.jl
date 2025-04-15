function ADMM!(results::Dict, ADMM::Dict, data::Dict, agents::Dict)
    """
    ADMM!(results::Dict, ADMM::Dict, data::Dict, agents::Dict)

    Solves the equilibrium model using the ADMM algorithm. This function iterates until convergence is achieved or the maximum number of iterations is reached.

    # Arguments
    - `results::Dict`: Dictionary to store the results.
    - `ADMM::Dict`: Dictionary containing ADMM parameters and residuals.
    - `data::Dict`: Dictionary containing input data.
    - `agents::Dict`: Dictionary containing agent models.
    """
    convergence = 0
    if data["printoutlevel"] > 0 
        iterations = ProgressBar(1:data["max_iter"])
    else 
        iterations = 1:data["max_iter"]
    end
    nAgents = length(keys(agents))

    for iter in iterations
        if convergence == 0
            # Multi-threaded version
            @sync for (agent, model) in agents
                # created subroutine to allow multi-threading to solve agents' decision problems
                @spawn ADMM_subroutine!(model, data, results, ADMM, agent, nAgents)
            end

            # Update imbalances
            update_imbalances!(results, ADMM, agents)

            # Extract common parameters
            ADMM_ρ = ADMM["ρ"]
            primal_ets_residual = ADMM["Residuals"]["Primal"]["ETS"]
            primal_product_residual = ADMM["Residuals"]["Primal"]["product"]

            # Primal residuals
            push!(primal_ets_residual, sqrt(sum(ADMM["Imbalances"]["ETS"][end].^2)))
            push!(primal_product_residual, sqrt(sum(ADMM["Imbalances"]["product"][end].^2)))

            # Dual residuals
            if iter > 1
                push!(ADMM["Residuals"]["Dual"]["ETS"], 
                    sqrt(sum(sum((ADMM_ρ["ETS"][end]*((results["b"][m][end]-sum(results["b"][mstar][end] for mstar in keys(agents))./nAgents) 
                    - (results["b"][m][end-1]-sum(results["b"][mstar][end-1] for mstar in keys(agents))./nAgents))).^2 for m in keys(agents))))
                )
                push!(ADMM["Residuals"]["Dual"]["product"], 
                    sqrt(sum(sum((ADMM_ρ["product"][end]*((results["g"][m][end]-sum(results["g"][mstar][end] for mstar in keys(agents))./nAgents) 
                    - (results["g"][m][end-1]-sum(results["g"][mstar][end-1] for mstar in keys(agents))./nAgents))).^2 for m in keys(agents))))
                )
            end

            # Price updates 
            update_prices!(results, ADMM)
 
            # Update ρ-values
            update_rho!(ADMM, iter)

            # Update progress bar description (only every few iterations for efficiency)
            if iter % 1 == 0 && data["printoutlevel"] > 0
                set_description(
                    iterations,
                    string(@sprintf("ΔETS: %.3f -- Δproduct: %.3f ",
                        primal_ets_residual[end] / ADMM["Tolerance"]["ETS"],
                        primal_product_residual[end] / ADMM["Tolerance"]["product"]))
                )
            end            
            # Check convergence: primal and dual satisfy tolerance 
            if primal_ets_residual[end] <= ADMM["Tolerance"]["ETS"] && ADMM["Residuals"]["Dual"]["ETS"][end] <= ADMM["Tolerance"]["ETS"] && 
                primal_product_residual[end] <= ADMM["Tolerance"]["product"] && ADMM["Residuals"]["Dual"]["product"][end] <= ADMM["Tolerance"]["product"]
                    convergence = 1
            end
            # store number of iterations
        else
            ADMM["n_iter"] = copy(iter)
        end
    end
end

function ADMM_rolling_horizon!(results::Dict, ADMM::Dict, data::Dict, agents::Dict)
    """
    ADMM_rolling_horizon!(results::Dict, ADMM::Dict, data::Dict, agents::Dict)

    Solves the equilibrium model using the ADMM algorithm with a rolling horizon approach.

    # Arguments
    - `results::Dict`: Dictionary to store the results.
    - `ADMM::Dict`: Dictionary containing ADMM parameters and residuals.
    - `data::Dict`: Dictionary containing input data.
    - `agents::Dict`: Dictionary containing agent models.
    """
    ADMM[:isRollingHorizon] = true

    # Define start and end indices for ETS and product horizons
    set_masks!(agents, ADMM, data)
 
    set_lookahead_window!(agents)

    while ADMM[:end] < data["nyears"] || ADMM[:end_product] < data["nyears"]
        ADMM!(results, ADMM, data, agents)
        move_lookahead_window!(agents,ADMM)
        #println("Move window to " * string(ADMM[:start]) * ":" * string(ADMM[:end]))
    end
    ADMM!(results, ADMM, data, agents)
end

function ADMM_dual_rolling_horizon!(results::Dict, ADMM::Dict, data::Dict)
    """
    dual_rolling_horizon!(results::Dict, ADMM::Dict, data::Dict, agents::Dict)

    Solves the equilibrium model using the dual rolling horizon approach.

    # Arguments
    - `results::Dict`: Dictionary to store the results.
    - `ADMM::Dict`: Dictionary containing ADMM parameters and residuals.
    - `data::Dict`: Dictionary containing input data.
    - `agents::Dict`: Dictionary containing agent models.
    """
    ADMM[:isRollingHorizon] = true
    ADMM[:isDualRollingHorizon] = true

    # Initialize λ_ETS as a circular buffer
    global λ_ETS = CircularBuffer{Array{Float64,2}}(data["CircularBufferSize"]) 
    push!(λ_ETS, zeros(data["nyears"],data["nsamples"]))  # Initialize with zeros (or another default value)

    iter = 1

    ϵ = 0.1

    agents = Dict()

    
    while iter == 1 || sqrt(mean((λ_ETS[end] - λ_ETS[end-1]).^2)) > ϵ
        # Rebuild model
        agents["trader"] = build_stochastic_liquidity_constraint_trader!( Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))), data)
        agents["fringe"] = build_stochastic_competitive_fringe!( Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))), data)
        for (route, dict) in data["technologies"]
            agents[route] = build_stochastic_producer!( Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))), data, route)
        end

        ADMM_rolling_horizon!(results, ADMM, data, agents)


        push!(λ_ETS, results["λ"]["ETS"][end])

        println("Residual = ", sqrt(mean((λ_ETS[end] - λ_ETS[end-1]).^2)))

        iter += 1
        # Extract ETS Price
        #define_results_hot_start!(data,results,ADMM)
    end

    return agents, results
end

function ADMM_subroutine!(mod::Model, data::Dict, results::Dict, ADMM::Dict, agent::String, nAgents::Int)
    """
    ADMM_subroutine!(mod::Model, data::Dict, results::Dict, ADMM::Dict, agent::String, nAgents::Int)

    Subroutine to solve the decision problem for each agent in parallel.

    # Arguments
    - `mod::Model`: The optimization model for the agent.
    - `data::Dict`: Dictionary containing input data.
    - `results::Dict`: Dictionary to store the results.
    - `ADMM::Dict`: Dictionary containing ADMM parameters and residuals.
    - `agent::String`: The agent identifier.
    - `nAgents::Int`: The number of agents.
    """
    # Extract common parameters

    results_b = results["b"][agent]
    results_g = results["g"][agent]
    results_e = results["e"][agent]
    results_λ = results["λ"]
    ADMM_ρ = ADMM["ρ"]

    # Execute common price updates
    mod.ext[:parameters][:b_bar] = results_b[end] + 1/nAgents*ADMM["Imbalances"]["ETS"][end]
    mod.ext[:parameters][:λ_ets] = results_λ["ETS"][end]
    mod.ext[:parameters][:ρ_ets] = ADMM_ρ["ETS"][end]

    mod.ext[:parameters][:g_bar] = results_g[end] + 1/nAgents*ADMM["Imbalances"]["product"][end]
    mod.ext[:parameters][:λ_product] = results_λ["product"][end]
    mod.ext[:parameters][:ρ_product] = ADMM_ρ["product"][end]

    if agent == "fringe"
        solve_stochastic_competitive_fringe!(mod, data)
    elseif agent == "trader"
        solve_stochastic_trader!(mod, data)
    else
        solve_stochastic_producer!(mod)
    end
    
    # Query results
    if agent == "fringe"
        push!(results_e, value.(mod.ext[:variables][:e]))
        push!(results["π_MAC"]["fringe"], value.(mod.ext[:expressions][:π_MAC]))
    else 
        push!(results_e, value.(mod.ext[:expressions][:netto_emiss]))
    end
    push!(results_b, value.(mod.ext[:variables][:b]))
    push!(results_g, value.(mod.ext[:variables][:g]))
end


function update_imbalances!(results::Dict, ADMM::Dict, agents::Dict)
    """
    update_imbalances!(results::Dict, ADMM::Dict, agents::Dict)

    Updates the imbalances for ETS and product based on the current results.

    # Arguments
    - `results::Dict`: Dictionary to store the results.
    - `ADMM::Dict`: Dictionary containing ADMM parameters and residuals.
    - `agents::Dict`: Dictionary containing agent models.
    """
    if is_rolling_horizon(ADMM)
        push!(ADMM["Imbalances"]["ETS"], (results["s"].-sum(results["b"][m][end] for m in keys(agents))).* ADMM[:mask])
        push!(ADMM["Imbalances"]["product"], (results["D"].-sum(results["g"][m][end] for m in keys(agents))).* ADMM[:mask_product])
    else
        push!(ADMM["Imbalances"]["ETS"], results["s"].-sum(results["b"][m][end] for m in keys(agents)))
        push!(ADMM["Imbalances"]["product"], results["D"].-sum(results["g"][m][end] for m in keys(agents)))
    end
    return results
end


function update_rho!(ADMM::Dict, iter::Int64)
    """
    update_rho!(ADMM::Dict, iter::Int64)

    Dynamically adjusts the 'rho' tuning parameter in ADMM to potentially lead to faster convergence.

    # Arguments
    - `ADMM::Dict`: Dictionary containing ADMM parameters and residuals.
    - `iter::Int64`: The current iteration number.
    """
    if mod(iter,1) == 0
        # ρ-updates following Boyd et al.  
        if ADMM["Residuals"]["Primal"]["ETS"][end]> 2*ADMM["Residuals"]["Dual"]["ETS"][end]
            push!(ADMM["ρ"]["ETS"], minimum([1000,1.1*ADMM["ρ"]["ETS"][end]]))
        elseif ADMM["Residuals"]["Dual"]["ETS"][end] > 2*ADMM["Residuals"]["Primal"]["ETS"][end]
            push!(ADMM["ρ"]["ETS"], 1/1.1*ADMM["ρ"]["ETS"][end])
        end
        if ADMM["Residuals"]["Primal"]["product"][end]> 2*ADMM["Residuals"]["Dual"]["product"][end]
            push!(ADMM["ρ"]["product"], minimum([1000,1.1*ADMM["ρ"]["product"][end]]))
        elseif ADMM["Residuals"]["Dual"]["product"][end] > 2*ADMM["Residuals"]["Primal"]["product"][end]
            push!(ADMM["ρ"]["product"], 1/1.1*ADMM["ρ"]["product"][end])
        end
    end
end

function update_prices!(results::Dict, ADMM::Dict)
    """
    update_prices!(results::Dict, ADMM::Dict)

    Updates the prices for ETS and product based on the current imbalances and rho values.

    # Arguments
    - `results::Dict`: Dictionary to store the results.
    - `ADMM::Dict`: Dictionary containing ADMM parameters and residuals.
    """
    if is_rolling_horizon(ADMM)
        push!(results["λ"]["ETS"], max.(0,results[ "λ"]["ETS"][end] - (ADMM["ρ"]["ETS"][end]*ADMM["Imbalances"]["ETS"][end]/1000).* ADMM[:mask] ))
        push!(results["λ"]["product"], results["λ"]["product"][end] + (ADMM["ρ"]["product"][end]*ADMM["Imbalances"]["product"][end]/1000).* ADMM[:mask_product])
    else
        push!(results["λ"]["ETS"], max.(0,results[ "λ"]["ETS"][end] - ADMM["ρ"]["ETS"][end]*ADMM["Imbalances"]["ETS"][end]/1000))
        push!(results["λ"]["product"], results["λ"]["product"][end] + ADMM["ρ"]["product"][end]*ADMM["Imbalances"]["product"][end]/1000) 
    end
end