# Contains all functionality needed to solve the equilibrium model
function ADMM!(results::Dict,ADMM::Dict,data::Dict,sector::String,agents::Dict)
    convergence = 0
    iterations = ProgressBar(1:data["max_iter"])
    nAgents = length(keys(agents))

    for iter in iterations
        if convergence == 0
            # Multi-threaded version
            @sync for (agent,model) in agents
                # created subroutine to allow multi-treading to solve agents' decision problems
                @spawn ADMM_subroutine!(model,data,results,ADMM,agent,nAgents)
            end

            update_imbalances!(results,ADMM,agents)

            # Primal residuals
            push!(ADMM["Residuals"]["Primal"]["ETS"], sqrt(sum(ADMM["Imbalances"]["ETS"][end].^2)))
            push!(ADMM["Residuals"]["Primal"]["product"], sqrt(sum(ADMM["Imbalances"]["product"][end].^2)))

            # Dual residuals
            if iter > 1
                push!(ADMM["Residuals"]["Dual"]["ETS"], 
                    sqrt(sum(sum((ADMM["ρ"]["ETS"][end]*((results["b"][m][end]-sum(results["b"][mstar][end] for mstar in keys(agents))./nAgents) 
                    - (results["b"][m][end-1]-sum(results["b"][mstar][end-1] for mstar in keys(agents))./nAgents))).^2 for m in keys(agents))))
                )
                push!(ADMM["Residuals"]["Dual"]["product"], 
                    sqrt(sum(sum((ADMM["ρ"]["product"][end]*((results["g"][m][end]-sum(results["g"][mstar][end] for mstar in keys(agents))./nAgents) 
                    - (results["g"][m][end-1]-sum(results["g"][mstar][end-1] for mstar in keys(agents))./nAgents))).^2 for m in keys(agents))))
                )
            end

            # Price updates 
            update_prices!(results,ADMM)
 
            # Update ρ-values
            update_rho!(ADMM,iter)

            # Progress bar
            set_description(iterations, string(@sprintf("ΔETS: %.3f -- Δproduct: %.3f ",  ADMM["Residuals"]["Primal"]["ETS"][end]/ADMM["Tolerance"]["ETS"], ADMM["Residuals"]["Primal"]["product"][end]/ADMM["Tolerance"]["product"])))
            # Check convergence: primal and dual satisfy tolerance 
            if ADMM["Residuals"]["Primal"]["ETS"][end] <= ADMM["Tolerance"]["ETS"] && ADMM["Residuals"]["Dual"]["ETS"][end] <= ADMM["Tolerance"]["ETS"] && 
               ADMM["Residuals"]["Primal"]["product"][end] <= ADMM["Tolerance"]["product"] && ADMM["Residuals"]["Dual"]["product"][end] <= ADMM["Tolerance"]["product"]
                    convergence = 1
            end
            # store number of iterations
        else
            ADMM["n_iter"] = copy(iter)
        end
    end
end

function ADMM_rolling_horizon!(results::Dict,ADMM::Dict,data::Dict,sector::String,agents::Dict)
    ADMM[:isRollingHorizon] = true

    ADMM[:start] = 1
    ADMM[:end]  = min(ADMM[:start] + data["horizon_ets"]-1, data["nyears"])
    mask = zeros(data["nyears"])
    mask[ADMM[:start]:ADMM[:end]] .= 1
    ADMM[:mask] = mask
 
    set_lookahead_window!(agents,ADMM)

    while ADMM[:end] < data["nyears"]
        ADMM!(results,ADMM,data,sector,agents)
        move_lookahead_window!(agents,ADMM)
        #println("Move window to " * string(ADMM[:start]) * ":" * string(ADMM[:end]))
    end
    ADMM!(results,ADMM,data,sector,agents)
end

function ADMM_subroutine!(mod::Model,data::Dict,results::Dict,ADMM::Dict,agent::String,nAgents::Int)
    # Execute common price updates
    mod.ext[:parameters][:b_bar] = results["b"][agent][end] + 1/nAgents*ADMM["Imbalances"]["ETS"][end]
    mod.ext[:parameters][:λ_ets] = results["λ"]["ETS"][end]
    mod.ext[:parameters][:ρ_ets] = ADMM["ρ"]["ETS"][end]

    mod.ext[:parameters][:g_bar] = results["g"][agent][end] + 1/nAgents*ADMM["Imbalances"]["product"][end]
    mod.ext[:parameters][:λ_product] = results["λ"]["product"][end]
    mod.ext[:parameters][:ρ_product] = ADMM["ρ"]["product"][end]

    if agent == "fringe"
        if is_stochastic(mod)
            solve_stochastic_competitive_fringe!(mod)
        else
            solve_competitive_fringe!(mod)
        end
    elseif agent == "trader"
        if is_stochastic(mod)
            solve_stochastic_trader!(mod,data)
        else
            solve_trader!(mod)
        end
    else
        if is_stochastic(mod)
            solve_stochastic_producer!(mod)
        else
            solve_producer!(mod)
        end
    end
    
    # Query results
    push!(results["b"][agent], collect(value.(mod.ext[:variables][:b])))
    push!(results["e"][agent], collect(value.(mod.ext[:expressions][:netto_emiss])))
    push!(results["g"][agent], collect(value.(mod.ext[:variables][:g])))
end

function update_imbalances!(results::Dict,ADMM::Dict,agents::Dict)
    if is_rolling_horizon(ADMM)
        push!(ADMM["Imbalances"]["ETS"], (results["s"].-sum(results["b"][m][end] for m in keys(agents))).* ADMM[:mask])
        push!(ADMM["Imbalances"]["product"], (results["D"].-sum(results["g"][m][end] for m in keys(agents))).* ADMM[:mask])
    else
        push!(ADMM["Imbalances"]["ETS"], results["s"].-sum(results["b"][m][end] for m in keys(agents)))
        push!(ADMM["Imbalances"]["product"], results["D"].-sum(results["g"][m][end] for m in keys(agents)))
    end
    return results
end

function update_rho!(ADMM::Dict, iter::Int64)
    #=
    Implements dynamic adjustments of 'rho' tuning parameter in ADMM. 
    This should in theory lead to faster convergence.
    =#
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

function update_prices!(results::Dict,ADMM::Dict)
    if is_rolling_horizon(ADMM)
        push!(results["λ"]["ETS"], max.(0,results[ "λ"]["ETS"][end] - (ADMM["ρ"]["ETS"][end]*ADMM["Imbalances"]["ETS"][end]/10).* ADMM[:mask] ))
        push!(results["λ"]["product"], results["λ"]["product"][end] + (ADMM["ρ"]["product"][end]*ADMM["Imbalances"]["product"][end]/10).* ADMM[:mask])
    else
        push!(results["λ"]["ETS"], max.(0,results[ "λ"]["ETS"][end] - ADMM["ρ"]["ETS"][end]*ADMM["Imbalances"]["ETS"][end]/10))
        push!(results["λ"]["product"], results["λ"]["product"][end] + ADMM["ρ"]["product"][end]*ADMM["Imbalances"]["product"][end]/10) 
    end
end