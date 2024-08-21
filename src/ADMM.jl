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
                @spawn ADMM_subroutine!(data,results,ADMM,agent,model,nAgents)
            end
            push!(results["s"],copy(data["S"][:]))
            push!(results["D"],copy(data["D"][:]))

            # Imbalances 
                push!(ADMM["Imbalances"]["ETS"], results["s"][end]-sum(results["b"][m][end] for m in keys(agents)))
                push!(ADMM["Imbalances"]["product"], results["D"][end]-sum(results["g"][m][end] for m in keys(agents)))
                
            # Primal residuals
                push!(ADMM["Residuals"]["Primal"]["ETS"], sqrt(sum(ADMM["Imbalances"]["ETS"][end].^2)))
                push!(ADMM["Residuals"]["Primal"]["product"], sqrt(sum(ADMM["Imbalances"]["product"][end].^2)))

            # Dual residuals
            if iter > 1
                push!(ADMM["Residuals"]["Dual"]["ETS"], 
                    sqrt(sum(sum((ADMM["ρ"]["EUA"][end]*((results["b"][m][end]-sum(results["b"][mstar][end] for mstar in keys(agents))./nAgents) 
                    - (results["b"][m][end-1]-sum(results["b"][mstar][end-1] for mstar in keys(agents))./nAgents))).^2 for m in keys(agents))))
                )
                push!(ADMM["Residuals"]["Dual"]["product"], 
                    sqrt(sum(sum((ADMM["ρ"]["product"][end]*((results["g"][m][end]-sum(results["g"][mstar][end] for mstar in keys(agents))./nAgents) 
                    - (results["g"][m][end-1]-sum(results["g"][mstar][end-1] for mstar in keys(agents))./nAgents))).^2 for m in keys(agents))))
                )
            end

            # Price updates 
            push!(results["λ"]["EUA"], [data["P_2021"]; results[ "λ"]["EUA"][end][2:end] - ADMM["ρ"]["EUA"][end]*ADMM["Imbalances"]["ETS"][end][2:end]])
            push!(results["λ"]["product"], results["λ"]["product"][end] + ADMM["ρ"]["product"][end]*ADMM["Imbalances"]["product"][end])
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
            ADMM["n_iter"] = copy(iter)
        end
    end
end

function ADMM_subroutine!(data::Dict,results::Dict,ADMM::Dict,agent,mod,nAgents)

    mod.ext[:parameters][:b_bar] = results["b"][agent][end] + 1/nAgents*ADMM["Imbalances"]["ETS"][end]
    mod.ext[:parameters][:λ_EUA] = results["λ"]["EUA"][end] 
    mod.ext[:parameters][:ρ_EUA] = ADMM["ρ"]["EUA"][end]

    mod.ext[:parameters][:g_bar] = results["g"][agent][end] - 1/nAgents*ADMM["Imbalances"]["product"][end]
    mod.ext[:parameters][:λ_product] = results["λ"]["product"][end]
    mod.ext[:parameters][:ρ_product] = ADMM["ρ"]["product"][end]
    
    # Solve agents decision problems:
    if agent == "fringe"
        solve_competitive_fringe!(mod,data)
    else
        solve_producer!(mod,data,sector,agent)
    end
    
    # Query results
    push!(results["b"][agent], collect(value.(mod.ext[:variables][:b])))
    push!(results["e"][agent], collect(value.(mod.ext[:expressions][:netto_emiss])))
    push!(results["g"][agent], collect(value.(mod.ext[:variables][:g])))
end

function update_rho!(ADMM::Dict, iter::Int64)
    if mod(iter,1) == 0
        # ρ-updates following Boyd et al.  
        if ADMM["Residuals"]["Primal"]["ETS"][end]> 2*ADMM["Residuals"]["Dual"]["ETS"][end]
            push!(ADMM["ρ"]["EUA"], minimum([1000,1.1*ADMM["ρ"]["EUA"][end]]))
        elseif ADMM["Residuals"]["Dual"]["ETS"][end] > 2*ADMM["Residuals"]["Primal"]["ETS"][end]
            push!(ADMM["ρ"]["EUA"], 1/1.1*ADMM["ρ"]["EUA"][end])
        end
        if ADMM["Residuals"]["Primal"]["product"][end]> 2*ADMM["Residuals"]["Dual"]["product"][end]
            push!(ADMM["ρ"]["product"], minimum([1000,1.1*ADMM["ρ"]["product"][end]]))
        elseif ADMM["Residuals"]["Dual"]["product"][end] > 2*ADMM["Residuals"]["Primal"]["product"][end]
            push!(ADMM["ρ"]["product"], 1/1.1*ADMM["ρ"]["product"][end])
        end
    end
end