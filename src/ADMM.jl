# Contains all functionality needed to solve the equilibrium model
function ADMM!(results::Dict,ADMM::Dict,data::Dict,sector::String,agents::Dict)
    convergence = 0
    iterations = ProgressBar(1:data["max_iter"])
    nAgents = length(keys(agents))

    push!(results["s"],copy(data["S"][:]))
    push!(results["D"],copy(data["D"][:]))

    for iter in iterations
        if convergence == 0
            # Multi-threaded version
            @sync for (agent,model) in agents
                # created subroutine to allow multi-treading to solve agents' decision problems
                @spawn ADMM_subroutine!(model,data,results,ADMM,agent,nAgents)
            end

            # Imbalances 
            push!(ADMM["Imbalances"]["ETS"], results["s"][end]-sum(results["b"][m][end] for m in keys(agents)))
            push!(ADMM["Imbalances"]["product"], results["D"][end]-sum(results["g"][m][end] for m in keys(agents)))
            
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
            push!(results["λ"]["ETS"], results[ "λ"]["ETS"][end] - ADMM["ρ"]["ETS"][end]*ADMM["Imbalances"]["ETS"][end]/10)
            push!(results["λ"]["product"], results["λ"]["product"][end] + ADMM["ρ"]["product"][end]*ADMM["Imbalances"]["product"][end]/10)
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

function ADMM_subroutine!(mod::Model,data::Dict,results::Dict,ADMM::Dict,agent::String,nAgents::Int)
    # Execute common price updates
    mod.ext[:parameters][:b_bar] = results["b"][agent][end] + 1/nAgents*ADMM["Imbalances"]["ETS"][end]
    mod.ext[:parameters][:λ_ets] = results["λ"]["ETS"][end]
    mod.ext[:parameters][:ρ_ets] = ADMM["ρ"]["ETS"][end]

    mod.ext[:parameters][:g_bar] = results["g"][agent][end] + 1/nAgents*ADMM["Imbalances"]["product"][end]
    mod.ext[:parameters][:λ_product] = results["λ"]["product"][end]
    mod.ext[:parameters][:ρ_product] = ADMM["ρ"]["product"][end]

    # Myopic specific updates

    if is_myopic(mod)
        solve_myopic_agent!(mod,results,ADMM,agent,nAgents*data["nyears"])
    else
        if agent == "fringe"
            update_ind_emissions!(mod,data)
            solve_competitive_fringe!(mod,data)
        else
            solve_producer!(mod)
        end
    end
    
    # Query results
    push!(results["b"][agent], collect(value.(mod.ext[:variables][:b])))
    push!(results["e"][agent], collect(value.(mod.ext[:expressions][:netto_emiss])))
    push!(results["g"][agent], collect(value.(mod.ext[:variables][:g])))
end

function solve_myopic_agent!(mod::Model,results::Dict,ADMM::Dict,agent::String,nAgents::Int)
    @assert is_myopic(mod) "Agent is not myopic"

    mod.ext[:parameters][:g_bar_τ] = results["g_τ"][agent][end] + 1/nAgents*repeat(ADMM["Imbalances"]["product"][end],1,data["nyears"])

    if agent == "fringe"
        solve_myopic_competitive_fringe!()
    else
        solve_myopic_producer!(mod)
    end
    push!(results["g_τ"][agent], collect(value.(mod.ext[:variables_myopic][:g_τ])))

    return mod
end

function update_rho!(ADMM::Dict, iter::Int64)
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

function update_ind_emissions!(mod::Model,data::Dict)
    # Baseline emissions, corrected for share of industry in emissions 
    E_REF = data["E_ref"]*ones(data["nyears"],1)
    #E_REF[4:end] = E_REF[4:end] + data["E_ref_maritime"]*ones(data["nyears"]-3,1)

    for y = 1:data["nyears"]
        λ_nom = maximum([0,mod.ext[:parameters][:λ_ets][y]/(1+data["inflation"])^(y-1)]) # M€/MtCO2, discounted to 2021 values, limited to positive values
        mod.ext[:parameters][:e][y] = minimum([E_REF[y],maximum([0,E_REF[y] - (λ_nom/data["MAC"])^(1/data["gamma"])])]) # emissions according to MACC

        # Account for maximum price-induced change in emission allowance demand (as percentage of 2021 emissions):
        #if y > 1 && mod.ext[:parameters][:e][y-1] > data["max_em_change"]*mod.ext[:parameters][:e][1]/100 
        #    if mod.ext[:parameters][:e][y-1] - mod.ext[:parameters][:e][y] > data["max_em_change"]*mod.ext[:parameters][:e][1]/100
        #        mod.ext[:parameters][:e][y] = mod.ext[:parameters][:e][y-1]-data["max_em_change"]*mod.ext[:parameters][:e][1]/100
        #    end
        #end
        
        # Compute abatement cost
        #mod.ext[:parameters][:AC][y] = mod.ext[:parameters][:β]*(E_REF[y]-mod.ext[:parameters][:e][y])^(data["gamma"]+1)/(data["gamma"]+1)

        # Add covid or overlapping policy
       # mod.ext[:parameters][:e][y] = mod.ext[:parameters][:e][y]+data["Δe"][y]
    end

    return mod
end