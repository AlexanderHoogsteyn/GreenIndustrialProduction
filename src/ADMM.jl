# Contains all functionality needed to solve the equilibrium model

include("EmissionMarket.jl")


function ADMM_iteration(mod::Model, agents::Dict, ρ::Any,scenario)
    optimize!(mod)
    for agent in values(agents)
        optimize!(agent::Model);
    end;
    prices = zeros(0,1)
    update_abatement!(mod,agents)
    for agent in values(agents)
        prices_agent = update_prices!(agent,mod)
        prices = vcat(prices,prices_agent)
    end;
    # TO DO: implement penalty term
    apply_policies!(scenario)
    return prices
end

function ADMM(mod::Model, agents::Dict, ρ::Any, δ_stop::Any,scenario::Dict)
    prices = ADMM_iteration(mod,agents,ρ)
    δ = 1
    iter = 0
    while δ > δ_stop
        prices = hcat(prices, ADMM_iteration(mod,agents,ρ))
        δ = sum((prices[:,end-1] - prices[:,end]).^2)
        iter += 1 
        if iter == 100
            print("Did not converge after " * string(iter) * " runs     δ= " * string(δ) * "\n")
            print("--------------------------------------------------\n")
            break
        end;
        print("Run: " * string(iter) * "     δ= " * string(δ) * "\n")
        print("--------------------------------------------------\n")
    end
    # final optimization to circonvent OptimizeNotCalled()
    optimize!(mod);
    for agent in values(agents)    
        optimize!(agent::Model);
    end;
    return prices
end

function apply_policies!(scenario)
    # TO DO implement policies
end;

