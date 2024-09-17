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