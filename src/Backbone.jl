#-------------------------------------------------------------------------
# First estimation of saving potential of innovative green steel producer
#   Author: Alexander Hoogsteyn
#   Date: March 2022
#
#
#-------------------------------------------------------------------------
#import Pkg;
#Pkg.activate(@__DIR__) # @__DIR__ = directory this script is in
#Pkg.update()
#Pkg.instantiate()
#Pkg.add("GR")
#Pkg.add("CSV")
#Pkg.add("DataFrames")
#Pkg.add("YAML")
#Pkg.add("Plots")
#Pkg.add("Distributions")
#Pkg.add("Random")
#Pkg.add("StatsPlots")

using CSV
using DataFrames
using YAML
using Plots
using Statistics
using Distributions, Random
using StatsPlots

# Expects hydrogen price (or range) in €/kg
function needed_eua_price(hydrogen_price,benchmark_route::Dict,route::Dict)
    return 1000*hydrogen_price*route["Needs"]["Hydrogen"] / (benchmark_route["Needs"]["ETS"]-route["Needs"]["ETS"] )
end

function route_costs(commodityPrices::Dict{Any,Any}, route::Dict{Any,Any},policies::Dict{Any,Any})
    cost = 0;
    needs = copy(route)
    needs = delete!(needs,"CAPEX")
    needs = delete!(needs,"ETS")


    # Account for policies that reduce needs (e.g. grandfathering of ETS)
    #if "Grandfathering" in keys(policies)
    #    needs["ETS"] = needs["ETS"] - policies["Grandfathering"]["Benchmark"]
    #end

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

function vary_need(need::String,price_range,commodity_prices::Dict{String,Float64})
    commodity_prices_range = Dict{Int64,Dict{String,Float64}}();
    for (i, price) in enumerate(price_range)
        commodity_prices_range[i] = copy(commodity_prices);
        commodity_prices_range[i][need] = price;
    end
    return commodity_prices_range
end

function random_vary_need(need::String,avg,std,draws,commodity_prices::Dict{String,Float64})
    d = Normal(avg,std)
    price_range = rand(d,draws)
    return vary_need(need,price_range,commodity_prices)
end

function floor_revenue(policy::Dict{String, Float64}, commodity_prices::Dict{String,Float64},route::Dict{Any,Any})
    if policy["Strike"] > commodity_prices["ETS"]
        return (policy["Strike"] - commodity_prices["ETS"])*(policy["Benchmark"] - route["Needs"]["ETS"])
    else
        return 0
    end
end

# CCfD functionalities

function CCfD_revenue(policy::Dict{String, Float64}, commodity_prices::Dict{String,Float64},route::Dict{Any,Any})
        return (policy["Strike"] - commodity_prices["ETS"])*(policy["Benchmark"] - route["Needs"]["ETS"])
end


function optimal_strike_price(commodity_prices::Dict{String,Float64},route::Dict{Any,Any},benchmark_route::Dict{Any,Any})
    no_policy = Dict{String, Dict{String, Float64}}()
    subsidy = (route_costs(commodity_prices,route,no_policy) - route_costs(commodity_prices,benchmark_route,no_policy))/
    1.288
    return commodity_prices["ETS"] + subsidy    
end

function optimal_strike_price(commodidy_prices::Dict{Int64,Dict{String,Float64}}, route::Dict{Any,Any},benchmark_route::Dict{Any,Any})
    prices = Vector{Float64}(undef,length(commodidy_prices));
    for (key,dict) in commodidy_prices
        prices[key] = optimal_strike_price(dict, route, benchmark_route);
    end
    return prices
end

function optimal_CCfD(commodity_prices::Dict{String,Float64},route::Dict{Any,Any},benchmark_route::Dict{Any,Any})
    optimal_strike = optimal_strike_price(commodity_prices::Dict{String,Float64},route::Dict{Any,Any},benchmark_route::Dict{Any,Any})
    benchmark = benchmark_route["Needs"]["ETS"]
    return Dict{String, Dict{String, Float64}}("CCfD"=>Dict("Strike"=>optimal_strike,"Benchmark"=>1.288))
end

# Functiuanality to calculate social welfare

function windfall_profit(policies::Dict{String, Dict{String, Float64}}, commodity_prices::Dict{String,Float64},route::Dict{Any,Any}, benchmark_route::Dict{Any,Any})
    no_policies = Dict{String, Dict{String, Float64}}();
    return route_costs(commodity_prices,route,policies) - route_costs(commodity_prices,benchmark_route,no_policies)
end

function windfall_profit(policies::Dict{String, Dict{String, Float64}}, commodity_prices::Dict{Int64,Dict{String,Float64}},route::Dict{Any,Any}, benchmark_route::Dict{Any,Any})
    no_policies = Dict{String, Dict{String, Float64}}();
    profits = Vector{Float64}(undef,length(commodity_prices));
    for (key,dict) in commodity_prices
        profits[key] = route_costs(dict,route,policies) - route_costs(dict,benchmark_route,no_policies);
    end
    return mean(profits), profits
end