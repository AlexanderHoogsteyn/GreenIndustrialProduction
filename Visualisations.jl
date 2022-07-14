# Essentials to run, needed for every little demo below
assume = YAML.load_file(joinpath(@__DIR__, "data_assumptions.yaml"));
properties = YAML.load_file(joinpath(@__DIR__, "data_properties.yaml"));

SteelRoutes = Dict{String,Dict}(assume["SteelRoutes"])
CommodityPrices = Dict{String,Float64}(assume["CommodityPrices"])

# 1. Influence of hydrogen price on cost of steel
hydrogen_price_range_kg = 1:6 # €/kg
hydrogen_price_range = hydrogen_price_range_kg/properties["Hydrogen"]["EnergyDensity"]*1000 # convert to €/MWh

prices = vary_need("Hydrogen", hydrogen_price_range, CommodityPrices);
no_CCfD = Dict{String, Dict{String, Float64}}("CCfD"=>Dict("Strike"=>0.0, "Benchmark"=>1.850));

x = hydrogen_price_range_kg
y = [route_costs(prices,SteelRoutes["BF-BOF"],no_CCfD),
     route_costs(prices,SteelRoutes["BF-Blend20-BOF"],no_CCfD),
     route_costs(prices,SteelRoutes["HDRI-EAF"],no_CCfD)
    ]

plot(x, y, title="Cost of steel-making routes",
    label = ["BF-BOF" "BF-Blend20-BOF" "HDRI-EAF"],
    xlabel = "Hydrogen price (€/kg)",
    ylabel = "Cost of ton of steel(€)",
    lw = 2)

# 2. Influenece of ETS price on cost of steel
eua_price_range = 10:1:200 # €/kg


prices = vary_need("EUA", eua_price_range, CommodityPrices);

CCfD = Dict{String, Dict{String, Float64}}("CCfD"=>Dict("Strike"=>101.03,"Benchmark"=>1.850))

x = eua_price_range
y = [route_costs(prices,SteelRoutes["BF-BOF"],no_CCfD),
     route_costs(prices,SteelRoutes["BF-Blend20-BOF"],CCfD),
     route_costs(prices,SteelRoutes["HDRI-EAF"],CCfD)
    ]

plot(x, y, title="Cost of steel-making routes",
    label = ["BF-BOF" "BF-Blend20-BOF" "HDRI-EAF"],
    xlabel = "EU-ETS price (€/ton)",
    ylabel = "Cost of ton of steel(€)",
    legend=:topleft,
    lw = 2)



# 3. Optimal strike price

no_policy = Dict{String, Dict{String, Float64}}()
optimal_strike_price(CommodityPrices,SteelRoutes["HDRI-EAF"],SteelRoutes["BF-BOF"])

# 4. Calculationg Windfall profits for ETS price range
route = SteelRoutes["HDRI-EAF"]
benchmark_route = SteelRoutes["BF-BOF"]
CCfD = optimal_CCfD(CommodityPrices,route,benchmark_route)
windfall_profit(CCfD, CommodityPrices, route, benchmark_route)

eua_price_range = 10:1:200 # €/kg
prices = vary_need("EUA", eua_price_range, CommodityPrices);
ev, vector = windfall_profit(optimal_CCfD(CommodityPrices,route,benchmark_route), prices, route, benchmark_route)

x = eua_price_range
y = [route_costs(prices,route,no_policy),
     route_costs(prices,benchmark_route,CCfD),
     vector
    ]

plot(x, y, title="Cost of steel-making routes",
    label = ["BF-BOF" "HDRI-EAF" "windfall profit"],
    xlabel = "EU-ETS price (€/ton)",
    ylabel = "Cost of ton of steel(€)",
    legend=:topleft,
    lw = 2)

# 5. Calculating windfall profits for stocastic ETS price & electricity price
route = SteelRoutes["HDRI-EAF"]
benchmark_route = SteelRoutes["BF-BOF"]
CCfD_BF = Dict{String, Dict{String, Float64}}("CCfD"=>Dict("Strike"=>70,"Benchmark"=>1.850),
                                                "Grandfathering"=>Dict("Benchmark"=>1.2880))
CCfD_HDRI = optimal_CCfD(CommodityPrices,route,benchmark_route)
CCfD_HDRI["Grandfathering"] = Dict("Benchmark"=>1.288)

no_policy = Dict{String, Dict{String, Float64}}("Grandfathering"=>Dict("Benchmark"=>1.2880))

prices = random_vary_need("EUA", 70,10,1000, CommodityPrices);

y = DataFrame()
y."BF-BOF" = route_costs(prices,benchmark_route,no_policy)
y."HDRI-EAF" = route_costs(prices,route,no_policy)
#y."ETS-benchmark" = route_costs(prices,SteelRoutes["ETS-Benchmark"],no_policy)
y."HDRI-EAF CCfD" = route_costs(prices,route,CCfD_HDRI)
#CSV.write("costs_w&wo_CCfD.csv",y)
unpivot = stack(y);
y."sensitivity" .= "EUA";
#@df unpivot violin(string.(:variable), :value, fillalpha=0.75, linewidth=2, ylims=[0,700])
#@df unpivot boxplot(string.(:variable), :value, fillalpha=0.8, linewidth=2,
#label="EUA price sensitivity",ylabel="Cost (€/t)")

# Show that even with a perfectly priced CCfD if market conditions change they do not hedge all risk
prices = random_vary_need("Electricity", 40,10,1000, CommodityPrices);
z = DataFrame() 
z."BF-BOF" = route_costs(prices,benchmark_route,no_policy)
z."HDRI-EAF" = route_costs(prices,route,no_policy)
#z."Blend 20% H2 BF-BOF" = route_costs(prices,SteelRoutes["BF-Blend20-BOF"])
#z."ETS-benchmark" = route_costs(prices,SteelRoutes["ETS-Benchmark"],no_policy)
z."HDRI-EAF CCfD" = route_costs(prices,route,CCfD_HDRI)
#CSV.write("costs_w&wo_CCfD.csv",y)
unpivot = stack(z);
z."sensitivity" .= "Electricity";
append!(y,z)
unpivot = stack(y)
#@df unpivot violin(string.(:variable), :value, fillalpha=0.75, linewidth=2, ylims=[0,700])
@df unpivot groupedboxplot(string.(:variable), :value,group=:sensitivity, fillalpha=0.8, linewidth=2, label=["ETS-price sensitivity" "Electricity price sensitivity"])
#@df unpivot groupedviolin(string.(:variable), :value,group=:sensitivity, fillalpha=0.8, linewidth=2, label=["ETS-price sensitivity" "Electricity price sensitivity"])



# 6. Show cost of different routes component wise


    # Costs w/o CCfD
x = DataFrame()
for (component,amount) in CommodityPrices
    x_append = DataFrame()
    x_append."BF-BOF" = [cost_component(CommodityPrices,benchmark_route,component)]
    x_append."HDRI-EAF" = [cost_component(CommodityPrices,route,component)]
    x_append."Blend 20% H2 BF-BOF" = [cost_component(CommodityPrices,SteelRoutes["BF-Blend20-BOF"],component)]
    x_append."Component" = [component]
    append!(x,x_append)
end



x_trans = permutedims(x, 1)

    # Costs w CCfD
#y."HDRI-EAF CCfD" = route_costs(prices,route,CCfD_HDRI)
#y."BF-BOF CCfD" = route_costs(prices,benchmark_route,CCfD_BF)
CSV.write("C:\\\\Users\\u0137718\\OneDrive\\Documenten\\Julia\\GreenSteelAnalysis\\DesignCriteriaForCCfDs\\Data\\costs_components.csv",x)

unpivot = stack(x)

    # Plots
@df unpivot groupedbar(string.(:variable), :Component, bar_position=:value, fillalpha=0.8, linewidth=2) # label=collect(keys(CommodityPrices)))