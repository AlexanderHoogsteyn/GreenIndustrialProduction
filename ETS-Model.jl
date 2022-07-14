import Pkg;
#Pkg.add("JuMP")
#Pkg.add("Gurobi")
#Pkg.add("YAML")
#Pkg.add("Plots")
#modPkg.add("DataFrames")


using JuMP
using Gurobi
using YAML
using Plots
using DataFrames
using Plots
# If there exist different carbon reducing technologies with different CAPEX OPEX (efficiencies) what is an optimal policy

mod = Model(optimizer_with_attributes(Gurobi.Optimizer))

data = YAML.load_file(joinpath(@__DIR__, "data_assumptions_ets.yaml"));

# Parameters/variables ETS 
ETS = Dict()
function define_ETS_parameters!(ETS::Dict,data::Dict)
    # LRF 2017 - 2021, no Brexit, expansion of scope (aviation, maritime) or Green Deal
    ETS["LRF"] = zeros(data["nyears"],1);
    ETS["LRF"][1:4] = data["LRF_2017"]*ones(4,1); 
    ETS["LRF"][5:14] = ETS["LRF"][1]*data["LRF_2021"]/0.0174*ones(10,1);                            # 2021-2030
    ETS["LRF"][15:end] = ETS["LRF"][1]*data["LRF_2031"]/0.0174*ones(data["nyears"]-14,1);           # 2030 - end ETS
       
    # CAP
    ETS["CAP"] = zeros(data["nyears"],1);
    for y =1:data["nyears"]
        ETS["CAP"][y]= maximum([data["CAP_2016"]-sum(ETS["LRF"][1:y]) 0])
    end
    ETS["CAP"][1] = data["S_2017"]+data["EX_2016"] # real supply in 2017 and surplus in market at end 2016

    # Supply 
    ETS["S"] = copy(ETS["CAP"])

    return ETS
end
define_ETS_parameters!(ETS,data)

function define_sets_parameters!(mod::Model, data::Dict)
    
    # Define sets
    mod.ext[:sets] = Dict()
    mod.ext[:sets][:Y] = 1:data["nyears"]
    Y = mod.ext[:sets][:Y]

    # define parameters 
    mod.ext[:parameters] = Dict()

    # Emissions representative agents, bound to historical values in 2017-2019
    mod.ext[:parameters][:e] = ones(data["nyears"],1)*data["E_2019"]
    #mod.ext[:parameters][:e][1] = data["E_2017"]
    #mod.ext[:parameters][:e][2] = data["E_2018"]
    #mod.ext[:parameters][:e][3] = data["E_2019"]
    E = mod.ext[:parameters][:e]

    mod.ext[:parameters][:A] = ones(data["nyears"],1)
    for y in 4:data["nyears"]
        mod.ext[:parameters][:A][y] = 1/(1+data["discount_rate"])^(y-3);
    end
    A = mod.ext[:parameters][:A]
    
    # Price structure
    mod.ext[:parameters][:λ] = zeros(data["nyears"],1)

    # Define variables
    mod.ext[:variables] = Dict() 
    cap = mod.ext[:variables][:cap] = @variable(mod, [y=Y], lower_bound=0, base_name="capacity")
    #g = mod.ext[:variables][:g] = @variable(mod, [y=Y], lower_bound=0, base_name="generation")
    #buy = mod.ext[:variables][:buy] =  @variable(mod, [y=Y], lower_bound=0, base_name="allowances bougth")
    a = mod.ext[:variables][:a] =  @variable(mod, [y=Y], lower_bound=0, base_name="abatement")

    # Define exprsessions
    mod.ext[:expressions] = Dict()
    mod.ext[:expressions][:bank] = @expression(mod, [y=Y], sum(ETS["S"][1:y])+sum(a[1:y])-sum(E[1:y]))
    bank = mod.ext[:expressions][:bank]
    mod.ext[:expressions][:netto_emiss] = @expression(mod, [y=Y], E[y]-a[y])
    mod.ext[:expressions][:elek] = @expression(mod, [y=Y], a[y]*data["efficiency_PEM"])

    CAPEX = 1000 #
    OPEX = 1
    MCD = 100 # Marginal carbon damage
    MAC = 20 # Marginal abatement cost
    # Objective
    #mod.ext[:objective] = @objective(mod, Min,sum((CAPEX*cap[y] + OPEX*g[y]) for y in Y))
    mod.ext[:objective] = @objective(mod, Min, sum(MAC*a[y]^2*A[y] for y in Y))
 
    # Define constraint
    mod.ext[:constraints] = Dict()
    #mod.ext[:constraints][:buycons] = @constraint(mod,[y=Y], buy[y] <= ETS["S"][y])
    mod.ext[:constraints][:auctioncons] = @constraint(mod, [y=Y], sum(ETS["S"][1:y]) >= sum(E[1:y]) - sum(a[1:y]))
    mod.ext[:constraints][:nonnegemiss] = @constraint(mod, [y=Y], 0 <= E[y] - a[y] )
    #mod.ext[:constraints][:capacitycons] = @constraint(mod, [y=Y], sum(cap[1:y]) >= g[y])
    return mod
end
define_sets_parameters!(mod,data)

optimize!(mod)
@show objective_value(mod);


#buy_opt = convert(Array, JuMP.value.(mod.ext[:variables][:buy]))
abate_opt = convert(Array, JuMP.value.(mod.ext[:variables][:a]))
bank_opt = convert(Array, JuMP.value.(mod.ext[:expressions][:bank]))
emission_opt = mod.ext[:parameters][:e][:,1] - abate_opt 
supply = ETS["S"][:,1]
sol = DataFrame(Y=2019:2063,Supply=supply,Abate=abate_opt, Emmit=emission_opt,Bank=bank_opt)

plot(Matrix(sol)[:,2:end],labels=permutedims(names(sol[:,2:end])))

prices =  [dual(mod.ext[:constraints][:auctioncons][i]) for i in 1:45]