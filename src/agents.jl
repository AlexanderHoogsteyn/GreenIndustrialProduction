# Parameters/variables ETS 
function define_ETS_parameters!(data::Dict)
    # LRF 2017 - 2021, no Brexit, expansion of scope (aviation, maritime) or Green Deal
    data["LRF"] = zeros(data["nyears"],1);
    data["LRF"][1:4] = data["LRF_2024"]/100*ones(4,1);                            # 2021-2030
    data["LRF"][5:end] = data["LRF_2028"]/100*ones(data["nyears"]-4,1);           # 2030 - end ETS
       
    # CAP
    data["CAP"] = zeros(data["nyears"],1);
    for y =1:data["nyears"]
        data["CAP"][y]= maximum([data["CAP_2024"]*(1-sum(data["LRF"][1:y])) 0])
    end
    data["CAP"][1] = data["CAP_2024"] #+data["EX_2016"] # real supply in 2017 and surplus in market at end 2016

    # Supply 
    data["S"] = copy(data["CAP"])

    return data
end

function define_sector_parameters!(data::Dict,route::String)
    # Demand good
    data["D"] = ones(data["nyears"],1)*data["demand"][route];

    return data
end


# Stochastic optimization parameters
function define_stoch_parameters!(stoch::Dict,data::Dict)
    nsamples = data["nsamples"]
    avg = data["MAC"]
    std = data["std"]

    if data["use_cvar"]
        d = Normal(avg,std)
        nselect = floor(Int,data["cvar"]*nsamples)
        stoch["MAC"] = sort(rand(d,nsamples), rev=true)[1:nselect]
        stoch["nsamples"] = size(SO["MAC"])[1]
    else
        stoch["MAC"] = avg
        stoch["nsamples"] = 1
    end
    return stoch
end

# Model Representation of the agents
function build_agent!(agent::Model, data::Dict)
    # Build structure of agents which is shared among all types
    agent.ext[:sets] = Dict()
    Y = agent.ext[:sets][:Y] = 1:data["nyears"]
    agent.ext[:sets][:S] = 1:1

    # define parameters 
    agent.ext[:parameters] = Dict()
    agent.ext[:parameters][:λ_ets] = zeros(data["nyears"],1)
    agent.ext[:parameters][:λ_elec] = zeros(data["nyears"],1)
    agent.ext[:parameters][:λ_product] = zeros(data["nyears"],1)
    agent.ext[:parameters][:A] = ones(data["nyears"],1)
    agent.ext[:parameters][:i] = ones(data["nyears"],1)
    agent.ext[:parameters][:r_equity] = ones(data["nyears"],1)
    agent.ext[:parameters][:r_debt] = ones(data["nyears"],1)
    for y in 1:data["nyears"]
        agent.ext[:parameters][:A][y] = 1/(1+data["discount_rate"])^(y);
        agent.ext[:parameters][:r_debt][y] = 1/(1+data["debt_rate"])^(y);
        agent.ext[:parameters][:r_equity][y] = 1/(1+data["equity_rate"])^(y);
        agent.ext[:parameters][:i][y] = (1+data["inflation"])^(y);
    end

    # define expressions
    agent.ext[:expressions] = Dict()

    # define variables
    agent.ext[:variables] = Dict()
    b = agent.ext[:variables][:b] =  @variable(agent, [y=Y], base_name="allowances bougth")

    # Define constraint
    agent.ext[:constraints] = Dict()

    # Objective
    agent.ext[:objective] = @objective(agent, Min, 0)
end

function build_agent!(agent::Model,data::Dict,stoch::Dict)
    build_agent!(agent,data)
    agent.ext[:sets][:S] = 1:stoch["nsamples"]
end
