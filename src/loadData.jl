function get_solution(agents::Dict, results::Dict)
    frac_digit = 4
    λ_ets = round.(results["λ"]["ETS"][end], digits=frac_digit)
    λ_product = round.(results["λ"]["product"][end], digits=frac_digit)
    
    # Initialize the DataFrame with the year column
    sol = DataFrame(Y=2024:(2023 + size(λ_ets, 1)))
    
    # Add the ETS and product matrix as separate columns for each scenario
    for col in 1:size(λ_ets, 2)
        sol[!, Symbol("λ_ets_scenario_" * string(col))] = λ_ets[:, col]
        sol[!, Symbol("λ_product_scenario_" * string(col))] = λ_product[:, col]
    end

    for (key, agent) in agents
        for variable in keys(agent.ext[:variables])
            values = round.(convert(Array, JuMP.value.(agent.ext[:variables][variable])), digits=frac_digit)
            variable_name = Symbol(string(key) * "_" * string(variable))
            
            if ndims(values) == 1
                # If it's a vector, assign it directly to a column
                sol[!, variable_name] = values
            elseif ndims(values) == 2
                # If it's a matrix, store each column as a separate variable
                for col in 1:size(values, 2)
                    col_name = Symbol(string(variable_name) * "_scenario_" * string(col))
                    sol[!, col_name] = values[:, col]
                end
            else
                error("Unsupported variable dimension: ", ndims(values))
            end
        end
        
        for expression in keys(agent.ext[:expressions])
            values = round.(convert(Array, JuMP.value.(agent.ext[:expressions][expression])), digits=frac_digit)
            expression_name = Symbol(string(key) * "_" * string(expression))
            
            if ndims(values) == 1
                sol[!, expression_name] = values
            elseif ndims(values) == 2
                for col in 1:size(values, 2)
                    col_name = Symbol(string(expression_name) * "_scenario_" * string(col))
                    sol[!, col_name] = values[:, col]
                end
            else
                error("Unsupported expression dimension: ", ndims(values))
            end
        end
    end
    return sol
end

function get_solution_summarized(agents::Dict, results::Dict)
    frac_digit = 4
    λ_ets = round.(results["λ"]["ETS"][end], digits=frac_digit)
    λ_product = round.(results["λ"]["product"][end], digits=frac_digit)
    
    # Initialize the DataFrame with the year column
    sol = DataFrame(Y=2024:(2023 + size(λ_ets, 1)))
    
    # Add summarized metrics for λ_ets and λ_product in the new order: min, Q1, mean, Q3, max
    sol[!, :λ_ets_min] = minimum(λ_ets, dims=2)[:, 1]
    sol[!, :λ_ets_Q1] = [quantile(row, 0.25) for row in eachrow(λ_ets)]
    sol[!, :λ_ets_mean] = median(λ_ets, dims=2)[:, 1]
    sol[!, :λ_ets_Q3] = [quantile(row, 0.75) for row in eachrow(λ_ets)]
    sol[!, :λ_ets_max] = maximum(λ_ets, dims=2)[:, 1]
    
    sol[!, :λ_product_min] = minimum(λ_product, dims=2)[:, 1]
    sol[!, :λ_product_Q1] = [quantile(row, 0.25) for row in eachrow(λ_product)]
    sol[!, :λ_product_mean] = median(λ_product, dims=2)[:, 1]
    sol[!, :λ_product_Q3] = [quantile(row, 0.75) for row in eachrow(λ_product)]
    sol[!, :λ_product_max] = maximum(λ_product, dims=2)[:, 1]

    for (key, agent) in agents
        for variable in keys(agent.ext[:variables])
            values = round.(convert(Array, JuMP.value.(agent.ext[:variables][variable])), digits=frac_digit)
            variable_name = Symbol(string(key) * "_" * string(variable))
            
            if ndims(values) == 1
                # If it's a vector, assign it directly to a column
                sol[!, variable_name] = values
            elseif ndims(values) == 2
                # If it's a matrix, compute summary statistics in the new order: min, Q1, mean, Q3, max
                sol[!, Symbol(string(variable_name) * "_min")] = minimum(values, dims=2)[:, 1]
                sol[!, Symbol(string(variable_name) * "_Q1")] = [quantile(row, 0.25) for row in eachrow(values)]
                sol[!, Symbol(string(variable_name) * "_mean")] = median(values, dims=2)[:, 1]
                sol[!, Symbol(string(variable_name) * "_Q3")] = [quantile(row, 0.75) for row in eachrow(values)]
                sol[!, Symbol(string(variable_name) * "_max")] = maximum(values, dims=2)[:, 1]
            else
                error("Unsupported variable dimension: ", ndims(values))
            end
        end
        
        for expression in keys(agent.ext[:expressions])
            values = round.(convert(Array, JuMP.value.(agent.ext[:expressions][expression])), digits=frac_digit)
            expression_name = Symbol(string(key) * "_" * string(expression))
            
            if ndims(values) == 1
                sol[!, expression_name] = values
            elseif ndims(values) == 2
                # Compute summary statistics for expressions in the new order: min, Q1, mean, Q3, max
                sol[!, Symbol(string(expression_name) * "_min")] = minimum(values, dims=2)[:, 1]
                sol[!, Symbol(string(expression_name) * "_Q1")] = [quantile(row, 0.25) for row in eachrow(values)]
                sol[!, Symbol(string(expression_name) * "_mean")] = median(values, dims=2)[:, 1]
                sol[!, Symbol(string(expression_name) * "_Q3")] = [quantile(row, 0.75) for row in eachrow(values)]
                sol[!, Symbol(string(expression_name) * "_max")] = maximum(values, dims=2)[:, 1]
            else
                error("Unsupported expression dimension: ", ndims(values))
            end
        end
    end
    return sol
end

function define_results(data::Dict,agents::Dict) 
    results = Dict()
    ADMM = Dict()
    results["g"] = Dict()
    results["b"] = Dict()
    results["e"] = Dict()
    results["e"] = Dict()
    results["g_τ"] = Dict()
    results["π_MAC"] = Dict()

    # TO DO: incorporate setting the sector better
    sector = "steelmaking"


    for (agent,model) in agents
        results["b"][agent] = CircularBuffer{Array{Float64,1}}(data["CircularBufferSize"])  
        push!(results["b"][agent],zeros(data["nyears"]))
        results["e"][agent] = CircularBuffer{Array{Float64,1}}(data["CircularBufferSize"])  
        push!(results["e"][agent],zeros(data["nyears"]))
        results["g"][agent] = CircularBuffer{Array{Float64,1}}(data["CircularBufferSize"]) 
        push!(results["g"][agent],zeros(data["nyears"]))
        if is_myopic(model)
            results["g_τ"][agent] = CircularBuffer{Array{Float64,2}}(data["CircularBufferSize"]) 
            push!(results["g_τ"][agent],zeros(data["nyears"],data["nyears"]))
        end
    end

    results["e"]["fringe"] = CircularBuffer{Array{Float64,1}}(data["CircularBufferSize"])  
    push!(results["e"]["fringe"],zeros(data["nyears"]))

    results["s"] = data["S"][:]
    results["D"] = data["D"][:]
    
    results["λ"] = Dict()
    results["λ"]["ETS"] = CircularBuffer{Array{Float64,1}}(data["CircularBufferSize"]) 
    push!(results["λ"]["ETS"],ones(data["nyears"]))
    results["λ"]["product"] = CircularBuffer{Array{Float64,1}}(data["CircularBufferSize"]) 
    push!(results["λ"]["product"],ones(data["nyears"]))


    ADMM["Imbalances"] = Dict()
    ADMM["Imbalances"]["ETS"] = CircularBuffer{Array{Float64,1}}(data["CircularBufferSize"])  
    push!(ADMM["Imbalances"]["ETS"],zeros(data["nyears"]))
    ADMM["Imbalances"]["product"] = CircularBuffer{Array{Float64,1}}(data["CircularBufferSize"])
    push!(ADMM["Imbalances"]["product"],zeros(data["nyears"]))

    ADMM["Residuals"] = Dict()
    ADMM["Residuals"]["Primal"] = Dict()
    ADMM["Residuals"]["Primal"]["ETS"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Primal"]["ETS"],0)
    ADMM["Residuals"]["Primal"]["product"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Primal"]["product"],0)

    ADMM["Residuals"]["Dual"] = Dict()
    ADMM["Residuals"]["Dual"]["ETS"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Dual"]["ETS"],0)
    ADMM["Residuals"]["Dual"]["product"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Dual"]["product"],0)

    ADMM["Tolerance"] = Dict()
    ADMM["Tolerance"]["ETS"] = data["epsilon"]/100*maximum(data["CAP_2016"])*sqrt(data["nyears"])
    ADMM["Tolerance"]["product"] = data["epsilon"]/100*maximum(data["demand"][sector])*sqrt(data["nyears"])

    ADMM["ρ"] = Dict()
    ADMM["ρ"]["ETS"] = CircularBuffer{Float64}(data["CircularBufferSize"]) 
    push!(ADMM["ρ"]["ETS"],data["rho_ETS"])
    ADMM["ρ"]["product"] = CircularBuffer{Float64}(data["CircularBufferSize"]) 
    push!(ADMM["ρ"]["product"],data["rho_product"])

    ADMM["n_iter"] = 1 
    ADMM["walltime"] = 0
    ADMM[:start] = 1
    ADMM[:end] = data["nyears"]
    ADMM[:mask] = ones(data["nyears"])

    return results, ADMM
end

function define_results_stochastic(data::Dict,agents::Dict)
    results, ADMM = define_results(data,agents)

    for (agent,model) in agents
        results["b"][agent] = CircularBuffer{Array{Float64,2}}(data["CircularBufferSize"])  
        push!(results["b"][agent],zeros(data["nyears"],data["nsamples"]))
        results["e"][agent] = CircularBuffer{Array{Float64,2}}(data["CircularBufferSize"])  
        push!(results["e"][agent],zeros(data["nyears"],data["nsamples"]))
        results["g"][agent] = CircularBuffer{Array{Float64,2}}(data["CircularBufferSize"]) 
        push!(results["g"][agent],zeros(data["nyears"],data["nsamples"]))
        if is_myopic(model)
            results["g_τ"][agent] = CircularBuffer{Array{Float64,3}}(data["CircularBufferSize"]) 
            push!(results["g_τ"][agent],zeros(data["nyears"],data["nyears"],data["nsamples"]))
        end
    end

    results["e"]["fringe"] = CircularBuffer{Array{Float64,2}}(data["CircularBufferSize"])  
    push!(results["e"]["fringe"],zeros(data["nyears"],data["nsamples"]))


    results["λ"]["ETS"] = CircularBuffer{Array{Float64,2}}(data["CircularBufferSize"]) 
    push!(results["λ"]["ETS"],ones(data["nyears"],data["nsamples"]))
    results["λ"]["product"] = CircularBuffer{Array{Float64,2}}(data["CircularBufferSize"]) 
    push!(results["λ"]["product"],ones(data["nyears"],data["nsamples"]))

    results["π_MAC"]["fringe"] = CircularBuffer{Array{Float64,2}}(data["CircularBufferSize"])  
    push!(results["π_MAC"]["fringe"],zeros(data["nyears"],data["nsamples"]))

    ADMM["Imbalances"]["ETS"] = CircularBuffer{Array{Float64,2}}(data["CircularBufferSize"])  
    push!(ADMM["Imbalances"]["ETS"],ones(data["nyears"],data["nsamples"]))
    ADMM["Imbalances"]["product"] = CircularBuffer{Array{Float64,2}}(data["CircularBufferSize"])
    push!(ADMM["Imbalances"]["product"],zeros(data["nyears"],data["nsamples"]))

    return results, ADMM
end

function define_results_hot_start!(data::Dict,results::Dict,ADMM::Dict)

    sector = "steelmaking"

    # Convert Arrays to circular buffers
    for (r,arr) in results
        for (m,arr2) in arr
            cb = CircularBuffer(data["CircularBufferSize"])
            append!(cb,arr2)
            results[r][m] = cb
        end
    end
    for (r,arr) in ADMM["ρ"]
        cb = CircularBuffer(data["CircularBufferSize"])
        append!(cb, arr)
        ADMM["ρ"][r] = cb
    end
    for (r,arr) in ADMM["Residuals"]["Primal"]
        cb = CircularBuffer(data["CircularBufferSize"])
        append!(cb, arr)
        ADMM["Residuals"]["Primal"][r] = cb
    end
    for (r,arr) in ADMM["Residuals"]["Dual"]
        cb = CircularBuffer(data["CircularBufferSize"])
        append!(cb, arr)
        ADMM["Residuals"]["Dual"][r] = cb
    end
    for (r,arr) in ADMM["Imbalances"]
        cb = CircularBuffer(data["CircularBufferSize"])
        append!(cb, arr)
        ADMM["Imbalances"][r] = cb
    end

    ADMM["Tolerance"] = Dict()
    ADMM["Tolerance"]["ETS"] = data["epsilon"]/100*maximum(data["CAP_2016"])*sqrt(data["nyears"])
    ADMM["Tolerance"]["good"] = data["epsilon"]/100*maximum(data["demand"][sector])*sqrt(data["nyears"])

    ADMM["n_iter"] = 1 
    ADMM["walltime"] = 0

    return results, ADMM
end