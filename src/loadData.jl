using CSV
using DataFrames
using YAML
using Plots
using Statistics
using Distributions, Random
using StatsPlots

function get_solution(agents::Dict,results::Dict)
    frac_digit = 4
    λ_ets = round.(results["λ"]["ETS"][end],digits=frac_digit)
    λ_product = round.(results["λ"]["product"][end],digits=frac_digit)
    sol = DataFrame(Y=2024:(2023+length(λ_ets)),λ_ets=λ_ets, λ_product=λ_product)

    for (key,agent) in agents
        for variable in keys(agent.ext[:variables])
            values = round.(convert(Array, JuMP.value.(agent.ext[:variables][variable])),digits=frac_digit)
            variable_name = Symbol(string(key) * "_" * string(variable))
            sol[!,variable_name] = values 
        end
        for expression in keys(agent.ext[:expressions])
            values = round.(convert(Array, JuMP.value.(agent.ext[:expressions][expression])),digits=frac_digit)
            variable_name = Symbol(string(key) * "_" * string(expression))
            sol[!,variable_name] = values 
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
    results["g_τ"] = Dict()

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

    results["s"] = CircularBuffer{Array{Float64,1}}(data["CircularBufferSize"])  
    push!(results["s"],zeros(data["nyears"]))
    results["D"] = CircularBuffer{Array{Float64,1}}(data["CircularBufferSize"])  
    push!(results["s"],zeros(data["nyears"]))

    results["λ"] = Dict()
    results["λ"]["ETS"] = CircularBuffer{Array{Float64,1}}(data["CircularBufferSize"]) 
    push!(results["λ"]["ETS"],zeros(data["nyears"]))
    results["λ"]["product"] = CircularBuffer{Array{Float64,1}}(data["CircularBufferSize"]) 
    push!(results["λ"]["product"],zeros(data["nyears"]))

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