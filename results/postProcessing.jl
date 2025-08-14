using CSV
using DataFrames
using Statistics

# Load the scenarios CSV file
scenarios_filepath = joinpath(@__DIR__, "../data/scenarios.csv")
scenarios_df = CSV.read(scenarios_filepath, DataFrame)

# Function to process a single file and extract relevant data
function extract_ets_prices_and_mac(filepath::String)
    # Read the CSV file into a DataFrame
    df = CSV.read(filepath, DataFrame)
    
    # Filter for years 2025 and 2035
    filtered_df = df[(df.Y .== 2025) .| (df.Y .== 2035), :]
    
    # Extract ETS price statistics (mean, max, etc.) for λ
    ets_prices = filtered_df[!, r"λ.*"]
    ets_prices[!, :year] = Vector(filtered_df[!, :Y])  # Add the 'year' column explicitly as a vector

    # Calculate MAC statistics across all years
    total_mac_mean = mean(df[!, :fringe_π_MAC_mean])
    total_mac_min = minimum(df[!, :fringe_π_MAC_min])
    total_mac_max = maximum(df[!, :fringe_π_MAC_max])
    total_mac_q1 = quantile(df[!, :fringe_π_MAC_mean], 0.25)
    total_mac_q3 = quantile(df[!, :fringe_π_MAC_mean], 0.75)
    
    # Add these statistics as columns (repeated for all rows)
    ets_prices[!, :total_mac_mean] = fill(total_mac_mean, nrow(ets_prices))
    ets_prices[!, :total_mac_min] = fill(total_mac_min, nrow(ets_prices))
    ets_prices[!, :total_mac_max] = fill(total_mac_max, nrow(ets_prices))
    ets_prices[!, :total_mac_q1] = fill(total_mac_q1, nrow(ets_prices))
    ets_prices[!, :total_mac_q3] = fill(total_mac_q3, nrow(ets_prices))

    return ets_prices, total_mac_mean
end

##########################
#  Variable foresight   #
##########################

# File pattern and scenario numbers
file_pattern = "scenario_2_"
scenario_range = 12:21

# Initialise a DataFrame to hold all results
results = DataFrame()
scenario_mac_totals = Dict{Int, Float64}()

# Loop through each scenario file, process, and append results
for scenario in scenario_range
    filepath = joinpath(@__DIR__, "$file_pattern$scenario.csv")
    if isfile(filepath)
        ets_prices, total_mac_mean = extract_ets_prices_and_mac(filepath)
        
        # Store the total MAC for this scenario
        scenario_mac_totals[scenario] = total_mac_mean
        
        # Add scenario as a column
        ets_prices[!, :scenario] = fill(scenario, nrow(ets_prices))  # Create a vector of the same length
        
        # Add horizon_ets column based on the scenarios_df
        horizon_ets_value = scenarios_df[scenarios_df.scenario_id .== scenario, :horizon_ets][1]
        ets_prices[!, :horizon_ets] = fill(horizon_ets_value, nrow(ets_prices))  # Add the value to all rows for this scenario
        
        # Append to results DataFrame
        append!(results, ets_prices)
    end
end

# Calculate percentage changes relative to Scenario 20
baseline_mac = 38813.682552 #scenario_mac_totals[12]
results[!, :percentage_change_mac] = [((scenario_mac_totals[row[:scenario]] - baseline_mac) / baseline_mac) * 100 for row in eachrow(results)]

# Sort the DataFrame so rows for 2025 appear first, followed by 2035
results = sort(results, :year)
results[!, :index] = 1:nrow(results)

# Reshape the DataFrame to long format
results_long = stack(results, Not([:scenario, :year, :horizon_ets, :total_mac_mean, :total_mac_min, :total_mac_max, :total_mac_q1, :total_mac_q3, :percentage_change_mac]))

# Pivot the long DataFrame
results_pivoted = unstack(results_long, [:scenario, :year, :horizon_ets, :total_mac_mean, :total_mac_min, :total_mac_max, :total_mac_q1, :total_mac_q3, :percentage_change_mac], :variable, :value)

# Save the final DataFrame to a CSV for review
output_path = joinpath(@__DIR__, "ets_prices_summary_rolling_horizon.csv")
CSV.write(output_path, results_pivoted)

println("ETS prices summary saved to $output_path")

###############################
#  Variable liquidity factor  #
###############################

# File pattern and scenario numbers
scenario_range = 2:11

# Initialise a DataFrame to hold all results
results = DataFrame()
scenario_mac_totals = Dict{Int, Float64}()

# Loop through each scenario file, process, and append results
for scenario in scenario_range
    filepath = joinpath(@__DIR__, "$file_pattern$scenario.csv")
    if isfile(filepath)
        ets_prices, total_mac_mean = extract_ets_prices_and_mac(filepath)
        
        # Store the total MAC for this scenario
        scenario_mac_totals[scenario] = total_mac_mean
        
        # Add scenario as a column
        ets_prices[!, :scenario] = fill(scenario, nrow(ets_prices))  # Create a vector of the same length
        
        # Add liquidity_factor column based on the scenarios_df
        liquidity_factor_value = scenarios_df[scenarios_df.scenario_id .== scenario, :liquidity_factor][1]
        ets_prices[!, :liquidity_factor] = fill(liquidity_factor_value, nrow(ets_prices))  # Add the value to all rows for this scenario

        # Append to results DataFrame
        append!(results, ets_prices)
    end
end

# Calculate percentage changes relative to Scenario 2
baseline_mac = 38813.682552 #scenario_mac_totals[2]
results[!, :percentage_change_mac] = [((scenario_mac_totals[row[:scenario]] - baseline_mac) / baseline_mac) * 100 for row in eachrow(results)]

# Sort the DataFrame so rows for 2025 appear first, followed by 2035
results = sort(results, :year)
results[!, :index] = 1:nrow(results)

# Reshape the DataFrame to long format
results_long = stack(results, Not([:scenario, :year, :liquidity_factor, :total_mac_mean, :total_mac_min, :total_mac_max, :total_mac_q1, :total_mac_q3, :percentage_change_mac]))

# Pivot the long DataFrame
results_pivoted = unstack(results_long, [:scenario, :year, :liquidity_factor, :total_mac_mean, :total_mac_min, :total_mac_max, :total_mac_q1, :total_mac_q3, :percentage_change_mac], :variable, :value)

# Save the final DataFrame to a CSV for review
output_path = joinpath(@__DIR__, "ets_prices_summary_liquidity_constraint.csv")
CSV.write(output_path, results_pivoted)

println("ETS prices summary saved to $output_path")


###############################
#  Variable risk premium      #
###############################

# # File pattern and scenario numbers
# scenario_range = 22:31

# # Initialise a DataFrame to hold all results
# results = DataFrame()
# scenario_mac_totals = Dict{Int, Float64}()

# # Loop through each scenario file, process, and append results
# for scenario in scenario_range
#     filepath = joinpath(@__DIR__, "$file_pattern$scenario.csv")
#     if isfile(filepath)
#         ets_prices, total_mac_mean = extract_ets_prices_and_mac(filepath)
        
#         # Store the total MAC for this scenario
#         scenario_mac_totals[scenario] = total_mac_mean
        
#         # Add scenario as a column
#         ets_prices[!, :scenario] = fill(scenario, nrow(ets_prices))  # Create a vector of the same length
        
#         # Add liquidity_factor column based on the scenarios_df
#         liquidity_factor_value = scenarios_df[scenarios_df.scenario_id .== scenario, :liquidity_factor][1]
#         ets_prices[!, :liquidity_factor] = fill(liquidity_factor_value, nrow(ets_prices))  # Add the value to all rows for this scenario

#         # Append to results DataFrame
#         append!(results, ets_prices)
#     end
# end

# # Calculate percentage changes relative to Scenario 22
# baseline_mac = 38813.682552 #scenario_mac_totals[22]
# results[!, :percentage_change_mac] = [((scenario_mac_totals[row[:scenario]] - baseline_mac) / baseline_mac) * 100 for row in eachrow(results)]

# # Sort the DataFrame so rows for 2025 appear first, followed by 2035
# results = sort(results, :Y)
# results[!, :index] = 1:nrow(results)

# # Reshape the DataFrame to long format
# results_long = stack(results, Not([:scenario, :year, :liquidity_factor, :total_mac_mean, :total_mac_min, :total_mac_max, :total_mac_q1, :total_mac_q3, :percentage_change_mac]))

# # Pivot the long DataFrame
# results_pivoted = unstack(results_long, [:scenario, :year, :liquidity_factor, :total_mac_mean, :total_mac_min, :total_mac_max, :total_mac_q1, :total_mac_q3, :percentage_change_mac], :variable, :value)


# output_path = joinpath(@__DIR__, "ets_prices_summary_discounting.csv")
# CSV.write(output_path, results_pivoted)

# println("ETS prices summary saved to $output_path")
