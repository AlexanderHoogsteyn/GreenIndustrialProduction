using CSV, DataFrames, Glob

# Get current script directory
script_dir = @__DIR__

# Define required cost columns
const AGENT_COLUMNS = Dict(
    :bf_bof => "BF-BOF_cost_mean",
    :bf_bof_cc => "BF-BOF-CC_cost_mean",
    :hdri_eaf => "HDRI-EAF_cost_mean",
    :ng_dri_cc => "NG-DRI-CC_cost_mean",
    :fringe => "fringe_π_MAC_mean"
)

# Prepare output DataFrame
output_columns = [
    :sensitivity, :scenario,
    :total_bf_bof, :total_bf_bof_cc, :total_hdri_eaf, :total_ng_dri_cc, :total_fringe,
    :steel_industry_total, :system_total
]
output_df = DataFrame(; [col => [] for col in output_columns]...)

# Find all CSV files in script directory
files = glob("scenario_*.csv", script_dir)

if isempty(files)
    @warn "No scenario files found in directory: $script_dir"
else
    println("Found $(length(files)) scenario files in directory: $script_dir")
end

for file in files
    # Extract filename without path
    filename = basename(file)
    println("Processing: $filename")
    
    # Flexible filename pattern matching
    pattern_match = match(r"scenario_(\d+)_(\d+)\.csv$"i, filename)
    
    if pattern_match === nothing
        # Try alternative pattern if first fails
        pattern_match = match(r"scenario_(\d+)[^\d]*(\d+)\.csv$"i, filename)
    end
    
    if pattern_match === nothing
        @warn "  Skipping - Filename pattern not recognized"
        continue
    end
    
    # Parse sensitivity and scenario numbers
    sensitivity = tryparse(Int, pattern_match[1])
    scenario = tryparse(Int, pattern_match[2])
    
    if sensitivity === nothing || scenario === nothing
        @warn "  Skipping - Invalid numbers in filename"
        continue
    end

    try
        # Read CSV file
        df = CSV.read(file, DataFrame)
        
        # Check for required columns
        missing_cols = String[]
        for col in values(AGENT_COLUMNS)
            if !(col in names(df))
                push!(missing_cols, col)
            end
        end
        
        if !isempty(missing_cols)
            @warn "  Skipping - Missing columns: $(join(missing_cols, ", "))"
            continue
        end
        
        # Calculate total costs
        totals = Dict{Symbol, Float64}()
        for (agent, col) in AGENT_COLUMNS
            totals[agent] = sum(skipmissing(df[!, col]))
        end
        
        # Calculate aggregates
        steel_agents = [:bf_bof, :bf_bof_cc, :hdri_eaf, :ng_dri_cc]
        steel_total = sum(totals[agent] for agent in steel_agents)
        system_total = steel_total + totals[:fringe]
        
        # Add to results
        push!(output_df, (
            sensitivity, scenario,
            totals[:bf_bof], totals[:bf_bof_cc], totals[:hdri_eaf], totals[:ng_dri_cc], totals[:fringe],
            steel_total, system_total
        ))
        println("  Processed successfully - Sensitivity: $sensitivity, Scenario: $scenario")
        
    catch e
        @warn "  Error processing file: $(sprint(showerror, e))"
    end
end

# Save results
if nrow(output_df) > 0
    # Round all cost columns to 1 decimal place
    cost_cols = setdiff(propertynames(output_df), [:sensitivity, :scenario])
    for col in cost_cols
        output_df[!, col] = round.(output_df[!, col], digits=1)
    end
    
    output_path = joinpath(script_dir, "aggregated_costs_results.csv")
    CSV.write(output_path, output_df)
    println("\nSuccess! Processed $(nrow(output_df)) files")
    println("Results saved to: $output_path")
    
    # Show preview
    println("\nOutput preview (first 5 rows):")
    if nrow(output_df) > 5
        display(first(output_df, 5))
    else
        display(output_df)
    end
else
    @warn "\nNo valid files processed. Output file not created."
end