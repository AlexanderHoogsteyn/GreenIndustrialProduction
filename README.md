# GreenIndustrialProduction
This Julia project is designed to perform agent-based optimization modeling capturing the interaction of different steelmaking decarbonization routes using Gurobi as the optimization solver. The decarbonization of the steel sector is governed by a rising carbon price, which we consider endogenously by modeling the EU emission trading scheme (ETS). The model integrates several components to simulate different agents' behavior, such as myopic foresight, leveraging the ADMM (Alternating Direction Method of Multipliers) algorithm to solve the model. It can be easily adapted to model a different sector than steelmaking.

## Model Overview
This project implements a partial equilibrium model to simulate the interaction of various agents in the EU Emission Trading Scheme (ETS). The model captures the behavior of producers, financial traders, and a competitive fringe, considering factors such as policy risk, myopic behavior, and financial constraints. The model is designed to evaluate the impact of these factors on carbon pricing and decarbonization pathways in the steelmaking sector.

Key features of the model include:
- **Multi-agent Modeling**: Producers, traders, and fringe agents interact within the ETS framework.
- **Stochastic Optimization**: Incorporates uncertainty in key parameters such as emissions and abatement costs.
- **Rolling Horizon Optimization**: Models agents with limited foresight to simulate myopic behavior.
- **Financial Constraints**: Captures the impact of liquidity constraints on agent decisions.
- **Policy Risk**: Evaluates the effect of policy uncertainty on carbon pricing and investment decisions.

### Key Agents

The model includes the following types of agents:

1. **Producers**:
   - Represent steelmaking firms with specific production routes.
   - Optimize production and investment decisions under financial constraints and policy uncertainty.
   - Can operate under deterministic or stochastic settings.

2. **Financial Traders**:
   - Trade allowances in the ETS to maximize profits.
   - Can bank allowances for future use or sale.
   - Operate under liquidity constraints in some scenarios.

3. **Competitive Fringe**:
   - Represents compliance actors with actual emissions in sectors not explicitly modeled.
   - Ensures the ETS market clears by balancing supply and demand.

4. **Agents**:
   - General agent type that contains functions common to all agents.

## Requirements
### Julia Version
- Julia 1.10 or higher

### Julia Packages
The script utilizes the following Julia packages:
- `JuMP`: A modeling language for optimization problems.
- `Gurobi`: An optimization solver (requires a valid license).
- `DataFrames`: For handling data in tabular form.
- `CSV`: For reading and writing CSV files.
- `YAML`: For reading and writing configuration files in YAML format.
- `DataStructures`: Provides data structure utilities.
- `ProgressBars`: Adds progress bars for long-running operations.
- `Printf`: For formatted output.
- `JLD2`: For saving and loading Julia data.
- `Base.Threads`: For threading and parallel execution.

### Solver
- **Gurobi Solver**: Requires installation and a valid Gurobi license. The script also configures the Gurobi environment with specific parameters such as suppressing output (`OutputFlag = 0`) and limiting computation time (`TimeLimit = 300 seconds`).
- **(Alternative) HiGHS**: Open-source alternative to Gurobi. The script could be easily adapted to make use of any alternative solvers such as this one - as long as it is capable of solving quadratic optimization problems. E.g.: [HiGHS](https://highs.dev/)

## Files and Directory Structure

- `agents.jl`: Defines the common structure and parameters for all agent types.
- `loadData.jl`: Loads and preprocesses data, including scenarios and assumptions.
- `backbone.jl`: Core functions and utilities for the model backbone.
- `ADMM.jl`: Implements the ADMM optimization algorithm for solving the model.
- `producer.jl`: Defines producer agents, including their production capacity, emissions, and financial constraints.
- `fringe.jl`: Models the competitive fringe, representing compliance actors with actual emissions.
- `traders.jl`: Defines financial trading agents that can bank allowances and interact with the ETS.
- `postProcessing.jl`: Processes and analyzes model outputs, including ETS price summaries and MAC statistics.
- `MAIN.jl`: The main script to run the model for multiple scenarios and save results.
- `demos/`: Contains example scripts demonstrating the model's usage.
- `results/`: Directory where model outputs are saved, including scenario results and detailed outputs.
- `data/`: Contains input data files, including:
  - `scenarios.csv`: Defines the scenarios to be modeled.
  - `sensetivities.csv`: Contains sensitivity analysis parameters.
  - `assumptions.yaml`: Stores common assumptions and parameters.

### Data Files
- `../data/scenarios.csv`: Defines the different scenarios to be modeled. Each row represents a scenario, and columns represent the parameters for each scenario.
- `../data/assumptions.yaml`: Contains assumptions and parameters common to all scenarios, including commodity prices.

## Usage

### Running the Script
The script in `MAIN.JL` can be used to run the scenarios you defined. The variant `MAIN_single.JL` runs a single scenario, which you can optionally parse as an argument in the command line as follows:
```bash
julia MAIN_single.jl --scenario 1
```

### Customization

- **Scenario Configuration**: Scenarios can be customized by modifying the `scenarios.csv` file. Each row represents a scenario, and columns represent the parameters for each scenario.
- **Common Parameters**: Modify the common parameters, including commodity prices, in the `assumptions.yaml` file.

## Output

- The script will generate CSV files containing the solution for each scenario. These files will be saved in the `results/` directory. Results will be saved as `scenario_X.csv` where `X` is the scenario number.

## License

This project is licensed under the MIT License. Please refer to the `LICENSE` file for more details.

## Contact

For questions or contributions, please contact: alexander.hoogsteyn@kuleuven.be

## Citing This Work

If you use this model in your research, please cite the following working paper:

Hoogsteyn A., Meus J., Bruninx K., Delarue E. "Barriers to efficient carbon pricing: policy risk, myopic behavior, and financial constraints." 2025. Working paper.

