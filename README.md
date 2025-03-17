# GreenIndustrialProduction
This Julia project is designed to perform agent-based optimization modeling capturing the interaction of different steelmaking decarbonization routes using Gurobi as the optimization solver. The decarbonization of the steel sector is governed by a rising carbon price, which we consider endogenously by modeling the EU emission trading scheme (ETS). The model integrates several components to simulate different agents' behavior, such as myopic foresight, leveraging the ADMM (Alternating Direction Method of Multipliers) algorithm to solve the model. It can be easily adapted to model a different sector than steelmaking.

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

- `agents.jl`: Defines the different agent types in the model.
- `loadData.jl`: Loads and preprocesses data for the model.
- `backbone.jl`: Core functions and utilities for the model backbone.
- `ADMM.jl`: Implements the ADMM optimization algorithm.
- `producer.jl`: Defines the producer agents for the sector being modeled.
- `fringe.jl`: Defines the competitive fringe in the emission trading scheme that represents compliance actors with actual emissions in sectors that are not explicitly modelled
- `traders.jl`: Defines financial trading agents that can bank allowances

### Data Files
- `../data/scenarios.csv`: Defines the different scenarios to be modeled. Each row represents a scenario, and columns represent the parameters for each scenario.
- `../data/assumptions.yaml`: Contains assumptions and parameters common to all scenarios, including commodity prices.

## Usage

### Running the Script
Examples of how to use the code are provided in the `demos` folder.

### Customization

- **Scenario Configuration**: Scenarios can be customized by modifying the `scenarios.csv` file. Each row represents a scenario, and columns represent the parameters for each scenario.
- **Common Parameters**: Modify the common parameters, including commodity prices, in the `assumptions.yaml` file.

## Output

- The script will generate CSV files containing the solution for each scenario. These files will be saved in the `results/` directory. Results will be saved as `scenario_X.csv` where `X` is the scenario number.

## License

This project is licensed under the MIT License. Please refer to the `LICENSE` file for more details.

## Contact

For questions or contributions, please contact: alexander.hoogsteyn@kuleuven.be
