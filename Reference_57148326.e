ERROR: LoadError: ArgumentError: Package Distributions not found in current path.
- Run `import Pkg; Pkg.add("Distributions")` to install the Distributions package.
Stacktrace:
 [1] macro expansion
   @ ./loading.jl:1163 [inlined]
 [2] macro expansion
   @ ./lock.jl:223 [inlined]
 [3] require(into::Module, mod::Symbol)
   @ Base ./loading.jl:1144
in expression starting at /vsc-hard-mounts/leuven-data/351/vsc35199/GreenIndustrialProduction/demos/rollingHorizonStochastic.jl:7
srun: error: r25i13n19: task 0: Exited with exit code 1
