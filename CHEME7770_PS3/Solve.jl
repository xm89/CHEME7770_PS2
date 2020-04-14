using GLPK

include("Flux.jl")
include("DataDictionary.jl")

# Outputs:
# `objective_value` - value of the objective function at the optimum
# `calculated_flux_array` - R x 1 flux array at the optimum
# `dual_value_array` - R x 1 dual values
# `uptake_array` - M x 1 array of S*v
# `exit_flag` = 0 if optimal
# `status_flag` = 5 if optimal

data_dictionary = build_data_dictionary()

# solve the lp problem -
(objective_value, flux_array, dual_array, uptake_array, exit_flag) = calculate_optimal_flux_distribution(data_dictionary)
