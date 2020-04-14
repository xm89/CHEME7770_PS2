using CSV
using LinearAlgebra
using DataFrames
# 1 Build stoichiometric_matrix
stoichiometric_matrix = DataFrame(CSV.File("Network.csv",header=false))
stoichiometric_matrix = convert(Matrix,stoichiometric_matrix)

# 2 Build element_matrix
element_matrix = DataFrame(CSV.File("Elements.csv",header=false))
element_matrix = convert(Matrix,element_matrix)

# 3 Check element_balance
element_balance = transpose(element_matrix )*stoichiometric_matrix
