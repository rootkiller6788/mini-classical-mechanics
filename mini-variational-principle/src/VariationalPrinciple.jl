module VariationalPrinciple
using LinearAlgebra
include("functional.jl")
include("euler_lagrange.jl")
include("constrained.jl")
include("second_variation.jl")
include("optimal_control.jl")
include("field_theory.jl")
include("hamilton_principle.jl")
include("numerical_variational.jl")
include("multiscale.jl")
end
