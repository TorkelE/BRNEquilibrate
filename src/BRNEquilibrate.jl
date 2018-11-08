
module BRNEquilibrate

using DiffEqBase
using LinearAlgebra
using DynamicPolynomials

include("equilibrate.jl")
include("bifurcations.jl")

BRN_test_func() = println("BRNEquilibrate is avaiable")

export BRN_test_func

end # module
