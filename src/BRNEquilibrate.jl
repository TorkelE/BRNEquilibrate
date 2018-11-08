
module BRNEquilibrate

using DiffEqBase
using HomotopyContinuation
using Linear Algebra
using DynamicPolynomials

include("equilibrate.jl")

BRN_test_func() = println("BRNEquilibrate is avaiable")

export BRN_test_func

end # module
