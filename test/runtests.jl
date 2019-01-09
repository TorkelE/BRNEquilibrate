### Tests all functionality provided by the oackage.

#Fetches the packages required to run the tests.
using BRNEquilibrate
using DifferentialEquations
using Random
using Test

#Runs the tests.
@time begin
  @time @testset "Model Macro" begin include("reaction_networks.jl") end
  @time @testset "Model Macro" begin include("test_equilibrate.jl") end
  @time @testset "Model Macro" begin include("test_bifurcations.jl") end
  @time @testset "Model Macro" begin include("test_plots.jl") end
end
