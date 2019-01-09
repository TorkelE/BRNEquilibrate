### Main Module for Biochemical Reaction Networks Equilibrate. Exports appropraite functions ###
module BRNEquilibrate

#Fetches required packages
using DiffEqBase
using HomotopyContinuation
using LinearAlgebra
using DynamicPolynomials
using Plots
using ProgressMeter

#Runs the various file contained in the package.
include("equilibrate.jl")
include("bifurcations.jl")
include("plotting.jl")

#Exports the functions that should be avaiable to the user.
export steady_states, solve_poly, bifurcation_diagram, bifurcation_diagram_2d, stability, steady_states_stabiity
export plot, plot!, plot_bifurcation_diagram, plot_bifurcation_diagram!, plot_bifurcation_diagram_2d, plot_bifurcation_diagram_2d!

end
