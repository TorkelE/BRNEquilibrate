import RecipesBase.plot, RecipesBase.plot!

struct bifurcation_diagram
    param::Symbol
    param_values::Vector{Float64}
    fixed_points::Array{Float64,2}
    stabilities::Array{Tuple{Bool,Bool},2}
    colors::Array{color,2}
end

struct bifurcation_diagram_2d
    param1::Symbol
    param2::Symbol
    param_values::Matrix{Float64}
    fixed_points::Matrix{Array{Float64,2}}
    stabilities::Matrix{Array{Tuple{Bool,Bool},2}}
    colors::Matrix{Array{color,2}}
end

function bifurcation_color(stability)
    stability[1]&&stability[2]&&(return :yellow)
    stability[1]&&(return :red)
    stability[2]&&(return :green)
    return :blue
end

function plot(bifurcation_diagram)
end
function plot!(bifurcation_diagram)
end
function plot(bifurcation_diagram_2d)
end
function plot!(bifurcation_diagram_2d)
end
