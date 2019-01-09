#Imports plot recipes (so that they can be expanded).
import RecipesBase.plot, RecipesBase.plot!

#Plots for single dimensional bifurcation diagrams.
function plot(bif_diagram::bif_diagram,var::Int64=1)
    plot(); return plot!(bif_diagram,var);
end
function plot!(bif_diagram::bif_diagram,var::Int64=1)
    for i = 1:length(bif_diagram.param_values), j = 1:length(bif_diagram.fixed_points[i])
        scatter!((bif_diagram.param_values[i],bif_diagram.fixed_points[i][j][var]),color=bif_diagram.colors[i][j],label="")
    end
    return plot!()
end

#Plots for two dimensional bifurcation diagrams.
function plot(bif_diagram_2d::bif_diagram_2d,var::Int64=1)
    plot(); return plot!(bif_diagram_2d,var);
end
function plot!(bif_diagram_2d::bif_diagram_2d,var::Int64=1)
    for i = 1:length(bif_diagram_2d.param_values1), j = 1:length(bif_diagram_2d.param_values2), k = 1:length(bif_diagram_2d.fixed_points[i,j])
        scatter!((bif_diagram_2d.param_values1[i],bif_diagram_2d.param_values2[j],bif_diagram_2d.fixed_points[i,j][k][var]),color=bif_diagram_2d.colors[i,j][k],label="")
    end
    return plot!()
end

#These bifurcation diagrams takes a variable input as a symbol, and plots that variable in the bifurcation diagram.
function plot(bif_diagram::bif_diagram,var::Symbol)
    !in(var,bif_diagram.reactant_syms)&&error("The reactant ",var," do not seem to be a part of the original reaction network.")
    plot(bif_diagram,findfirst(bif_diagram.reactant_syms.==var))
end
function plot!(bif_diagram::bif_diagram,var::Symbol)
    !in(var,bif_diagram.reactant_syms)&&error("The reactant ",var," do not seem to be a part of the original reaction network.")
    plot!(bif_diagram,findfirst(bif_diagram.reactant_syms.==var))
end
function plot(bif_diagram_2d::bif_diagram_2d,var::Symbol)
    !in(var,bif_diagram_2d.reactant_syms)&&error("The reactant ",var," do not seem to be a part of the original reaction network.")
    plot(bif_diagram_2d,findfirst(bif_diagram_2d.reactant_syms.==var))
end
function plot!(bif_diagram_2d::bif_diagram_2d,var::Symbol)
    !in(var,bif_diagram_2d.reactant_syms)&&error("The reactant ",var," do not seem to be a part of the original reaction network.")
    plot!(bif_diagram_2d,findfirst(bif_diagram_2d.reactant_syms.==var))
end

#Creates plot_bifurcation_diagram function, which both calculates and plots a bifurcation diagram.
plot_bifurcation_diagram(var::Union{Int64,Symbol},args...) = plot(bifurcation_diagram(args...),var)
plot_bifurcation_diagram!(var::Union{Int64,Symbol},args...) = plot!(bifurcation_diagram(args...),var)
plot_bifurcation_diagram_2d(var::Union{Int64,Symbol},args...) = plot(bifurcation_diagram_2d(args...),var)
plot_bifurcation_diagram_2d!(var::Union{Int64,Symbol},args...) = plot!(bifurcation_diagram_2d(args...),var)
