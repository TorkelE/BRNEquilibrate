### File for testing the plotting functions.

#Tests plottings for normal one dimensional bifurcation diagrams.
for rn in reaction_networks_standard
    p = 10*rand(length(rn.params))
    p_bif = rand(rn.params)
    p_range = (0.4:0.4:2.0)*p[1]
    bif = bifurcation_diagram(rn,p,p_bif,p_range)
    plot(bif)
    plot!(bif)
    plot(bif,rand(1:length(rn.syms)))
    plot!(bif,rand(1:length(rn.syms)))
    plot(bif,rand(rn.syms))
    plot!(bif,rand(rn.syms))
    plot_bifurcation_diagram(rand(1:length(rn.syms)),rn,p,p_bif,p_range)
    plot_bifurcation_diagram!(rand(1:length(rn.syms)),rn,p,p_bif,p_range)
    plot_bifurcation_diagram(rand(rn.syms),rn,p,p_bif,p_range)
    plot_bifurcation_diagram!(rand(rn.syms),rn,p,p_bif,p_range)
end

#Tests plottings fortwo dimensional bifurcation diagrams.
for rn in reaction_networks_standard
    p = 10*rand(length(rn.params))
    p_bif_1 = rand(rn.params)
    p_range_1 = (1.0:2.0)*p[1]
    p_bif_2 = rand(setdiff(rn.params,Set([p_bif_1])))
    p_range_2 = (1.0:2.0)*p[1]
    bif_2d = bifurcation_diagram_2d(rn,p,p_bif_1,p_bif_2,p_range_1,p_range_2)
    plot(bif_2d)
    plot!(bif_2d)
    plot(bif_2d,rand(1:length(rn.syms)))
    plot!(bif_2d,rand(1:length(rn.syms)))
    plot(bif_2d,rand(rn.syms))
    plot!(bif_2d,rand(rn.syms))
    plot_bifurcation_diagram_2d(rand(1:length(rn.syms)),rn,p,p_bif_1,p_bif_2,p_range_1,p_range_2)
    plot_bifurcation_diagram_2d!(rand(1:length(rn.syms)),rn,p,p_bif_1,p_bif_2,p_range_1,p_range_2)
    plot_bifurcation_diagram_2d(rand(rn.syms),rn,p,p_bif_1,p_bif_2,p_range_1,p_range_2)
    plot_bifurcation_diagram_2d!(rand(rn.syms),rn,p,p_bif_1,p_bif_2,p_range_1,p_range_2)
end
