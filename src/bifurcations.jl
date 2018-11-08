import RecipesBase.plot, RecipesBase.plot!

function bifurcation_color(stability)
    stability[1]&&stability[2]&&(return :yellow)
    stability[1]&&(return :red)
    stability[2]&&(return :green)
    return :blue
end
