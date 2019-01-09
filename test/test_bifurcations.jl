### File for testing the bifurcation functions.

#Tests on the standard networks.
for rn in reaction_networks_standard
    p_vals = map(amp->amp*rand(length(rn.params)),[1.,10.,100.])
    for p in p_vals
        p_bif = rand(rn.params)
        p_range = (0.4:0.4:2.0)*p[1]
        bif = bifurcation_diagram(rn,p,p_bif,p_range)
        for i = 1:5
            p[findfirst(p_bif.==rn.params)] = p_range[i]
            for j = 1:length(bif.fixed_points[i])
                bif.stabilities[i][j][1] ? (tend = 10) : (tend = -10)
                end_point = solve(ODEProblem(rn,bif.fixed_points[i][j],(0.,tend),p),Rosenbrock23())[end]
                @test maximum(end_point-bif.fixed_points[i][j])<0.0001
            end
        end
    end
end

#Tests on the networks containing hill functions.
for rn in reaction_networks_hill
    p_vals = map(amp->rand(1:amp,length(rn.params)),[3,5,10])
    for p in p_vals
        p_bif = rand(rn.params)
        p_range = 1:5
        bif = bifurcation_diagram(rn,p,p_bif,p_range)
        for i = 1:5
            p[findfirst(p_bif.==rn.params)] = p_range[i]
            for j = 1:length(bif.fixed_points[i])
                bif.stabilities[i][j][1] ? (tend = 10) : (tend = -10)
                end_point = solve(ODEProblem(rn,bif.fixed_points[i][j],(0.,tend),p),Rosenbrock23())[end]
                @test maximum(end_point-bif.fixed_points[i][j])<0.0001
            end
        end
    end
end

#Tests on the networks requiring fixed concentrations.
for i = 1:length(reaction_networks_fixed_conc)
    p_vals = map(amp->amp*rand(length(reaction_networks_fixed_conc[i].params)),[1.,10.,100.])
    for p in p_vals
        p_bif = rand(reaction_networks_fixed_conc[i].params)
        p_range = (0.4:0.4:2.0)*p[1]
        bif = bifurcation_diagram(reaction_networks_fixed_conc[i],p,p_bif,p_range,fixed_conc[i]...)
        for I = 1:5
            p[findfirst(p_bif.==reaction_networks_fixed_conc[i].params)] = p_range[I]
            for j = 1:length(bif.fixed_points[I])
                fps = steady_states(reaction_networks_fixed_conc[i],p,fixed_conc[i]...,stability_type=true)
                bif.stabilities[I][j][1] ? (tend = 10) : (tend = -10)
                end_point = solve(ODEProblem(reaction_networks_fixed_conc[i],bif.fixed_points[I][j],(0.,tend),p),Rosenbrock23())[end]
                @test maximum(end_point-bif.fixed_points[I][j])<0.0001
            end
        end
    end
end

#Test that the functionality to allow multiplication of original values workd properly (only tested on standard networks to avoid making tests run for to long).
for rn in reaction_networks_standard
    p_vals = map(amp->amp*rand(length(rn.params)),[1.,10.,100.])
    for p in p_vals
        p_bif = rn.params[1]
        p[1] = 2.0
        p_range_1 = [1.0,1.5,2.0,3.0,4.0]
        p_range_2 = [0.5,0.75,1.0,1.5,2.0]
        bif_1 = bifurcation_diagram(rn,p,p_bif,p_range_1)
        bif_2 = bifurcation_diagram(rn,p,p_bif,p_range_2,mult=true)
        @test bif_1.param_values==bif_2.param_values
        for i = 1:5
            fps1 = sort(bif_1.fixed_points[i])
            fps2 = sort(bif_1.fixed_points[i])
            @test length(fps1)==length(fps2)
            for j = 1:length(fps1)
                @test maximum(fps1[j]-fps2[j])<0.0001
            end
        end
    end
end

#Test to ensure the 2d bifurcation diagrams are found properly (only tested on standard networks to avoid making tests run for to long).
for rn in reaction_networks_standard
    p_vals = map(amp->amp*rand(length(rn.params)),[1.,10.,100.])
    for p in p_vals
        p_bif_1 = rand(rn.params)
        p_range_1 = (1.0:2.0)*p[1]
        p_bif_2 = rand(setdiff(rn.params,Set([p_bif_1])))
        p_range_2 = (1.0:2.0)*p[1]
        bif = bifurcation_diagram_2d(rn,p,p_bif_1,p_bif_2,p_range_1,p_range_2)
        for i1 = 1:2, i2 = 1:2
            p[findfirst(p_bif_1.==rn.params)] = p_range_1[i1]
            p[findfirst(p_bif_2.==rn.params)] = p_range_2[i2]
            for j = 1:length(bif.fixed_points[i1,i2])
                bif.stabilities[i1,i2][j][1] ? (tend = 10) : (tend = -10)
                end_point = solve(ODEProblem(rn,bif.fixed_points[i1,i2][j],(0.,tend),p),Rosenbrock23())[end]
                @test maximum(end_point-bif.fixed_points[i1,i2][j])<0.0001
            end
        end
    end
end
