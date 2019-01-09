### File for testing that the equilibrate function finds the proper steady states and stabilities.

#Tests on the standard networks.
for rn in reaction_networks_standard                                                #For every reaction network in the set.
    p_vals = map(amp->amp*rand(length(rn.params)),[1.,10.,100.])                    #Generates three random sets of parameter values.
    for p in p_vals                                                                 #For each parameter value set.
        fps = steady_states(rn,p,stability_type=true)                               #Finds the fixed points..
        stab = map(fp->stability(fp[1],rn,p), fps)                                  #Calculates stability using the stability calculation function
        @test [entry[2] for entry in fps]==stab                                     #Tests that the stability values corresponds.
        for fp in fps                                                               #For each fixed point found.
            fp[2][1] ? (tend = 10) : (tend = -10)                                   #Checks if fixed point is unstable, in that case we will run solver in negative time.
            end_point = solve(ODEProblem(rn,fp[1],(0.,tend),p),Rosenbrock23())[end] #Starts a simulation in the fixed point found.
            @test maximum(end_point-fp[1])<0.0001                                   #Checks that we are still at the same point.
        end
    end
end

#Tests on the networks containing hill functions.
for rn in reaction_networks_hill
    p_vals = map(amp->rand(1:amp,length(rn.params)),[3,5,10])
    for p in p_vals
        fps = steady_states(rn,p,stability_type=true)
        stab = map(fp->stability(fp[1],rn,p), fps)
        @test [entry[2] for entry in fps]==stab
        for fp in fps
            fp[2][1] ? (tend = 10) : (tend = -10)
            end_point = solve(ODEProblem(rn,fp[1],(0.,tend),p),Rosenbrock23())[end]
            @test maximum(end_point-fp[1])<0.0001
        end
    end
end

#Tests on the networks requiring fixed concentrations.
for i = 1:length(reaction_networks_fixed_conc)
    p_vals = map(amp->amp*rand(length(reaction_networks_fixed_conc[i].params)),[1.,10.,100.])
    for p in p_vals
        fps = steady_states(reaction_networks_fixed_conc[i],p,stability_type=true,fixed_conc[i]...)
        stab = map(fp->stability(fp[1],reaction_networks_fixed_conc[i],p), fps)
        @test [entry[2] for entry in fps]==stab
        for fp in fps
            fp[2][1] ? (tend = 10) : (tend = -10)
            end_point = solve(ODEProblem(reaction_networks_fixed_conc[i],fp[1],(0.,tend),p),Rosenbrock23())[end]
            @test maximum(end_point-fp[1])<0.0001
        end
    end
end
