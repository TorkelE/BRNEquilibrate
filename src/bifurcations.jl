#Creates a bifurcation diagram of a specific reaction network and for the given parameter and range.
function bifurcation_diagram(model::DiffEqBase.AbstractReactionNetwork, pOrg::AbstractVector, parameter::Symbol, range::AbstractVector, replacements::Expr...;mult=false::Bool,progress_update=2147483640.::Float64)
    p_idx = findfirst(model.params.==parameter)
    mult && (range*=pOrg[p_idx])
    params_const = Dict(zip(deleteat!(deepcopy(model.params),p_idx), deleteat!(deepcopy(pOrg),p_idx)))
    expr = meta_program(deepcopy(model.f_func), model.syms, params_const,replacements...,allowed_exp_params=Set{Symbol}([parameter]))
    jacobian = map(val->recursive_replace!(val,params_const),deepcopy(model.symjac))
    fps = Vector{Vector{Vector{Float64}}}(undef,length(range))
    stabilities = Vector{Vector{Tuple{Bool,Bool}}}(undef,length(range))
    colors = Vector{Vector{Symbol}}(undef,length(range))
    @showprogress progress_update for i = 1:length(range)
        fps[i] = eval(recursive_replace!(deepcopy(expr),Dict(parameter=>range[i])))
        stabilities[i] = map(fp->stability_internal(jacobian,Dict(zip(model.syms,fp)),Dict(parameter=>range[i])),fps[i])
        colors[i] = bifurcation_color.(stabilities[i])
    end
    return bif_diagram(parameter,range,fps,stabilities,colors,deepcopy(model.syms))
end

#Bifurcation diagram structure, contains all the information in a bifurcation diagram.
struct bif_diagram
    param::Symbol
    param_values::AbstractVector
    fixed_points::Vector{Vector{Vector{Float64}}}
    stabilities::Vector{Vector{Tuple{Bool,Bool}}}
    colors::Vector{Vector{Symbol}}
    reactant_syms::Vector{Symbol}
end

#Creates a two dimensional bifurcation diagram of a specific reaction network and for the given parameters and ranges.
function bifurcation_diagram_2d(model::DiffEqBase.AbstractReactionNetwork, pOrg::AbstractVector, parameter1::Symbol, parameter2::Symbol, range1::AbstractVector, range2::AbstractVector, replacements::Expr...;mult=false::Bool,progress_update=2147483640.::Float64)
    p_idx1 = findfirst(model.params.==parameter1); p_idx2 = findfirst(model.params.==parameter2);
    mult && (range1*=pOrg[p_idx1]; range2*=pOrg[p_idx2];)
    params_const = Dict(zip(deleteat!(deepcopy(model.params),sort([p_idx1,p_idx2])), deleteat!(deepcopy(pOrg),sort([p_idx1,p_idx2]))))
    expr = meta_program(deepcopy(model.f_func), model.syms, params_const,replacements...,allowed_exp_params=Set{Symbol}([parameter1,parameter2]))
    jacobian = map(val->recursive_replace!(val,params_const),deepcopy(model.symjac))
    fps = Matrix{Vector{Vector{Float64}}}(undef,length(range1),length(range2))
    stabilities = Matrix{Vector{Tuple{Bool,Bool}}}(undef,length(range1),length(range2))
    colors = Matrix{Vector{Symbol}}(undef,length(range1),length(range2))
    @showprogress progress_update for i = 1:length(range1), j = 1:length(range2)
        fps[i,j] = eval(recursive_replace!(deepcopy(expr),Dict(parameter1=>range1[i],parameter2=>range2[j])))
        stabilities[i,j] = map(fp->stability_internal(jacobian,Dict(zip(model.syms,fp)),Dict(parameter1=>range1[i],parameter2=>range2[j])),fps[i,j])
        colors[i,j] = bifurcation_color.(stabilities[i,j])
    end
    return bif_diagram_2d(parameter1,parameter2,range1,range2,fps,stabilities,colors,deepcopy(model.syms))
end

#Two dimensional bifurcation diagram structure, contains all the information in a two dimensional bifurcation diagram.
struct bif_diagram_2d
    param1::Symbol
    param2::Symbol
    param_values1::AbstractVector
    param_values2::AbstractVector
    fixed_points::Matrix{Vector{Vector{Float64}}}
    stabilities::Matrix{Vector{Tuple{Bool,Bool}}}
    colors::Matrix{Vector{Symbol}}
    reactant_syms::Array{Symbol}
end

#Function for determening stability of a reaction network. Only used internally by bifurcation diagram functions. Not to be used by the user.
function stability_internal(jacobian::Matrix{Expr},var_vals::Dict,par_vals::Dict)
    jac = map(val->recursive_replace!(recursive_replace!(val,var_vals),par_vals),deepcopy(jacobian))
    eigenvalues = eigvals(eval.(jac))
    return (maximum(real(eigenvalues))<1e-6,any(imag(eigenvalues).>1e-6))
end

#Determines the color of a point in the bifurcation diagram for the given stability information.
function bifurcation_color(stability::Tuple{Bool,Bool})
    stability[1]&&stability[2]&&(return :green)
    stability[1]&&(return :blue)
    stability[2]&&(return :yellow)
    return :red
end
