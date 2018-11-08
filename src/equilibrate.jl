
function steady_states(model::DiffEqBase.AbstractReactionNetwork)

end

function stability(fp::Vector{Flot64},model::DiffEqBase.AbstractReactionNetwork)
    jacobian = map(val->recursive_replace(val,Dict(zip(model.syms),fp)),deepcopy(model.symjac))
    eigenvalues = eigvals(jacobian)
    return (maximum(real(eigenvalues))>0,any(imag(eigenvalues))>0)
end
