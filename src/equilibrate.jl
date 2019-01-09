#Finds the steady states of a biochemical reaction network.
function steady_states(model::DiffEqBase.AbstractReactionNetwork,p::Vector,replacements::Expr...)
    expr = meta_program(deepcopy(model.f_func), model.syms, Dict(zip(model.params, p)),replacements...)
    return eval(expr)
end

#Finds the steady states of a biochemical reaction network and the corresponding stability, returns the result in an array.
function steady_states_stabiity(model::DiffEqBase.AbstractReactionNetwork,p::Vector,replacements::Expr...)
    expr = meta_program(deepcopy(model.f_func), model.syms, Dict(zip(model.params, p)),replacements...)
    return [(fp,stability(fp,model,p)) for fp in eval(expr)]
end

#Takes the reaction rates in a biochemical reaction networks and pieces together an expression which can be evaluated to generate the steady states.
function meta_program(eqs::Vector{Expr},vars::Vector{Symbol},pars::Dict,replacements::Expr...;allowed_exp_params::Set{Symbol}=Set{Symbol}())
    poly_vars = :(@polyvar)
    foreach(var -> push!(poly_vars.args,var), vars)
    polynomial = :(internal___var___poly = [])
    foreach(eq -> push!(polynomial.args[2].args,eq), eqs)
    my_replace!(polynomial,vars,replacements...)
    recursive_check_powers!(recursive_replace!(polynomial,pars),allowed_exp_params)
    return Expr(:block, poly_vars, polynomial, :(solve_poly(internal___var___poly)))
end

#Takes the reaction rates in a biochemical reaction networks and pieces together an expression which can be evaluated to generate the steady states.
function my_replace!(polynomial::Expr,syms::Vector{Symbol},replacements::Expr...)
    replaced = Set{Symbol}()
    for rep in replacements
        next_replace = recursive_find_next_replacement(rep,setdiff(syms,replaced))
        (next_replace!=nothing) ? push!(replaced,next_replace) : continue
        polynomial.args[2].args[findfirst(syms.==next_replace)] = balance_poly(deepcopy(rep))
    end
end

#Take a polynomial enterd with a '=' and makes it to a single side.
function balance_poly(poly::Expr)
    return :($(poly.args[1])-$(poly.args[2].args[2]))
end

#Solves a polynomial system using HomotopyContinuation, but discards solutions not of interest for BRNs.
function solve_poly(poly)
    in(:num,fieldnames(typeof(poly[1]))) ? (polynomial = [poly[1].num]) : (polynomial = [poly[1]])
    for i = 2:length(poly)
        in(:num,fieldnames(typeof(poly[i]))) ? (push!(polynomial, poly[i].num)) : (push!(polynomial, poly[i]))
    end
    results = HomotopyContinuation.solve(polynomial,report_progress=false)
    results = realsolutions(results; tol=1e-6, onlynonsingular=false, singulartol=1e10, onlyfinite=true)
    final_results = Vector{Vector{Float64}}()
    foreach(res -> (minimum(res)>-0.001) && push!(final_results,res), results)
    return final_results
end

#Parses an expression and inserts number instead of symbols.
function recursive_replace!(expr::Any,pars::Dict)
    haskey(pars,expr) && return pars[expr]
    (typeof(expr) != Expr) && return expr
    for i = 1:length(expr.args)
        expr.args[i] = recursive_replace!(expr.args[i], pars)
    end
    return expr
end

#Find the reactant in a fixed concentration equation for which to replace its reaction rate equation.
function recursive_find_next_replacement(expr::Any,syms::Array{Symbol,1})
    (typeof(expr)==Symbol)&&in(expr,syms)&&(return expr)
    (typeof(expr) != Expr)&&(return nothing)
    for arg in expr.args
        output = recursive_find_next_replacement(arg,syms)
        (output!=nothing) && (return output)
    end #Returns nothing if for loop passed.
end

#Checks all of the powers in an expressions, ensuring they have the correct form (i.e. not a float).
function recursive_check_powers!(expr::Any,allowed_params::Set{Symbol}=Set{Symbol}())
    !(typeof(expr) <: Expr) && return
    foreach(arg -> recursive_check_powers!(arg,allowed_params), expr.args)
    if (expr.args[1] == :^)
        exp = expr.args[3]
        if (typeof(exp)==Expr)&&(exp.head==:call)&&(exp.args[1]==:-)&&(length(exp.args)==2)&&(typeof(exp.args[2])<:Number)
            expr.args[3] = -exp.args[2]
            exp = expr.args[3]
        end
        in(exp,allowed_params) && return
        !(typeof(exp)<:Number) && error("Error: Bad exponent: ", exp)
        (exp%1.0 != 0) && error("Error: Non integer exponent.")
        expr.args[3] = Int64(exp)
    end
end

#Calculates the stability of a given fixed point for a given reaction network and parameter values.
function stability(fp::Vector{Float64},model::DiffEqBase.AbstractReactionNetwork,pars::Vector)
    jacobian = map(val->recursive_replace!(recursive_replace!(val,Dict(zip(model.syms,fp))),Dict(zip(model.params,pars))),deepcopy(model.symjac))
    eigenvalues = eigvals(eval.(jacobian))
    return (maximum(real(eigenvalues))<1e-6,any(imag(eigenvalues).>1e-6))
end
