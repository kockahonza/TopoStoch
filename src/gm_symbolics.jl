################################################################################
# Concretizing symbolic objects
################################################################################
get_variables(term) = Set{Num}(Symbolics.get_variables(term))
function get_variables(arr::AbstractArray)
    variables = Set{Num}()
    for term in arr
        for var in get_variables(term)
            push!(variables, Num(var))
        end
    end
    variables
end
get_variables(::Nothing) = Set{Num}()
function get_variables(gm::AbstractGraphModel{<:Num})
    throw(ErrorException(f"No method of \"get_variables\" was provided for type \"{typeof(gm)}\""))
end

"""
(Simple)Substitute, works almost like Symbolics substitute but
is defined on the custom types in this project and will always
use Num, not any other funny types
"""
ssubstitute(concrete_num::Number, _) = concrete_num
ssubstitute(num::Num, terms) = substitute(num, terms)
function ssubstitute(arr::AbstractArray, terms)
    new_arr = similar(arr)
    for i in eachindex(arr, new_arr)
        new_arr[i] = ssubstitute(arr[i], terms)
    end
    new_arr
end
function ssubstitute(gm::AbstractGraphModel{<:Num}, terms)
    throw(ErrorException(f"No method of \"ssubstitute\" was provided for type \"{typeof(gm)}\""))
end

"""
Substitute all variables in a symbolic object and make it concrete with real
floating point numbers only only (using Float64).
"""
function substitute_to_float_!(farr, arr, terms::Dict{Num,Float64})
    for i in eachindex(arr, farr)
        farr[i] = Float64(Symbolics.unwrap(substitute(arr[i], terms)))
    end
end
function substitute_to_float(arr::Array, terms::Dict{Num,Float64})
    farr = Array{Float64}(undef, size(arr))
    substitute_to_float_!(farr, arr, terms)
    farr
end
function substitute_to_float(arr::AbstractSparseMatrix, terms::Dict{Num,Float64})
    farr = similar(arr, Float64)
    substitute_to_float_!(farr, arr, terms)
    farr
end
function substitute_to_float(gm::AbstractGraphModel{<:Num}, terms::Dict{Num,Float64})
    throw(ErrorException(f"No method of \"substitute_to_float\" was provided for type \"{typeof(gm)}\""))
end

function substitute_test_terms(args...)
    terms = Dict{Num,Float64}()
    for var in get_variables(args...)
        terms[var] = 1.0
    end
    terms[Symbolics.variable(:μ)] = 1.0
    terms[Symbolics.variable(:kT)] = 1.0
    terms
end

"""
Achieves a similar goal to the above but much faster especially for multiple calls.
Each of these createsm a "factory" for making the given object given values for all
the variables in it. These are then very fast to use.
"""
function make_factory(arr, variables=Num.(get_variables(arr)); kwargs...)
    num_vars = length(variables)

    bfunc = build_function(arr, variables; expression=Val(false), kwargs...)[1]

    try
        bfunc([1.0 for _ in 1:num_vars])
    catch UndefVarError
        throw(ArgumentError(f"not all variables were provided, concrete factory cannot be made"))
    end

    function (args...)
        if length(args) != num_vars
            throw(ArgumentError(f"closure expected {num_vars} arguments but received {length(args)}"))
        end
        bfunc(SVector(args))
    end
end
function make_factory(gm::AbstractGraphModel{<:Num}, variables=Num.(get_variables(arr)))
    throw(ErrorException(f"No method of \"make_factory\" was provided for type \"{typeof(gm)}\""))
end

################################################################################
# Dealing with wolfram
################################################################################
includet(srcdir("letter_codes.jl"))

function substitute_safenames(obj)
    terms = Dict{Num,Num}()
    code = Code()
    for var in get_variables(obj)
        terms[var] = Symbolics.variable(repr(code))
        increment!(code)
    end
    rev_terms = map(reverse, collect(terms))
    ssubstitute(obj, terms), rev_terms
end

function towolfram(term::Num)
    SymbolicsMathLink.expr_to_mathematica(term)
end
function towolfram(term::AbstractMatrix)
    SymbolicsMathLink.expr_to_mathematica(collect(transpose(term)))
end

function prep_for_wolfram(obj)
    sobj, rt = substitute_safenames(obj)
    towolfram(sobj), rt
end

function make_safe(safe, obj)
    if safe
        sobj, rt = substitute_safenames(obj)
        sobj, (srlst -> ssubstitute(srlst, rt))
    else
        obj, identity
    end
end

# TODO: Change to full=false
function w_simplify(obj; full=true, safe=true)
    obj, desafe = make_safe(safe, obj)
    cmd = full ? "FullSimplify" : "Simplify"
    rslt = wcall.(cmd, obj)
    desafe(rslt)
end

function w_eigen(matrix; safe=true)
    matrix, desafe = make_safe(safe, matrix)
    rslt = wcall("Eigensystem", collect(transpose(matrix)))
    (evals=desafe(rslt[1]), evecs=desafe.(rslt[2]))
end

function w_test_eigen(matrix; eigen_func=w_eigen, safe=false)
    (evals, evecs) = eigen_func(matrix; safe)
    residues = []
    for (eval, evec) in zip(evals, evecs)
        coeffs = (matrix * evec) ./ evec
        push!(residues, w_simplify(coeffs .- eval; safe))
    end
    residues
end

################################################################################
# Util
################################################################################
"""Make a general symbolic version of a matrix - one variable for each nonzero component"""
function safemat_general(matrix)
    simple_matrix = similar(matrix, Num)
    rev_terms = Dict{Num,Num}()
    code = Code()
    for i in eachindex(matrix)
        v = matrix[i]
        if iszero(v)
            simple_matrix[i] = v
        else
            var = Symbolics.variable(repr(code))
            simple_matrix[i] = var
            increment!(code)
            rev_terms[var] = v
        end
    end
    simple_matrix, rev_terms
end

function get_rewriter()
    rule_divexp = @rule ~y / exp(~x) => ~y * exp(-~x)
    rules_explog = [
        (@rule exp(log(~z)) => ~z),
        (@rule exp(~y * log(~z)) => ~z^~y),
        (@rule exp(~x + ~y * log(~z)) => exp(~x) * ~z^~y)
    ]

    Rewriters.Prewalk(Rewriters.Chain([
        rule_divexp,
        rules_explog...,
        SymbolicUtils.default_simplifier()
    ]))
end
