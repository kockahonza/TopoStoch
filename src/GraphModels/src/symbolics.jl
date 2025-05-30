################################################################################
# Substitutions and concretizing symbolic objects
################################################################################
get_variables(::Nothing) = Num[]
get_variables(term::Num) = unique(Num.(Symbolics.get_variables(term)))
get_variables(term::Symbolics.SymbolicUtils.Symbolic) = unique(Num.(Symbolics.get_variables(term)))
function get_variables(arr::Union{AbstractArray,Tuple})
    variables = Num[]
    for term in arr
        for var in get_variables(term)
            if all(x -> !isequal(x, var), variables)
                push!(variables, Num(var))
            end
        end
    end
    variables
end
function get_variables(gm::AbstractGraphModel{<:Num})
    throw(ErrorException(f"No method of \"get_variables\" was provided for type \"{typeof(gm)}\""))
end
export get_variables

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
export ssubstitute

"""
Substitute all variables in a symbolic object and make it concrete with real
floating point numbers only only (using Float64).
"""
function substitute_to_float(anything, terms::Dict{Num,Float64})
    if eltype(anything) <: Number
        map(t -> Float64(Symbolics.unwrap(substitute(t, terms))), anything)
    else
        map(t -> substitute_to_float(t, terms), anything)
    end
end
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
function substitute_to_float(gm::AbstractGraphModel{<:Num}, vals::AbstractVector)
    substitute_to_float(gm, Dict(get_variables(gm) .=> Float64.(vals)))
end
export substitute_to_float

function terms_full_test(args...)
    terms = Dict{Num,Float64}()
    for var in get_variables(args...)
        terms[var] = 1.0
    end
    terms
end
export terms_full_test

"""
Achieves a similar goal to the above but much faster especially for multiple calls.
Each of these createsm a "factory" for making the given object given values for all
the variables in it. These are then very fast to use. Returned closure internally
converts all arguments to Float64 as otherwise I get runtime errors.
"""
function make_factory(obj, variables=get_variables(obj); kwargs...)
    num_vars = length(variables)

    bfunc_ = build_function(obj, variables; expression=Val(false), kwargs...)
    bfunc = isa(obj, AbstractArray) ? bfunc_[1] : bfunc_

    try
        bfunc([1.0 for _ in 1:num_vars])
    catch UndefVarError
        throw(ArgumentError(f"not all variables were provided, concrete factory cannot be made"))
    end

    let num_vars = num_vars, bfunc = bfunc
        function (args...)
            if length(args) != num_vars
                throw(ArgumentError(f"closure expected {num_vars} arguments but received {length(args)}"))
            end
            bfunc(SVector(convert.(Float64, args)))
        end
    end
end
function make_factory(gm::AbstractGraphModel{<:Num}, variables=get_variables(arr))
    throw(ErrorException(f"No method of \"make_factory\" was provided for type \"{typeof(gm)}\""))
end
export make_factory

################################################################################
# Util
################################################################################
"""
Make a general symbolic version of a matrix - one variable for each nonzero component
"""
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
export safemat_general

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
export substitute_safenames

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
export get_rewriter

function convert(::Type{T}, num::Num) where {T<:Number}
    up = Symbolics.unwrap(num)
    convert(T, up)
end

# Helpers for defining variables
function tuplevariables(args...)
    tuple(variables(args...)...)
end
# NOTE: These are definitely not type stable!
function savariables(symbol, args...)
    lenghts = filter(x -> x != 1, length.(args))
    SArray{Tuple{lenghts...},Num,length(lenghts)}(variables(symbol, args...)...)
end
function listsavariables(symbols...)
    SVector(variable.(symbols)...)
end
export tuplevariables, savariables, listsavariables
