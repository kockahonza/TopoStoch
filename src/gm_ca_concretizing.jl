################################################################################
# "Concretizing" - plugging values for symbolic terms and making everything Float
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
function get_variables(em::EnergyMatrices)
    union(get_variables(em.monomer), get_variables(em.interactions))
end
function get_variables(ca::ComplexAllosteryGM{S,Num}) where {S}
    variables = Set{Num}()
    for edge in edges(ca.graph)
        for var in Symbolics.get_variables(weight(edge))
            push!(variables, var)
        end
    end
    union(variables, get_variables(ca.energy_matrices))
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
function ssubstitute(em::EnergyMatrices, terms)
    EnergyMatrices(ssubstitute(em.monomer, terms), ssubstitute(em.interactions, terms))
end
function ssubstitute(ca::ComplexAllosteryGM{S,F}, terms::Dict) where {S,F}
    new_energy_matrices = ssubstitute(ca.energy_matrices, terms)
    new_graph = SimpleWeightedDiGraph(
        ssubstitute(adjacency_matrix(ca.graph), terms)
    )
    new_metadata = Dict{Any,Any}(terms)
    new_metadata["old_metadata"] = ca.metadata
    ComplexAllosteryGM(ca.N, ca.C, ca.B;
        symmetry=S(),
        numtype=Val(F),
        energy_matrices=new_energy_matrices,
        graph=new_graph,
        metadata=new_metadata
    )
end

"""
Substitute all variables in a symbolic Array, EnergyMatrices or
ComplexAllosteryGM making it a concrete Float64 version with numbers only.
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
function substitute_to_float(em::EnergyMatrices, terms::Dict{Num,Float64})
    EnergyMatrices{Float64}(substitute_to_float(em.monomer, terms), substitute_to_float(em.interactions, terms))
end
function substitute_to_float(ca::ComplexAllosteryGM{S}, terms::Dict{Num,Float64}) where {S}
    new_energy_matrices = substitute_to_float(ca.energy_matrices, terms)
    new_graph = SimpleWeightedDiGraph(
        substitute_to_float(adjacency_matrix(ca.graph), terms)
    )
    new_metadata = Dict{Any,Any}(terms)
    new_metadata["old_metadata"] = ca.metadata
    ComplexAllosteryGM(ca.N, ca.C, ca.B;
        symmetry=S(),
        numtype=Val(Float64),
        energy_matrices=new_energy_matrices,
        graph=new_graph,
        environment=(terms[Symbolics.variable(:μ)], terms[Symbolics.variable(:kT)]),
        metadata=new_metadata
    )
end

function get_test_terms(args...)
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
function make_factory(em::EnergyMatrices, variables=Num.(get_variables(em)); kwargs...)
    mono_fact = make_factory(em.monomer, variables; kwargs...)
    int_fact = make_factory(em.interactions, variables; kwargs...)

    function (args...)
        EnergyMatrices(
            Matrix{Float64}(mono_fact(args...)),
            Matrix{Float64}(int_fact(args...))
        )
    end
end
function make_factory(ca::ComplexAllosteryGM{S}, variables=Num.(get_variables(ca)); kwargs...) where {S}
    em_fact = make_factory(ca.energy_matrices, variables; kwargs...)
    adjmat_fact = make_factory(adjacency_matrix(ca.graph), variables; kwargs...)

    function (args...)
        ComplexAllosteryGM(ca.N, ca.C, ca.B;
            symmetry=S(),
            numtype=Val(Float64),
            energy_matrices=em_fact(args...),
            graph=SimpleWeightedDiGraph(adjmat_fact(args...))
        )
    end
end
