using DrWatson
@quickactivate "TopoStochSim"

include(srcdir("graphmodels.jl"))

using Graphs, SimpleWeightedGraphs
using GLMakie, GraphMakie, NetworkLayout
using StatsBase
using Symbolics
using StaticArrays, SparseArrays, DataFrames
using PrettyTables, Latexify
using PyFormattedStrings
using Base.Threads
using MathLink, SymbolicsMathLink

import SymbolicUtils
import SymbolicUtils.Rewriters

import Base: copy, broadcastable, show

################################################################################
# Needed combinatoric utilities
################################################################################
"""
Returns true if `perm` has length `N` and all its elements are between `base`(inclusive)
and `base` + `E`(exclusive).
"""
function validperm(perm::AbstractVector{<:Integer}, E, N; base=1)
    if length(perm) != N
        return false
    end
    for e in perm
        if (e < base) || (e >= base + E)
            return false
        end
    end
    return true
end

"""
`perm` is assumed to be a valid permutation and `E` is as in `validperm`.
Returns in the index of the given permutation.
"""
function permtoi(perm::AbstractVector{<:Integer}, E; base=1)
    1 + sum((ei - base) * E^(i - 1) for (i, ei) in enumerate(perm))
end

"""
Returns the `i`th permutation of length `N` with each element being an
integer between 1 and `E`. The ordering is the same order as in `permtoi`.
"""
function itoperm(i, E, N; base=1)
    perm = Vector{Int}(undef, N)
    for j in N:-1:1
        factor = E^(j - 1)
        perm[j] = div(i - 1, factor)
        i -= perm[j] * factor
        perm[j] += base
    end
    perm
end

################################################################################
# Complex Allostery types
################################################################################
struct CAState{T}
    conformations::Vector{T}
    occupations::Vector{T}
end
function copy((; conformations, occupations)::CAState)
    CAState(copy(conformations), copy(occupations))
end
show(io::IO, st::CAState) = print(io, f"({st.conformations}, {st.occupations})")

abstract type Symmetry end
struct Chain <: Symmetry end
struct Loop <: Symmetry end
broadcastable(s::Symmetry) = Ref(s)

struct EnergyMatrices{F<:Number}
    monomer::Matrix{F}
    interactions::Matrix{F}
end
broadcastable(em::EnergyMatrices) = Ref(em)
function copy((; monomer, interactions)::EnergyMatrices)
    EnergyMatrices(copy(monomer), copy(interactions))
end
function show(io::IO, mime::MIME"text/plain", em::EnergyMatrices)
    print(io, "monomer energies ")
    show(io, mime, em.monomer)
    print(io, "\ninteraction energies ")
    show(io, mime, em.interactions)
end

mutable struct ComplexAllosteryGM{S<:Symmetry,F<:Number,T<:Integer} <: AbstractGraphModel
    N::T # Number of monomers
    C::T # Number of possible conformational states per monomer
    B::T # Number of ligand binding sites per monomer

    energy_matrices::Union{Nothing,EnergyMatrices{F}}

    graph::SimpleWeightedDiGraph{T,F}

    version::Union{Nothing,Float64}
    environment::Union{Nothing,Tuple{F,F}} # Pair of μ and kT
    metadata::Union{Nothing,Dict}
    function ComplexAllosteryGM(N::T, C::T, B::T;
        symmetry::S=Chain(),
        numtype::Val{F}=Val(Num),
        energy_matrices=nothing,
        graph=nothing,
        version=nothing,
        environment=nothing,
        metadata=nothing
    ) where {S,F,T<:Integer}
        if !isnothing(energy_matrices) && !validEM(energy_matrices, C, B)
            throw(ArgumentError(f"the provided energy matrices do not have the right shape for the given system"))
        end

        if isnothing(graph)
            graph = SimpleWeightedDiGraph{T,F}(numstates(N, C, B))
        elseif nv(graph) != numstates(N, C, B)
            throw(ArgumentError(f"the provided graph does not have the right number of vertices for the given system"))
        end

        new{S,F,T}(N, C, B, energy_matrices, graph, version, environment, metadata)
    end
end
broadcastable(ca::ComplexAllosteryGM) = Ref(ca)
function copy(ca::ComplexAllosteryGM{S,F}) where {S,F}
    ComplexAllosteryGM(ca.N, ca.C, ca.B;
        symmetry=S(),
        numtype=Val(F),
        energy_matrices=isnothing(ca.energy_matrices) ? nothing : copy(ca.energy_matrices),
        graph=copy(ca.graph),
        environment=ca.environment,
        metadata=copy(ca.metadata)
    )
end
function show(io::IO, mime::MIME"text/plain", ca::ComplexAllosteryGM{S,F,T}) where {F,S,T}
    println(f"CAGM{{{S}, {F}, {T}}}")
    println(f"N={ca.N}, C={ca.C}, B={ca.B}")
    println(f"energy_matrices:")
    show(io, mime, ca.energy_matrices)
    print(f"\ngraph: ")
    show(io, mime, ca.graph)
    print(f"\nenvironment: ")
    show(io, mime, ca.environment)
    print(f"\nversion: ")
    show(io, mime, ca.version)
    if !isnothing(ca.metadata)
        println(f"\nmetadata:")
        show(io, mime, ca.metadata)
    end
end
graph(ca::ComplexAllosteryGM) = ca.graph

# Small utility functions
function validEM(em::EnergyMatrices, C, B)
    (size(em.monomer) == (C, B + 1)) && (size(em.interactions) == (C, C))
end
validEM(em::EnergyMatrices, ca::ComplexAllosteryGM) = validEM(em, ca.C, ca.B)

function numstates(N, C, B)
    (C * (B + 1))^N
end
numstates(ca::ComplexAllosteryGM) = numstates(ca.N, ca.C, ca.B)

function validstate(st::CAState, N, C, B)
    validperm(st.conformations, C, N) && validperm(st.occupations, B + 1, N; base=0)
end
validstate(st::CAState, ca::ComplexAllosteryGM) = validstate(st, ca.N, ca.C, ca.B)

function statetoi(st::CAState, N, C, B)
    permtoi(st.conformations, C) + (permtoi(st.occupations, B + 1; base=0) - 1) * (C^N)
end
statetoi(st::CAState, ca::ComplexAllosteryGM) = statetoi(st, ca.N, ca.C, ca.B)

function itostate(i, N, C, B)
    (iocc, iconf) = divrem(i - 1, C^N)
    CAState(itoperm(iconf + 1, C, N), itoperm(iocc + 1, B + 1, N; base=0))
end
itostate(i, ca::ComplexAllosteryGM) = itostate(i, ca.N, ca.C, ca.B)

allstates(N, C, B) = (itostate(i, N, C, B) for i in 1:numstates(N, C, B))
allstates(ca::ComplexAllosteryGM) = allstates(ca.N, ca.C, ca.B)

################################################################################
# Dealing with energy and equilibrium probability calculations
################################################################################
# Setting up the energy matrices
function make_EM_sym_general(C, B)
    monomer = Matrix{Num}(undef, C, B + 1)
    monomer[1, 1] = 0
    monomer[2:C, 1] = Symbolics.variables(:εr, 1:C-1)
    monomer[:, 2:B+1] = Symbolics.variables(:εl, 1:C, 1:B)

    interactions = Matrix{Num}(undef, C, C)
    terms = Symbolics.variables(:εb, 1:div((C - 1) * C, 2))
    term_i = 1
    for i in 1:C
        interactions[i, i] = 0
        for j in i+1:C
            interactions[i, j] = terms[term_i]
            interactions[j, i] = terms[term_i]
            term_i += 1
        end
    end

    EnergyMatrices(monomer, interactions)
end

function make_EM_sym_C2(B)
    @variables ε_t, Δε_r, ε_b

    monomer = Matrix{Num}(undef, 2, B + 1)
    monomer[1, 1] = 0
    monomer[2, 1] = Δε_r
    monomer[1, 2:B+1] .= ε_t .* collect(1:B)
    monomer[2, 2:B+1] .= monomer[1, 2:B+1] .- Δε_r

    interactions = [0 ε_b; ε_b 0]
    EnergyMatrices(monomer, interactions)
end

# Simple calculator functions
function get_neighbors(st::CAState, i, ::Chain)
    N = length(st.conformations)
    if N == 1
        return (nothing, nothing)
    end

    if i == 1
        (st.conformations[i+1], nothing)
    elseif i == N
        (st.conformations[N-1], nothing)
    else
        (st.conformations[i-1], st.conformations[i+1])
    end
end
function get_neighbors(st::CAState, i, ::Loop)
    N = length(st.conformations)
    if N == 1
        return (nothing, nothing)
    end

    if i == 1
        (st.conformations[N], st.conformations[i+1],)
    elseif i == N
        (st.conformations[N-1], st.conformations[1])
    else
        (st.conformations[i-1], st.conformations[i+1])
    end
end
function get_neighbors(st::CAState, i, _::ComplexAllosteryGM{S}) where {S}
    get_neighbors(st, i, S())
end

function calc_monomer_energy(st::CAState, i, em::EnergyMatrices)
    em.monomer[st.conformations[i], st.occupations[i]+1]
end
calc_monomer_energy(st::CAState, i, ca::ComplexAllosteryGM) = calc_monomer_energy(st, i, ca.energy_matrices)

function calc_interaction_energy(st::CAState, i, em::EnergyMatrices, s::Symmetry)
    (left, right) = get_neighbors(st, i, s)
    energy = 0
    if !isnothing(left)
        energy += em.interactions[left, st.conformations[i]]
    end
    if !isnothing(right)
        energy += em.interactions[st.conformations[i], right]
    end
    0.5 * energy
end
calc_interaction_energy(st::CAState, i, ca::ComplexAllosteryGM{S}) where {S} = calc_interaction_energy(st, i, ca.energy_matrices, S())

function calc_energy(st::CAState, em::EnergyMatrices, s::Symmetry)
    energy = 0
    for i in 1:length(st.conformations)
        energy += calc_monomer_energy(st, i, em)
        energy += calc_interaction_energy(st, i, em, s)
    end
    energy
end
function calc_energy(st::CAState, ca::ComplexAllosteryGM{S}) where {S}
    calc_energy(st, ca.energy_matrices, S())
end

calc_numligands(st::CAState) = sum(st.occupations)

calc_numofconf(st::CAState, conf_state=1) = count(x -> x == conf_state, st.conformations)

function calc_numboundaries(st::CAState, args...)
    total = 0.0
    for i in 1:length(st.conformations)
        neighbors = filter(!isnothing, get_neighbors(st, i, args...))
        total += 0.5 * count(x -> x != st.conformations[i], neighbors)
    end
    Int(total)
end

function calc_gibbs_factor(st::CAState, em::EnergyMatrices,
    s::Symmetry, μ=Symbolics.variable(:μ), kT=Symbolics.variable(:kT))
    exp(-(calc_energy(st, em, s) - μ * calc_numligands(st)) / kT)
end
function calc_gibbs_factor(st::CAState, ca::ComplexAllosteryGM{S}, args...) where {S}
    env = if length(args) > 0
        args
    elseif !isnothing(ca.environment)
        ca.environment
    else
        ()
    end
    calc_gibbs_factor(st, ca.energy_matrices, S(), env...)
end

# Calculating equilibrium observables
function calc_all_gibbs_factors(ca::ComplexAllosteryGM{S}) where {S}
    calc_gibbs_factor.(allstates(ca), ca)
end
function calc_partition_function(ca::ComplexAllosteryGM)
    sum(calc_all_gibbs_factors(ca))
end
function calc_avg_numligands(ca::ComplexAllosteryGM)
    factors = calc_all_gibbs_factors(ca)
    sum(calc_numligands.(allstates(ca)) .* factors) / sum(factors)
end
function calc_avg_energy(ca::ComplexAllosteryGM)
    factors = calc_all_gibbs_factors(ca)
    sum(calc_energy.(allstates(ca), ca) .* factors) / sum(factors)
end

################################################################################
# Dealing with transitions
################################################################################
# Adding some basic links without meaningful weights
function add_edges_base!(ca::ComplexAllosteryGM)
    # Just use a constant for all the rates, disregarding DB and all else too
    universal_rate = 0.75 / (ca.C + 1) # This means that the chance of not doing anything is ~0.25
    for vertex in 1:numstates(ca)
        state = itostate(vertex, ca)
        # Add conformational change transitions
        for i in 1:ca.N
            for new_c in 1:ca.C
                if new_c != state.conformations[i]
                    new_state = copy(state)
                    new_state.conformations[i] = new_c
                    add_edge!(ca.graph, vertex, statetoi(new_state, ca), universal_rate)
                end
            end
        end
        # Add ligand (de)binding rates
        for i in 1:ca.N
            cur_occ = state.occupations[i]
            if cur_occ > 0
                new_state = copy(state)
                new_state.occupations[i] -= 1
                add_edge!(ca.graph, vertex, statetoi(new_state, ca), universal_rate)
            end
            if cur_occ < ca.B
                new_state = copy(state)
                new_state.occupations[i] += 1
                add_edge!(ca.graph, vertex, statetoi(new_state, ca), universal_rate)
            end
        end
    end
end

function add_edges_sym!(ca::ComplexAllosteryGM)
    rates_i = 1
    function add_new_rate(v1, v2)
        if (v1 < v2)
            rate = Symbolics.variable(:r, rates_i)
            add_edge!(ca.graph, v1, v2, rate)
            add_edge!(ca.graph, v2, v1, rate)
            rates_i += 1
        end
    end
    for vertex in 1:numstates(ca)
        state = itostate(vertex, ca)
        # Add conformational change transitions
        for i in 1:ca.N
            for new_c in 1:ca.C
                if new_c != state.conformations[i]
                    new_state = copy(state)
                    new_state.conformations[i] = new_c
                    add_new_rate(vertex, statetoi(new_state, ca))
                end
            end
        end
        # Add ligand (de)binding rates
        for i in 1:ca.N
            cur_occ = state.occupations[i]
            if cur_occ > 0
                new_state = copy(state)
                new_state.occupations[i] -= 1
                add_new_rate(vertex, statetoi(new_state, ca))
            end
            if cur_occ < ca.B
                new_state = copy(state)
                new_state.occupations[i] += 1
                add_new_rate(vertex, statetoi(new_state, ca))
            end
        end
    end
end

function add_edges_balanced_universtal_scale!(ca::ComplexAllosteryGM)
    universal_scale = Symbolics.variable(:r)
    gibbs_factors = calc_gibbs_factor.(allstates(ca), ca)

    function add_new_rate(v1, v2)
        if (v1 < v2)
            divp = gibbs_factors[v2] / gibbs_factors[v1]
            expfactor = sqrt(divp)
            add_edge!(ca.graph, v1, v2, simplify(universal_scale * expfactor))
            add_edge!(ca.graph, v2, v1, simplify(universal_scale / expfactor))
        end
    end

    for vertex in 1:numstates(ca)
        state = itostate(vertex, ca)
        # Add conformational change transitions
        for i in 1:ca.N
            for new_c in 1:ca.C
                if new_c != state.conformations[i]
                    new_state = copy(state)
                    new_state.conformations[i] = new_c
                    add_new_rate(vertex, statetoi(new_state, ca))
                end
            end
        end
        # Add ligand (de)binding rates
        for i in 1:ca.N
            cur_occ = state.occupations[i]
            if cur_occ > 0
                new_state = copy(state)
                new_state.occupations[i] -= 1
                add_new_rate(vertex, statetoi(new_state, ca))
            end
            if cur_occ < ca.B
                new_state = copy(state)
                new_state.occupations[i] += 1
                add_new_rate(vertex, statetoi(new_state, ca))
            end
        end
    end
end

function check_detailed_balance(ca::ComplexAllosteryGM; sim=false, fil=true)
    gfs = calc_gibbs_factor.(allstates(ca), ca)

    conditions = []
    for i in 1:numstates(ca)
        for j in i:numstates(ca)
            wij = get_weight(ca.graph, i, j)
            wji = get_weight(ca.graph, j, i)
            zero1 = iszero(wij)
            zero2 = iszero(wji)
            if !zero1 && !zero2
                push!(conditions, (gfs[i] * wij) - (gfs[j] * wji))
            elseif !(zero1 && zero2)
                @error f"detailed balance is not satisfied as some non-zero transitions do not have a reverse transition! (i, j)={(i, j)}"
            end
        end
    end

    if sim
        for i in eachindex(conditions)
            kaka = simplify(conditions[i]; rewriter=get_rewriter())
            if !iszero(kaka)
                println(kaka)
            end
        end
    end

    if fil
        conditions = filter(!iszero, conditions)
    end

    conditions
end

################################################################################
# "Concretizing" - plugging values for symbolic terms and making everything Float
################################################################################
function get_variables(arr::Array{Num})
    variables = Set{Num}()
    for term in arr
        for var in Symbolics.get_variables(term)
            push!(variables, var)
        end
    end
    variables
end
get_variables(::Nothing) = Set{Num}()
function get_variables(em::EnergyMatrices{Num})
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
Substitute some variables in a symbolic Array, EnergyMatrices or
ComplexAllosteryGM keeping it as symbolic.
"""
function substitute_partial(arr::AbstractArray, terms)
    new_arr = similar(arr)
    for i in eachindex(arr, new_arr)
        new_arr[i] = substitute(arr[i], terms)
    end
    new_arr
end
function substitute_partial(em::EnergyMatrices, terms)
    EnergyMatrices(substitute_partial(em.monomer, terms), substitute_partial(em.interactions, terms))
end
function substitute_partial(ca::ComplexAllosteryGM{S,F}, terms::Dict) where {S,F}
    new_energy_matrices = substitute_partial(ca.energy_matrices, terms)
    new_graph = SimpleWeightedDiGraph(
        substitute_partial(adjacency_matrix(ca.graph), terms)
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
        farr[i] = Float64(Symbolics.symbolic_to_float(substitute(arr[i], terms)))
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

################################################################################
# ca/graph manipulations
################################################################################
function filter_edges!(graph, threshold)
    for e in edges(graph)
        if weight(e) < threshold
            rem_edge!(graph, e)
        end
    end
end
function filter_edges!(ca::ComplexAllosteryGM{S,F}, args...) where {S,F}
    if F == Num
        throw(ArgumentError(f"function filter_edges! is only valid for ComplexAllosteryGM types with a concrete weight type - aka F != Num"))
    end
    filter_edges!(ca.graph, args...)
end

function keep_best_only!(graph, next_vertex_func=next_vertex_choose_max)
    for vertex in 1:nv(graph)
        best = next_vertex_func(graph, vertex)
        kaka = copy(neighbors(graph, vertex))
        for neighbor in kaka
            if neighbor != best
                rem_edge!(graph, vertex, neighbor)
            end
        end
    end
end
function keep_best_only!(ca::ComplexAllosteryGM{S,F}, args...) where {S,F}
    if F == Num
        throw(ArgumentError(f"function keep_best_only! is only valid for ComplexAllosteryGM types with a concrete weight type - aka F != Num"))
    end
    keep_best_only!(ca.graph, args...)
end

################################################################################
# ca/graph/matrix symbolics manipulation using Wolfram
################################################################################
function substitute_safenames(obj)
    terms = Dict{Num,Num}()
    letter = 'a'
    for var in get_variables(obj)
        terms[var] = Symbolics.variable(letter)
        letter += 1
        if letter > 'z'
            throw(ErrorException("object has too many variables, ran out of alphabet"))
        end
    end
    rev_terms = map(reverse, collect(terms))
    substitute_partial(obj, terms), rev_terms
end

function make_safe(safe, obj)
    if safe
        sobj, rt = substitute_safenames(obj)
        sobj, (srlst -> substitute_partial(srlst, rt))
    else
        obj, identity
    end
end

function w_simplify(obj; full=true, safe=false)
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

################################################################################
# Running, plotting and saving functions
################################################################################
@inline function p_get_adjmat_symsafe(ca::ComplexAllosteryGM{S,F}) where {S,F}
    if F == Num
        map(x -> if iszero(x)
                0.0
            else
                1.0
            end, adjacency_matrix(ca.graph))
    else
        adjacency_matrix(ca.graph)
    end
end

# Plotting, also used by simGM
function plotGM(ca::ComplexAllosteryGM{S,F}, args...;
    ax=nothing,
    layout=(:NRbs3,),
    jitter_devs=0.01,
    node_size_scale=10.0,
    interactive=:auto,
    colorbar=:auto,
    colormap=:dense,
    edges_filter=nothing,
    symca_labels=:auto,
    kwargs...
) where {S,F}
    if !isnothing(edges_filter)
        ca = copy(ca)
        filter_edges!(ca, edges_filter)
    end

    dim = nothing
    axis_labels = nothing
    if isa(layout, NetworkLayout.AbstractLayout)
        layout = layout(p_get_adjmat_symsafe(ca))
        dim = length(layout[1])
    else
        if isa(layout, Symbol)
            layout = (layout,)
        end
        (type, layout_args...) = layout
        add_jitter = false
        dim = 3 # mostly
        if type == :partlyc3d
            sub_layout, z_scale = layout_args
            layout = [
                (x, y, z * z_scale) for ((x, y), z) in zip(sub_layout(p_get_adjmat_symsafe(ca)), calc_numligands.(allstates(ca)))
            ]
        elseif type == :NEgf3
            z_scale = layout_args[1]
            layout = [
                (calc_numligands(st), calc_energy(st, ca), z_scale * log(calc_gibbs_factor(st, ca))) for st in allstates(ca)
            ]
        elseif type == :NEbs3
            layout = [
                (calc_numligands(st), calc_energy(st, ca), calc_numboundaries(st, ca)) for st in allstates(ca)
            ]
            maxs = maximum(layout)
            axis_labels = (f"N({maxs[1]})", f"E({maxs[2]})", f"boundaries({maxs[3]})")
            add_jitter = true
        elseif type == :NRbs3
            layout = [
                Point{3,Float64}(calc_numligands(st), calc_numofconf(st), calc_numboundaries(st, ca)) for st in allstates(ca)
            ]
            maxs = maximum(layout)
            axis_labels = (f"N({maxs[1]})", f"#relaxed({maxs[2]})", f"#boundaries({maxs[3]})")
            add_jitter = true
        elseif type == :NE2
            layout = [
                (Float64(calc_numligands(st)), calc_energy(st, ca)) for st in allstates(ca)
            ]
            maxs = maximum(layout)
            axis_labels = (f"N({maxs[1]})", f"E({maxs[2]})")
            dim = 2
            add_jitter = true
        elseif type == :Nbs2
            layout = [
                (Float64(calc_numligands(st)), Float64(calc_numboundaries(st, ca))) for st in allstates(ca)
            ]
            maxs = maximum(layout)
            axis_labels = (f"N({maxs[1]})", f"boundaries({maxs[2]})")
            dim = 2
            add_jitter = true
        elseif type == :NR2
            layout = [
                (Float64(calc_numligands(st)), Float64(calc_numofconf(st))) for st in allstates(ca)
            ]
            maxs = maximum(layout)
            axis_labels = (f"N({maxs[1]})", f"Number of relaxed({maxs[2]})")
            dim = 2
            add_jitter = true
        elseif type == :N1
            layout = [
                (Float64(calc_numligands(st)), 0.0) for st in allstates(ca)
            ]
            maxs = maximum(layout)
            axis_labels = (f"N({maxs[1]})", "")
            dim = 2
            add_jitter = true
        else
            throw(ArgumentError(f"layout of {layout} is not recognised"))
        end

        if type in [:NEgf3, :NEbs3]
            ranges = (maximum(layout) .- minimum(layout))
            layout = [pos ./ ranges for pos in layout]
        end

        if add_jitter
            if length(jitter_devs) == 1
                jitter_devs = Tuple(jitter_devs for _ in 1:length(layout[1]))
            end
            for i in eachindex(layout)
                offsets = Tuple(randn() * dev for dev in jitter_devs)
                layout[i] = layout[i] .+ offsets
            end
        end
    end

    auto_kwargs = Dict{Symbol,Any}()

    node_color = [:snow3 for _ in 1:numstates(ca)]
    if F != Num
        node_size = node_size_scale * sqrt(nv(ca.graph)) * calc_gibbs_factor.(allstates(ca), ca) / calc_partition_function(ca)
        auto_kwargs[:edge_color] = weight.(edges(ca.graph))
        auto_kwargs[:edge_attr] = (; colormap=colormap, colorrange=(0.0, maximum(weight.(edges(ca.graph)))))
        auto_kwargs[:arrow_attr] = auto_kwargs[:edge_attr]
    else
        node_size = [10.0 for _ in 1:numstates(ca)]
    end

    if symca_labels == :auto
        symca_labels = (F == Num) && (nv(ca.graph) < 50)
    end
    if symca_labels
        auto_kwargs[:nlabels] = repr.(allstates(ca))
        auto_kwargs[:elabels] = repr.(weight.(edges(ca.graph)))
        auto_kwargs[:elabels_shift] = 0.4
    end

    if isnothing(ax)
        (fig, ax, plot) = graphplot(ca.graph, args...;
            layout,
            node_color,
            node_size,
            auto_kwargs...,
            kwargs...
        )
    else
        fig = ax.parent
        plot = graphplot!(ax, ca.graph, args...;
            layout,
            node_color,
            node_size,
            auto_kwargs...,
            kwargs...
        )
    end

    if interactive == :auto
        interactive = (dim == 2)
    end
    if interactive
        make_interactive(ax, plot)
    end

    if !isnothing(axis_labels)
        if dim == 3
            xlabel!(ax.scene, axis_labels[1])
            ylabel!(ax.scene, axis_labels[2])
            zlabel!(ax.scene, axis_labels[3])
        elseif dim == 2
            ax.xlabel[] = axis_labels[1]
            ax.ylabel[] = axis_labels[2]
        end
    end

    if colorbar == :auto
        colorbar = (F != Num) && (dim == 3)
    end
    if colorbar
        edgecolorbar(fig, plot)
    end

    Makie.FigureAxisPlot(fig, ax, plot)
end

function plot_CA_sym(ca::ComplexAllosteryGM, args...; kwargs...)
    plotGM(ca, args...;
        nlabels=repr.(allstates(ca)),
        elabels=repr.(weight.(edges(ca.graph))),
        elabels_shift=0.4,
        kwargs...
    )
end

function eq_stats_plot(ca::ComplexAllosteryGM{S,Num}) where {S}
    fig = Figure()

    variables = Num.(get_variables(ca.energy_matrices))
    @variables μ, kT
    append!(variables, μ)

    sliders = [(label=repr(var), range=0.0:0.01:10.0, startvalue=1.0) for var in variables]
    sg = SliderGrid(fig[1, 1], sliders...)
    variable_observables = [x.value for x in sg.sliders]

    overall_stats = GridLayout()
    fig[1, 2] = overall_stats
    colsize!(fig.layout, 2, 30)

    partition_function_f = build_function(substitute(calc_partition_function(ca), (kT => 1)), variables; expression=Val(false))

    avg_numligands_f = build_function(substitute(calc_avg_numligands(ca), (kT => 1)), variables; expression=Val(false))
    avg_numligands_o = lift((args...) -> f"N={avg_numligands_f(args):.3g}", variable_observables...)
    Label(overall_stats[1, 1], avg_numligands_o)

    avg_energy_f = build_function(substitute(calc_avg_energy(ca), (kT => 1)), variables; expression=Val(false))
    avg_energy_o = lift((args...) -> f"E={avg_energy_f(args):.3g}", variable_observables...)
    Label(overall_stats[2, 1], avg_energy_o)

    # Do the barplot
    state_reordering = sortperm(calc_numligands.(allstates(ca)); alg=Base.Sort.DEFAULT_STABLE)
    state_labels = repr.(allstates(ca))[state_reordering]
    ax = Axis(fig[2, :];
        xticks=(1:numstates(ca), state_labels),
        xticklabelrotation=pi / 5
    )
    ylims!(ax, 0.0, 0.5)

    gfs = substitute.(calc_gibbs_factor.(allstates(ca), ca), (kT => 1.0))
    gfs_fs = build_function.(gfs, Ref(variables); expression=Val(false))

    probabilities = lift(variable_observables...) do (vars...)
        pfval = partition_function_f(vars)
        [gfs_fs[i](vars) / pfval for i in 1:numstates(ca)][state_reordering]
    end

    bp = barplot!(ax, 1:numstates(ca), probabilities)

    Makie.FigureAxisPlot(fig, ax, bp)
end

function eq_coopbinding_plot(ca::ComplexAllosteryGM{S,Num}) where {S}
    fig = Figure()

    energy_vars = collect(get_variables(ca.energy_matrices))
    @variables μ, εsol, kT, c

    slider_vars = [energy_vars; εsol]

    sliders = [(label=repr(var), range=-5.0:0.1:10.0, startvalue=1.0) for var in slider_vars]
    sg = SliderGrid(fig[1, 1], sliders...)
    slider_observables = [x.value for x in sg.sliders]

    ax = Axis(fig[2, :])

    prepped_expression = substitute(calc_avg_numligands(ca), Dict(
        kT => 1,
        μ => εsol + log(c)
    ))
    avg_numligands_f = build_function(prepped_expression, [energy_vars; εsol; c]; expression=Val(false))

    conc_range = LinRange(0, 200, 5000)

    avg_numligands = lift(slider_observables...) do (vars...)
        [avg_numligands_f((vars..., cval)) for cval in conc_range]
    end

    sc = scatter!(ax, conc_range, avg_numligands)

    Makie.FigureAxisPlot(fig, ax, sc)
end

# Random util
transition_matrix(ca::ComplexAllosteryGM) = Matrix(transpose(adjacency_matrix(ca.graph)))

# Saving internal data
function save_transmat(ca::ComplexAllosteryGM; name=savename("transmat", ca, "table"), short=false)
    file = open(plotsdir(name), "w")

    row_names = repr.(1:numstates(ca)) .* " / " .* repr.(allstates(ca))
    am = transmat(ca; mat=true)
    if short
        header = [""; repr.(1:numstates(ca))]
        am = map(x -> if iszero(x)
                0
            else
                1
            end, am)
    else
        header = [""; row_names]
    end
    modified_matrix = [row_names am]

    pretty_table(file, modified_matrix; header)
    close(file)
end

function save_gibbs_factors(ca::ComplexAllosteryGM; name=savename("adjmat", ca))
    file = open(plotsdir(savename("gibbsfs", ca, "table")), "w")
    data = [repr.(allstates(ca)) repr.(calc_gibbs_factor.(allstates(ca), ca))]
    pretty_table(file, data; header=["state", "gibbs factor"])
    close(file)
end

# Symbolics helper functions
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
