using DrWatson
@quickactivate "TopoStochSim"

include(srcdir("graphmodels.jl"))

using Graphs, MetaGraphsNext, SimpleWeightedGraphs
using GLMakie, GraphMakie, NetworkLayout
using StatsBase
using Symbolics
using StaticArrays, SparseArrays, DataFrames
using PrettyTables, Latexify
using PyFormattedStrings
using Base.Threads

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

    environment::Union{Nothing,Tuple{F,F}} # Pair of μ and kT
    metadata
    function ComplexAllosteryGM(N::T, C::T, B::T;
        symmetry::S=Chain(),
        numtype::Val{F}=Val(Num),
        energy_matrices=nothing,
        graph=nothing,
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

        new{S,F,T}(N, C, B, energy_matrices, graph, environment, metadata)
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
    println(f"\ngraph:")
    show(io, mime, ca.graph)
    println(f"\nenvironment:")
    show(io, mime, ca.environment)
    if !isnothing(ca.metadata)
        println(f"\nmetadata:")
        show(io, mime, ca.metadata)
    end
end
graph(ca::ComplexAllosteryGM) = ca.graph

# Utils state related functions
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

function validEM(em::EnergyMatrices, C, B)
    (size(em.monomer) == (C, B + 1)) && (size(em.interactions) == (C, C))
end
validEM(em::EnergyMatrices, ca::ComplexAllosteryGM) = validEM(em, ca.C, ca.B)

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

function calc_neighbors_energy(st::CAState, i, em::EnergyMatrices, s::Symmetry)
    (left, right) = get_neighbors(st, i, s)
    energy = 0
    if !isnothing(left)
        energy += em.interactions[left, st.conformations[i]]
    end
    if !isnothing(right)
        energy += em.interactions[st.conformations[i], right]
    end
    energy
end
calc_neighbors_energy(st::CAState, i, ca::ComplexAllosteryGM{S}) where {S} = calc_neighbors_energy(st, i, ca.energy_matrices, S())

function calc_energy(st::CAState, em::EnergyMatrices, s::Symmetry)
    energy = 0
    for i in 1:length(st.conformations)
        energy += em.monomer[st.conformations[i], st.occupations[i]+1]
        energy += 0.5 * calc_neighbors_energy(st, i, em, s)
    end
    energy
end
function calc_energy(st::CAState, ca::ComplexAllosteryGM{S}) where {S}
    calc_energy(st, ca.energy_matrices, S())
end

calc_numligands(st::CAState) = sum(st.occupations)

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
    # Just use a constant for all the rates,disregarding detailed balance
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
        for j in 1:numstates(ca)
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
        # FIX: This is not 100% working idk why but not work the time...
        @threads for i in eachindex(conditions)
            kaka = simplify(conditions[i]; rewriter=get_rewriter())
            if !iszero(kaka)
                println(kaka)
            end
        end
    end

    if fil
        conditions = filter(!iszero, conditions)
    end

    nothing
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


function substitute_to_float_(arr, farr, terms::Dict{Num,Float64})
    for i in eachindex(arr, farr)
        farr[i] = Float64(Symbolics.symbolic_to_float(substitute(arr[i], terms)))
    end
end
function substitute_to_float(arr::Array, terms::Dict{Num,Float64})
    farr = Array{Float64}(undef, size(arr))
    substitute_to_float_(arr, farr, terms)
    farr
end
function substitute_to_float(arr::AbstractSparseMatrix, terms::Dict{Num,Float64})
    farr = similar(arr, Float64)
    substitute_to_float_(arr, farr, terms)
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
    ComplexAllosteryGM(ca.N, ca.C, ca.B;
        symmetry=S(),
        numtype=Val(Float64),
        energy_matrices=new_energy_matrices,
        graph=new_graph,
        environment=(terms[Symbolics.variable(:μ)], terms[Symbolics.variable(:kT)]),
        metadata=terms
    )
end

################################################################################
# The simple case, C=2 and B mostly 1
################################################################################
function make_simple(N, B; edge_t=:EL1)
    ca = ComplexAllosteryGM(N, 2, B; energy_matrices=make_EM_sym_C2(B))

    if edge_t == :EL1
        add_edges_simple_EL1!(ca)
    elseif edge_t == :EL2
        add_edges_simple_EL2!(ca)
    else
        throw(ArgumentError(f"edge_t \"{edge_t}\" not recognised"))
    end

    ca
end

"""
Energy landscape inspired symbolic edges that use energy height barriers
and and some physicsy arguments done on pen and paper. Maybe works?
"""
function add_edges_simple_EL1!(ca::ComplexAllosteryGM)
    rate_occ_change = Symbolics.variable(:r_B)
    rate_conf_change = Symbolics.variable(:r_C)
    @variables ε_t, Δε_r, ε_b, μ, kT
    fact_1_dec = exp(-(-ε_t + μ) / kT)
    fact_2_inc_0 = exp(-(-Δε_r) / kT)
    fact_2_inc = exp(-(Δε_r) / kT)
    fact_2_dec = exp(-(-ε_t + μ + Δε_r) / kT)

    for vertex in 1:numstates(ca)
        state = itostate(vertex, ca)
        # Add conformational change transitions
        for i in 1:ca.N
            for new_c in 1:ca.C
                if new_c != state.conformations[i]
                    new_state = copy(state)
                    new_state.conformations[i] = new_c
                    exponent_term = Num(0)
                    if state.occupations[i] != 0
                        exponent_term += ε_t
                        if state.conformations[i] == 2
                            exponent_term -= Δε_r
                        end
                    else
                        if state.conformations[i] == 2
                            exponent_term += Δε_r
                        end
                    end
                    exponent_term += calc_neighbors_energy(state, i, ca)
                    rate = rate_conf_change / exp(-exponent_term / kT)
                    add_edge!(
                        ca.graph,
                        vertex,
                        statetoi(new_state, ca),
                        rate
                    )
                end
            end
        end

        # Add ligand (de)binding rates
        for i in 1:ca.N
            cur_occ = state.occupations[i]
            cur_con = state.conformations[i]
            if cur_con == 1
                if cur_occ < ca.B
                    new_state = copy(state)
                    new_state.occupations[i] += 1
                    add_edge!(
                        ca.graph,
                        vertex,
                        statetoi(new_state, ca),
                        rate_occ_change
                    )
                end
                if cur_occ > 0
                    new_state = copy(state)
                    new_state.occupations[i] -= 1
                    add_edge!(
                        ca.graph,
                        vertex,
                        statetoi(new_state, ca),
                        rate_occ_change * fact_1_dec
                    )
                end
            elseif cur_con == 2
                if cur_occ == 0
                    new_state = copy(state)
                    new_state.occupations[i] += 1
                    add_edge!(
                        ca.graph,
                        vertex,
                        statetoi(new_state, ca),
                        rate_occ_change * fact_2_inc_0
                    )
                elseif cur_occ < ca.B
                    new_state = copy(state)
                    new_state.occupations[i] += 1
                    add_edge!(
                        ca.graph,
                        vertex,
                        statetoi(new_state, ca),
                        rate_occ_change * fact_2_inc
                    )
                end
                if cur_occ > 0
                    new_state = copy(state)
                    new_state.occupations[i] -= 1
                    add_edge!(
                        ca.graph,
                        vertex,
                        statetoi(new_state, ca),
                        rate_occ_change * fact_2_dec
                    )
                end
            else
                throw(ErrorException("kaka"))
            end
        end
    end
end

"""
Energy landscape inspired symbolic edges that use energy height barriers
and the actual gibbs factors for getting the rates. This however could very
well lead to rates > 1. Hence not ready!
"""
function add_edges_simple_EL2!(ca::ComplexAllosteryGM)
    gfs = calc_gibbs_factor.(allstates(ca), ca)

    # Label occupancy changes using r_(ci)+-
    rate_occ_change = Symbolics.variable(:r)
    @variables ε_t, Δε_r, ε_b, μ, kT

    # Label conformational changes via increasing indices or ri
    rates_i = 1
    function add_new_indexed_rate(v1, v2)
        if (v1 < v2)
            rate = Symbolics.variable(:ri, rates_i)
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
                    # add_new_indexed_rate(vertex, statetoi(new_state, ca))
                end
            end
        end

        # Add ligand (de)binding rates
        for i in 1:ca.N
            cur_occ = state.occupations[i]
            if cur_occ < ca.B
                new_state = copy(state)
                new_state.occupations[i] += 1
                rate = rate_occ_change * exp(-cur_occ * (ε_t - μ) / kT) / gfs[vertex]
                add_edge!(
                    ca.graph,
                    vertex,
                    statetoi(new_state, ca),
                    simplify(rate; rewriter=get_rewriter())
                )
            end
            if cur_occ > 0
                new_state = copy(state)
                new_state.occupations[i] -= 1
                rate = rate_occ_change * exp(-(cur_occ - 1) * (ε_t - μ) / kT) / gfs[vertex]
                add_edge!(
                    ca.graph,
                    vertex,
                    statetoi(new_state, ca),
                    simplify(rate; rewriter=get_rewriter())
                )
            end
        end
    end
end

us_calc_numtense(st::CAState) = count(x -> x == 1, st.conformations)
us_calc_numrelaxed(st::CAState) = count(x -> x == 2, st.conformations)

function terms_simple(mu, et, der, eb, rB, rC, kT_=1.0)
    @variables μ, ε_t, Δε_r, ε_b, r_B, r_C, kT
    Dict(
        μ => mu,
        ε_t => et,
        Δε_r => der,
        ε_b => eb,
        r_B => rB,
        r_C => rC,
        kT => kT_
    )
end
terms_simple0() = terms_simple(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

function terms_simple1(et, eb_option, mu_option)
    eb = if eb_option == 1
        0.0
    elseif eb_option == 2
        0.5 * et
    else
        throw(ArgumentError("eb_option must be 1 or 2"))
    end
    mu = if mu_option == 1
        0.3 * et
    elseif mu_option == 2
        0.8 * et
    elseif mu_option == 3
        1.2 * et
    else
        throw(ArgumentError("mu_option must be 1, 2 or 3"))
    end
    terms_simple(mu, et, 0.2 * et, eb, exp(-et), exp(-2.5 * et), 1.0)
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
    layout=(:c1, Shell(), 1.0),
    jitter_devs=0.01,
    node_size_scale=10.0,
    interactive=:auto,
    colorbar=:auto,
    colormap=:dense,
    edges_filter=nothing,
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
        (type, layout_args...) = layout
        add_jitter = false
        dim = 3 # mostly
        if type == :partlyc3d
            sub_layout, z_scale = layout_args
            layout = [
                (x, y, z * z_scale) for ((x, y), z) in zip(sub_layout(p_get_adjmat_symsafe(ca)), calc_numligands.(allstates(ca)))
            ]
        elseif type == :c3
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
                Point{3,Float64}(calc_numligands(st), us_calc_numrelaxed(st), calc_numboundaries(st, ca)) for st in allstates(ca)
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
                (Float64(calc_numligands(st)), Float64(us_calc_numrelaxed(st))) for st in allstates(ca)
            ]
            maxs = maximum(layout)
            axis_labels = (f"N({maxs[1]})", f"Number of relaxed({maxs[2]})")
            dim = 2
            add_jitter = true
        else
            throw(ArgumentError(f"layout of {layout} is not recognised"))
        end

        if type in [:c3, :c4]
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

    if isnothing(layout)
    elseif isa(layout, NetworkLayout.AbstractLayout)
        layout = layout(p_get_adjmat_symsafe(ca))
    end

    node_color = [:snow3 for _ in 1:numstates(ca)]
    if F != Num
        node_size = node_size_scale * sqrt(nv(ca.graph)) * calc_gibbs_factor.(allstates(ca), ca) / calc_partition_function(ca)
        edge_color = weight.(edges(ca.graph))
    else
        node_size = [10.0 for _ in 1:numstates(ca)]
        edge_color = [0.0 for _ in 1:ne(ca.graph)]
    end

    edge_attr = arrow_attr = (; colormap=colormap, colorrange=(0.0, maximum(weight.(edges(ca.graph)))))

    fap = graphplot(ca.graph, args...;
        layout,
        node_color,
        node_size,
        edge_color,
        edge_attr,
        arrow_attr,
        kwargs...
    )

    if interactive == :auto
        interactive = (dim == 2)
    end
    if interactive
        make_interactive(fap)
    end

    if !isnothing(axis_labels)
        if dim == 3
            xlabel!(fap.axis.scene, axis_labels[1])
            ylabel!(fap.axis.scene, axis_labels[2])
            zlabel!(fap.axis.scene, axis_labels[3])
        elseif dim == 2
            fap.axis.xlabel[] = axis_labels[1]
            fap.axis.ylabel[] = axis_labels[2]
        end
    end

    if colorbar == :auto
        colorbar = (F != Num) && (dim == 3)
    end
    if colorbar
        edgecolorbar(fap)
    end

    fap
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

# Saving internal data
function save_adj_matrix(ca::ComplexAllosteryGM; name=savename("adjmat", ca), short=false)
    file = open(plotsdir(savename("adjmat", ca, "table")), "w")

    row_names = repr.(1:numstates(ca)) .* " / " .* repr.(allstates(ca))
    am = Matrix{Any}(adjacency_matrix(ca.graph))
    if short
        header = [""; repr.(1:numstates(ca))]
        am = map(x -> if iszero(x)
                0
            else
                1
            end, am)
    else
        header = [""; repr.(allstates(ca))]
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