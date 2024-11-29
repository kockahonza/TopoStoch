module ComplexAllostery

using Reexport

@reexport using GraphModels # This reexports many of its dependencies

using GLMakie
using PrettyTables

import GraphModels: graph, numstates, allstates
import Base: copy, broadcastable, show, display

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
export CAState

abstract type Symmetry end
struct Chain <: Symmetry end
struct Loop <: Symmetry end
broadcastable(s::Symmetry) = Ref(s)
export Chain, Loop
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
export EnergyMatrices

mutable struct ComplexAllosteryGM{S<:Symmetry,F<:Number,T<:Integer} <: AbstractGraphModel{F}
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
        energy_matrices=copy(ca.energy_matrices),
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
export ComplexAllosteryGM

# Small utility functions
function validEM(em::EnergyMatrices, C, B)
    (size(em.monomer) == (C, B + 1)) && (size(em.interactions) == (C, C))
end
validEM(em::EnergyMatrices, ca::ComplexAllosteryGM) = validEM(em, ca.C, ca.B)
export validEM

function validstate(st::CAState, N, C, B)
    validperm(st.conformations, C, N) && validperm(st.occupations, B + 1, N; base=0)
end
validstate(st::CAState, ca::ComplexAllosteryGM) = validstate(st, ca.N, ca.C, ca.B)
export validstate

function statetoi(st::CAState, N, C, B)
    permtoi(st.conformations, C) + (permtoi(st.occupations, B + 1; base=0) - 1) * (C^N)
end
statetoi(st::CAState, ca::ComplexAllosteryGM) = statetoi(st, ca.N, ca.C, ca.B)
export statetoi

function itostate(i, N, C, B)
    (iocc, iconf) = divrem(i - 1, C^N)
    CAState(itoperm(iconf + 1, C, N), itoperm(iocc + 1, B + 1, N; base=0))
end
itostate(i, ca::ComplexAllosteryGM) = itostate(i, ca.N, ca.C, ca.B)
export itostate

function numstates(N, C, B)
    (C * (B + 1))^N
end
numstates(ca::ComplexAllosteryGM) = numstates(ca.N, ca.C, ca.B)

allstates_complexallostery(N, C, B) = (itostate(i, N, C, B) for i in 1:numstates(N, C, B))
allstates(ca::ComplexAllosteryGM) = allstates_complexallostery(ca.N, ca.C, ca.B)

labelstate_NT(st::CAState) = f"{calc_numligands(st)}, {calc_numofconf(st)}"
labelstate_NR(st::CAState) = f"{calc_numligands(st)}, {calc_numofconf(st, 2)}"
export labelstate_NT, labelstate_NR

################################################################################
# Include all the other bits
################################################################################
include("physics.jl")
include("symbolics.jl")
include("plots_and_io.jl")

include("Model2_5_C2.jl")

end # module ComplexAllostery
