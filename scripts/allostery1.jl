using DrWatson
@quickactivate "TopoStochSim"

using Graphs, SimpleWeightedGraphs
# using GLMakie, GraphMakie, NetworkLayout
using PyFormattedStrings
using StaticArrays
using Symbolics

abstract type AbstractSymmetryType end
struct Chain <: AbstractSymmetryType end
struct Loop <: AbstractSymmetryType end

get_mvect_field_type(::MVector{N,T}) where {N,T} = T

"""
N is the number of monomers
C is the number of possible conformational states per monomer
B is the number of ligand binding sites per monomer
"""
struct AState{N,C,B,S<:AbstractSymmetryType,T<:Integer}
    conformations::MVector{N,T}
    occupations::MVector{N,T}
    function AState(N, C, B; symmetry::S=Chain(), conformations=nothing, occupations=nothing) where {S}
        if N < 1
            throw(ArgumentError("N must be greater than 0"))
        end
        if C < 1
            throw(ArgumentError("C must be greater than 0"))
        end
        if B < 0
            throw(ArgumentError("B must at least 0"))
        end
        if isnothing(conformations)
            conformations = @MVector fill(1, N)
        end
        if isnothing(occupations)
            occupations = @MVector fill(0, N)
        end
        field_type = get_mvect_field_type(conformations)
        if field_type != get_mvect_field_type(occupations)
            throw(ArgumentError("conformations and occupations mvectors have to contain elements of the same type"))
        end
        new{N,C,B,S,field_type}(conformations, occupations)
    end
end
conformations((; conformations, occupations)::AState) = conformations
occupations((; conformations, occupations)::AState) = occupations
function Base.getindex((; conformations, occupations)::AState{N,C,B}, i) where {N,C,B}
    (conformations[i], occupations[i])
end
function setconformation((; conformations, occupations)::AState{N,C,B}, conf, i) where {N,C,B}
    if (conf < 1) || (conf > C)
        throw(ArgumentError(f"only conformations between 1 and C={C} (inclusive) are allowed"))
    end
    conformations[i] = conf
end
function setoccupation((; conformations, occupations)::AState{N,C,B}, occ, i) where {N,C,B}
    if (occ < 0) || (occ > B)
        throw(ArgumentError(f"only occupations between 0 and B={B} (inclusive) are allowed"))
    end
    occupations[i] = occ
end
function Base.setindex!(astate::AState{N,C,B}, val, i) where {N,C,B}
    if length(val) != 2
        throw(ArgumentError("setindex! on AState can only take containers of lenght 2 to set both conformation and occupation"))
    end
    setconformation(astate, val[1], i)
    setoccupation(astate, val[2], i)
end

function Base.reverse((; conformations, occupations)::AState{N,C,B,S})::AState{N,C,B,S} where {N,C,B,S}
    AState(N, C, B; symmetry=S(), conformations=reverse(conformations), occupations=reverse(occupations))
end

function same_under_symmetry(as1::AState{N,C,B,Chain}, as2::AState{N,C,B,Chain}) where {N,C,B}
    (as1 == as2) || (as1 == reverse(as2))
end
