module GraphModels

using Reexport

using LetterCodes

@reexport using DrWatson
@reexport using PyFormattedStrings, Printf
@reexport using StaticArrays, ElasticArrays, SparseArrays
@reexport using LinearAlgebra
@reexport using StatsBase
@reexport using Graphs, SimpleWeightedGraphs, NetworkLayout
@reexport using Makie, GraphMakie, Colors
@reexport using Symbolics
@reexport using Symbolics: variable, variables
@reexport using JLD2
@reexport import PlotUtils
using HDF5
import Dates

import Base: copy, broadcastable, show, display, convert
import Base: push!, getindex, length, iterate, view
copy(::Nothing) = nothing

################################################################################
# Core abstract types
################################################################################
"""
Pretty much synonymous to SimpleWeightedDiGraph but with a clear name.
Core API:
- graph
- numstates
- allstates
Fancy plotting API:
- p_named_layouts
- p_do_layout
- plotgm_layout_default
- plotgm_kwarg_defaults
Symbolics API (mainly for F <: Num):
- get_variables
- ssubstitute
- substitute_to_float
- make_factory
"""
abstract type AbstractGraphModel{F} end
function graph(gm::AbstractGraphModel)
    throw(ErrorException(f"No method of \"graph\" was provided for type \"{typeof(gm)}\""))
end

AbstractFloatGraph = AbstractSimpleWeightedGraph{T,<:AbstractFloat} where {T}
AbstractNumGraph = AbstractSimpleWeightedGraph{T,<:Num} where {T}
export AbstractGraphModel, AbstractFloatGraph, AbstractNumGraph

numstates(gm::AbstractGraphModel) = nv(graph(gm))
allstates(gm::AbstractGraphModel) = throw(ErrorException(f"No method of \"allstates\" was provided for type \"{typeof(gm)}\""))
export graph, numstates, allstates

################################################################################
# Any other types
################################################################################
struct Eigensystem{F}
    dim::Int
    evals::Vector{F}
    evecs::ElasticMatrix{F} # each column in an eigenvector, easily expandable
    function Eigensystem(evals::AbstractVector{F}, evecs::AbstractMatrix{F}, dim=nothing) where {F}
        if isnothing(dim)
            if length(evecs) < 1
                throw(ArgumentError(f"cannot infer dim from no eigenvectors, need to provide dim explicitly"))
            end
            dim = length(evecs[begin])
        end
        new{F}(dim, convert(Vector{F}, evals), convert(ElasticMatrix{F}, evecs))
    end
end
function Eigensystem(dim; numtype::Val{F}=Val(Float64)) where {F}
    Eigensystem(Vector{F}(undef, 0), ElasticMatrix{F}(undef, dim, 0), dim)
end
show(io::IO, es::Eigensystem{F}) where {F} = print(io, f"Eigensystem{{{F}}}(dim={es.dim}, {length(es.evals)} items)")
export Eigensystem

# somewhat general
function col(A::AbstractMatrix, i)
    @view A[:, i]
end
col(es::Eigensystem, i) = col(es.evecs, i)
function row(A::AbstractMatrix, i)
    @view A[i, :]
end
row(es::Eigensystem, i) = row(es.evecs, i)
evec(es::Eigensystem, i) = col(es, i)
export col, row, evec

# Eigensystem methods
function push!(es::Eigensystem{F}, (eval, evec)::Tuple{<:F,<:AbstractVector{F}}) where {F}
    if length(evec) != es.dim
        throw(ArgumentError(f"cannot append eigenvector of length {length(evec)} to eigensystem with dimension {es.dim}"))
    end
    push!(es.evals, eval)
    append!(es.evecs, evec)
end
getindex(es::Eigensystem, i) = (es.evals[i], collect(col(es, i)))
view(es::Eigensystem, i) = ((@view es.evals[i]), col(es, i))
length(es::Eigensystem) = length(es.evals)
function iterate(es::Eigensystem, i=1)
    if (1 <= i <= length(es.evals))
        ((es.evals[i], evec(es, i)), i + 1)
    else
        nothing
    end
end

subes(es::Eigensystem, is::AbstractVector) = Eigensystem(es.evals[is], col(es, is))
export subes

################################################################################

include("util.jl")     # Has bits of GM API
include("graphs_extra.jl")     # Has more custom pure graph related functions
include("manipulations.jl")
include("manipulations_sccs.jl")
include("symbolics.jl")     # Has bits of GM API
include("plots_and_io.jl")            # Has bits of GM API

end # module GraphModels
