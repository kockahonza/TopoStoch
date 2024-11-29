module GraphModels

using Reexport

@reexport using PyFormattedStrings

@reexport using DrWatson
@reexport using LinearAlgebra, StatsBase, StaticArrays, SparseArrays
@reexport using Graphs, SimpleWeightedGraphs, NetworkLayout
@reexport using Makie, GraphMakie, Colors
@reexport using Symbolics, MathLink, SymbolicsMathLink
@reexport using JLD2

import Dates
import PlotUtils

import Base: copy, broadcastable, display, convert
copy(::Nothing) = nothing

# Exports
export AbstractGraphModel, AbstractFloatGraph, AbstractNumGraph
export graph, numstates, allstates

################################################################################
# Core abstract types
################################################################################
"""
Pretty much synonymous to SimpleWeightedDiGraph but with a clear name.
Core API:
- graph
- numstates
- allstates
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

################################################################################

numstates(gm::AbstractGraphModel) = nv(graph(gm))
allstates(gm::AbstractGraphModel) = throw(ErrorException(f"No method of \"allstates\" was provided for type \"{typeof(gm)}\""))

include("util.jl")     # Has bits of GM API
include("symbolics.jl")     # Has bits of GM API
include("manipulations.jl")
include("plots_and_io.jl")            # Has bits of GM API
include("BasicSim.jl")

end # module GraphModels
