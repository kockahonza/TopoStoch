using DrWatson
@quickactivate "TopoStochSim"

using Revise

using PyFormattedStrings

using LinearAlgebra, StatsBase, StaticArrays, SparseArrays
using Graphs, SimpleWeightedGraphs, NetworkLayout
using GLMakie, GraphMakie, Colors
using Symbolics, MathLink, SymbolicsMathLink
using JLD2

import Dates
import PlotUtils

import Base: copy, broadcastable, show, display, convert
copy(::Nothing) = nothing

timestamp() = Dates.format(Dates.now(), "yymmdd_HMS")

################################################################################
# Base abstract types
################################################################################
"""
Pretty much synonymous to SimpleWeightedDiGraph but with a clear name.
The API is: graph(gm), plotgm(gm) and optionally allstates for all and for
symbolic gms (F <: Num) also get_variables, ssubstitute, substitute_to_float
and make_factory as defined in gm_symbolics.jl.
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

includet("gm_symbolics.jl")     # Has bits of GM API
includet("gm_manipulations.jl")
includet("gm_io.jl")            # Has bits of GM API

################################################################################
# Running a graph model
################################################################################
function next_vertex_random_proportional(graph, v)
    options = neighbors(graph, v)
    sample(options, Weights([get_weight(graph, v, o) for o in options]))
end

function next_vertex_choose_max(graph, v)
    options = neighbors(graph, v)
    weights = [get_weight(graph, v, o) for o in options]
    options[findmax(weights)[2]]
end

# TODO: This could use an update and a change of how highlighting is done
function simgm(gm::AbstractGraphModel{<:AbstractFloat}, steps;
    initial_vertices=nothing,
    num_vertices=1, next_vertex_func=next_vertex_random_proportional, delay=0.5,
    plot=true,
    print=false,
    kwargs...
)
    vertex = isnothing(initial_vertices) ? rand(1:nv(graph(gm))) : initial_vertices

    if isnothing(initial_vertices)
        vertices = [rand(1:nv(graph(gm))) for _ in 1:num_vertices]
    end

    if plot
        f, a, p = plotgm(gm; kwargs...)
        display(f)
        base_node_color = copy(p.node_color[])
        base_node_size = copy(p.node_size[])
        function highlight_node(vertex)
            p.node_color[][vertex] = :red
            p.node_size[][vertex] = 15
        end
        function unhighlight_node(vertex)
            p.node_color[][vertex] = base_node_color[vertex]
            p.node_size[][vertex] = base_node_size[vertex]
        end
        function update_observables()
            p.node_color = p.node_color[]
            p.node_size = p.node_size[]
        end
        for vertex in vertices
            highlight_node(vertex)
        end
        update_observables()
    end

    try
        for _ in 1:steps
            sleep(delay)
            if plot
                for vertex in vertices
                    unhighlight_node(vertex)
                end
            end
            if print
                for vertex in vertices
                    println(f"{vertex:<10d}")
                end
            end
            for i in eachindex(vertices)
                vertices[i] = next_vertex_func(graph(gm), vertices[i])
            end
            if plot
                for vertex in vertices
                    highlight_node(vertex)
                end
                update_observables()
            end
        end
    catch InterruptException
        try
            GLMakie.closeall()
        catch TaskFailedException
        end
    end

    vertices
end
