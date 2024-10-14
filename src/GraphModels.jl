using DrWatson
@quickactivate "TopoStochSim"

using Graphs, SimpleWeightedGraphs
using GLMakie
using StatsBase
using StaticArrays
using PyFormattedStrings


################################################################################
# Abstract type and API definitions
################################################################################
abstract type AbstractGraphModel end
function graph(gm::AbstractGraphModel)
    throw(ErrorException(f"No method of \"graph\" was provided for type \"{typeof(gm)}\""))
end
function plotGM(gm::AbstractGraphModel, args...; kwargs...)
    throw(ErrorException(f"No method of \"plotGM\" was provided for type \"{typeof(gm)}\""))
end

abstract type AbstractStatefulGraphModel{N,T<:Integer} <: AbstractGraphModel end
struct State{N,T<:Integer}
    extstate::Union{SVector{N,T},MVector{N,T}}
    substate::T
end
function states(gm::AbstractStatefulGraphModel)
    throw(ErrorException(f"No method of \"states\" was provided for type \"{typeof(gm)}\""))
end
function statetoi(gm::AbstractStatefulGraphModel{N,T}, state::State{N,T}) where {N,T}
    findfirst(x -> x == state, states(gm))
end
function itostate(gm::AbstractStatefulGraphModel, i)
    states(gm)[i]
end

################################################################################
# Running a graph model
################################################################################
function find_next_vertex(graph, v)
    options = neighbors(graph, v)
    sample(options, Weights([get_weight(graph, v, o) for o in options]))
end

function simGM(gm::AbstractGraphModel, steps; initial_vertex=nothing, delay=0.5, plot=true, print=false, kwargs...)
    vertex = isnothing(initial_vertex) ? rand(1:nv(graph(gm))) : initial_vertex

    if plot
        f, a, p = plotGM(gm; kwargs...)
        display(f)
        p.node_color[][vertex] = :red
        p.node_color = p.node_color[]
        p.node_size[][vertex] = 20
        p.node_size = p.node_size[]
    end

    try
        for _ in 1:steps
            sleep(delay)
            if plot
                p.node_color[][vertex] = :snow3
                p.node_size[][vertex] = 10
            end
            if print
                println(f"{vertex:<10d}")
                # println(f"{vertex:<10d} - {itostate(gm, vertex)}")
            end
            vertex = find_next_vertex(graph(gm), vertex)
            if plot
                p.node_color[][vertex] = :red
                p.node_color = p.node_color[]
                p.node_size[][vertex] = 20
                p.node_size = p.node_size[]
            end
        end
    catch InterruptException
        try
            GLMakie.closeall()
        catch TaskFailedException
        end
    end

    vertex
end

################################################################################
# Utility functions
################################################################################
function edgecolorbar((f, a, p))
    edgeploti = findfirst(typeof(plot) <: Plot{GraphMakie.edgeplot} for plot in p.plots)
    Colorbar(f[1, 2], p.plots[edgeploti])
end

function make_interactive((f, a, p))
    deregister_interaction!(a, :rectanglezoom)
    function node_drag_action(state, idx, event, axis)
        p[:node_pos][][idx] = event.data
        p[:node_pos][] = p[:node_pos][]
    end
    ndrag = NodeDragHandler(node_drag_action)
    register_interaction!(a, :ndrag, ndrag)
end
