using DrWatson
@quickactivate "TopoStochSim"

using Graphs, SimpleWeightedGraphs
using GLMakie
using StatsBase
using StaticArrays
using PyFormattedStrings

include("general.jl")


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
function next_vertex_random_proportional(graph, v)
    options = neighbors(graph, v)
    sample(options, Weights([get_weight(graph, v, o) for o in options]))
end

function next_vertex_choose_max(graph, v)
    options = neighbors(graph, v)
    weights = [get_weight(graph, v, o) for o in options]
    findmax(weights)[2]
end

function simGM(gm::AbstractGraphModel, steps;
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
        f, a, p = plotGM(gm; kwargs...)
        display(f)
        for vertex in vertices
            p.node_color[][vertex] = :red
            p.node_color = p.node_color[]
            p.node_size[][vertex] = 20
            p.node_size = p.node_size[]
        end
    end

    try
        for _ in 1:steps
            sleep(delay)
            if plot
                for vertex in vertices
                    p.node_color[][vertex] = :snow3
                    p.node_size[][vertex] = 10
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
                    p.node_color[][vertex] = :red
                    p.node_color = p.node_color[]
                    p.node_size[][vertex] = 20
                    p.node_size = p.node_size[]
                end
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
