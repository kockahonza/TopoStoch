using DrWatson
@quickactivate "TopoStochSim"

using Graphs, SimpleWeightedGraphs
using GLMakie
using StatsBase
using StaticArrays, SparseArrays
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
    options[findmax(weights)[2]]
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

################################################################################
# Custom graph analysis
################################################################################
# TODO: FIX: This is not finished!!
function find_cycles(graph, next_vertex_func=next_vertex_choose_max)
    cycles = []
    cycle_map = Vector{Union{Nothing,Int}}(nothing, nv(graph))
    for vertex in 1:nv(graph)
        if isnothing(cycle_map[vertex])
            cycle = []
            true_cycle_start_i = findfirst(x -> x == vertex, cycle)
            while isnothing(true_cycle_start_i)
                push!(cycle, vertex)
                vertex = next_vertex_func(graph, vertex)
                true_cycle_start_i = findfirst(x -> x == vertex, cycle)
            end
            for v in cycle
                cycle_map[v] = length(cycles) + 1
            end
            cycle = cycle[true_cycle_start_i:end]
            push!(cycles, cycle)
        end
    end
    (cycles, cycle_map)
end

################################################################################
# Utility functions
################################################################################
"""
The standard Markov chain/Master equation transition matrix where the
ij th element corresponds to a transition/rate from j to i.
"""
function transmat(gm::AbstractGraphModel; mat=false)
    tm = transpose(adjacency_matrix(graph(gm)))
    mat ? Matrix(tm) : sparse(tm)
end

# Plotting
function edgecolorbar(fig, plot)
    edgeploti = findfirst(typeof(plot_) <: Plot{GraphMakie.edgeplot} for plot_ in plot.plots)
    Colorbar(fig[1, 2], plot.plots[edgeploti])
end
edgecolorbar((fig, ax, plot)) = edgecolorbar(fig, plot)

function make_interactive(ax, plot)
    deregister_interaction!(ax, :rectanglezoom)
    function node_drag_action(state, idx, event, axis)
        plot[:node_pos][][idx] = event.data
        plot[:node_pos][] = plot[:node_pos][]
    end
    ndrag = NodeDragHandler(node_drag_action)
    register_interaction!(ax, :ndrag, ndrag)
end
make_interactive((_, ax, plot)) = make_interactive(ax, plot)
