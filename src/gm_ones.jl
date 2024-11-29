using DrWatson
@quickactivate "TopoStochSim"

using Revise

using GraphModels

import GraphModels: graph, numstates, allstates
import GraphModels: p_named_layouts, p_do_layout, plotgm_kwarg_defaults
import Base: copy, broadcastable, show

################################################################################
# Base setup
################################################################################
abstract type Symmetry end
struct Chain <: Symmetry end
struct Loop <: Symmetry end
broadcastable(s::Symmetry) = Ref(s)

struct OnesGM{S<:Symmetry,F} <: AbstractGraphModel{F}
    N::Int
    graph::SimpleWeightedDiGraph{Int,F}
    metadata::Union{Nothing,Dict}
    function OnesGM(N;
        symmetry::S=Loop(),
        numtype::Val{F}=Val(Float64),
        graph=nothing,
        metadata=nothing
    ) where {S,F}
        if isnothing(graph)
            graph = SimpleWeightedDiGraph{Int,F}(numstates(N))
        elseif nv(graph) != numstates(N)
            throw(ArgumentError("invalid graph passed, does not have the correct number of nodes"))
        end

        new{S,F}(N, graph, metadata)
    end
end
function copy(ogm::OnesGM{S,F}) where {S,F}
    OnesGM(ogm.N; symmetry=S(), numtype=Val(F), graph=copy(ogm.graph), metadata=copy(ogm.metadata))
end
function show(io::IO, mime::MIME"text/plain", ogm::OnesGM{S,F}) where {S,F}
    println(f"OnesGM{{{S},{F}}}(N={ogm.N})")
    print(f"graph: ")
    show(io, mime, ogm.graph)
    if !isnothing(ogm.metadata)
        print(f"\nmetadata: ")
        show(io, mime, ogm.metadata)
    end
end
graph(ogm::OnesGM) = ogm.graph

numstates_ones(N) = 2^N
numstates(ogm::OnesGM) = numstates_ones(ogm.N)

itostate(i, N) = digits(i - 1; base=2, pad=N)
itostate(i, ogm::OnesGM) = itostate(i, ogm.N)

statetoi(state) = 1 + sum(state[i] * 2^(i - 1) for i in 1:length(state))

allstates_ones(N) = itostate.(1:numstates_ones(N), N)
allstates(ogm::OnesGM) = allstates_ones(ogm.N)

function findall_subvects(_::Chain, vect::Vector, subvect::Vector)
    subvect_len = length(subvect)
    matches = []
    for i in 1:(length(vect)-subvect_len+1)
        if (@view vect[i:i+subvect_len-1]) == subvect
            push!(matches, i)
        end
    end
    matches
end
function findall_subvects(_::Loop, vect::Vector, subvect::Vector)
    vect_len = length(vect)
    subvect_len = length(subvect)
    matches = []
    for i in 1:vect_len
        if (@view vect[mod1.(i:i+subvect_len-1, vect_len)]) == subvect
            push!(matches, i)
        end
    end
    matches
end

function calc_numboundaries(_::Chain, vect::Vector)
    vect_len = length(vect)
    numboundaries = 0
    for i in 1:(vect_len-1)
        if vect[i] != vect[i+1]
            numboundaries += 1
        end
    end
    numboundaries
end
function calc_numboundaries(_::Loop, vect::Vector)
    numboundaries = calc_numboundaries(Chain(), vect)
    if vect[end] != vect[begin]
        numboundaries += 1
    end
    numboundaries
end

function subin!(vect::Vector, i, subvect::Vector)
    for j in 1:length(subvect)
        vect[mod1(i + j - 1, length(vect))] = subvect[j]
    end
end

################################################################################
# Plotting
################################################################################
function p_named_layouts(ogm::OnesGM{S}, layout_name, layout_args) where {S}
    try
        invoke(p_named_layouts, Tuple{supertype(typeof(ogm)),Any,Any}, ogm, layout_name, layout_args)
    catch ArgumentError
        def_roffsets = false
        axis_labels = nothing

        if layout_name == :tree
            dim = 2
            states = allstates(ogm)
            groups = [[] for _ in 1:ogm.N+1]
            for i in 1:numstates(ogm)
                push!(groups[1+count(x -> x == 1, states[i])], i)
            end
            for g in groups
                sort!(g; by=i -> states[i], rev=true)
            end
            layout = Vector{Any}(undef, numstates(ogm))
            for (gi, group) in enumerate(groups)
                mid = length(group) / 2
                for (i, reali) in enumerate(group)
                    layout[reali] = (i - mid, float(gi))
                end
            end

            axis_labels = ("sort", "num ones")
        elseif layout_name == :treeB
            dim = 3
            states = allstates(ogm)
            groups = [[] for _ in 1:ogm.N+1]
            for i in 1:numstates(ogm)
                push!(groups[1+count(x -> x == 1, states[i])], i)
            end
            for g in groups
                sort!(g; by=i -> states[i], rev=true)
            end
            layout = Vector{Any}(undef, numstates(ogm))
            for (gi, group) in enumerate(groups)
                mid = length(group) / 2
                for (i, reali) in enumerate(group)
                    layout[reali] = (i - mid, float(gi), calc_numboundaries(S(), states[reali]))
                end
            end

            axis_labels = ("sort", "num ones", "num boundaries")
        elseif layout_name == :hamc
            dim = 3
            if ogm.N == 3
                layout = allstates(ogm)
            elseif ogm.N == 4
                layout = [(x + off, y + off, z + 1.5 * off) for (x, y, z, off) in allstates(ogm)]
            else
                throw(ArgumentError(":hamc layout can only be used with N=3"))
            end
        else
            rethrow()
        end

        dim, layout, def_roffsets, axis_labels
    end
end
function p_do_layout(ogm::OnesGM, layout=nothing, roffset_devs=nothing)
    if isnothing(layout)
        layout = :tree
    end
    invoke(p_do_layout, Tuple{AbstractGraphModel,Any,Any}, ogm, layout, roffset_devs)
end

plotgm_kwarg_defaults(_::OnesGM) = (;
    fnlabels=:repr,
    felabels=false,
    n_size=10.0,
    e_color=:dcurrents
)

"""
Alternative plotting option which makes low rate edges partially transparent
"""
function plot_ogm_min(args...; ecutoff=1.1, amin=0.2, amax=1.0, kwargs...)
    plotgm(args...;
        layout=makeprefunclayout(Spring(), filter_edges!, ecutoff),
        felabels=false,
        c_colormap=make_linalpha_cmap(:dense; amin, amax),
        kwargs...
    )
end

################################################################################
# General edge adding functions
################################################################################
function add_edges_cycle!(ogm::OnesGM{S}, cycle; weight=1.0) where {S}
    for i in 1:numstates(ogm)
        state = itostate(i, ogm)
        matchindices = findall_subvects(S(), state, cycle[1])
        for starti in matchindices
            state1 = copy(state)
            state2 = copy(state)
            for ci in 1:(length(cycle)-1)
                subin!(state1, starti, cycle[ci])
                subin!(state2, starti, cycle[ci+1])
                inc_edge!(ogm.graph, statetoi(state1), statetoi(state2), weight)
            end
            subin!(state1, starti, cycle[end])
            subin!(state2, starti, cycle[begin])
            inc_edge!(ogm.graph, statetoi(state1), statetoi(state2), weight)
        end
    end
end

function add_edges_wcycle!(ogm::OnesGM{S}, wcycle) where {S}
    for i in 1:numstates(ogm)
        state = itostate(i, ogm)
        matchindices = findall_subvects(S(), state, wcycle[1][1])
        for starti in matchindices
            state1 = copy(state)
            state2 = copy(state)
            for ci in 1:(length(wcycle)-1)
                subin!(state1, starti, wcycle[ci][1])
                subin!(state2, starti, wcycle[ci+1][1])
                inc_edge!(ogm.graph, statetoi(state1), statetoi(state2), wcycle[ci][2])
            end
            subin!(state1, starti, wcycle[end][1])
            subin!(state2, starti, wcycle[begin][1])
            inc_edge!(ogm.graph, statetoi(state1), statetoi(state2), wcycle[end][2])
        end
    end
end

"""
This should really accomplish the same as add_edges_cycle! just scaled
by the length of the cycle but its good to have this as a test.
"""
function add_edges_cycleall!(ogm::OnesGM, cycle; kwargs...)
    cycle_len = length(cycle)
    for si in 1:cycle_len
        cycle_ = cycle[mod1.(si:si+cycle_len-1, cycle_len)]
        add_edges_cycle!(ogm, cycle_; kwargs...)
    end
end

################################################################################
# Making OnesGMs
################################################################################
# Basic equal rate FC makers
function make_single_cycle(N, cycle, symmetry::Symmetry=Loop(); kwargs...)
    ogm = OnesGM(N;
        symmetry,
        metadata=Dict("chash" => hash(cycle), "ctype" => "simple")
    )
    add_edges_cycle!(ogm, cycle; kwargs...)
    ogm
end

function make_multi_cycle(N, symmetry::Symmetry, cyclesandweights...; kwargs...)
    ogm = OnesGM(N;
        symmetry,
        metadata=Dict("chash" => hash(cyclesandweights), "ctype" => "simple")
    )
    if isodd(length(cyclesandweights))
        throw(ArgumentError("cyclesandweights is not even"))
    end
    for i in 1:div(length(cyclesandweights), 2)
        cycle = cyclesandweights[2*i-1]
        weight = cyclesandweights[2*i]
        add_edges_cycle!(ogm, cycle; weight)
    end
    ogm
end

# Same but have different rates
function make_single_wcycle(N, wcycle, symmetry::Symmetry=Loop(); kwargs...)
    ogm = OnesGM(N;
        symmetry,
        metadata=Dict("chash" => hash(wcycle), "ctype" => "simple")
    )
    add_edges_wcycle!(ogm, wcycle; kwargs...)
    ogm
end

"""
Get any of the interesting 2 digit cycles
"""
function get_c2d_int(i)
    xx = [
        [[0, 0], [1, 0], [1, 1], [0, 1]],
        [[0, 0], [1, 0], [0, 1], [1, 1]],
        [[0, 0], [1, 0], [0, 1]],
        [[0, 0], [1, 0], [1, 1]]
    ]
    xx[i]
end

################################################################################
# Util
################################################################################
plotgm_save_kwargs() = (;
    n_size=10.0,
    nlabels_fontsize=6,
    flabels=true,
    fnlabels=:repr
)
