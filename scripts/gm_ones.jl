using DrWatson
@quickactivate "TopoStochSim"

using Revise

includet(srcdir("gm.jl"))

################################################################################
# Base setup
################################################################################
struct OnesGM{F} <: AbstractGraphModel{F}
    N::Int
    graph::SimpleWeightedDiGraph{Int,F}
    metadata::Union{Nothing,Dict}
    function OnesGM(N;
        numtype::Val{F}=Val(Float64),
        graph=nothing,
        metadata=nothing
    ) where {F}
        if isnothing(graph)
            graph = SimpleWeightedDiGraph{Int,F}(numstates(N))
        elseif nv(graph) != numstates(N)
            throw(ArgumentError("invalid graph passed, does not have the correct number of nodes"))
        end

        new{F}(N, graph, metadata)
    end
end
function copy(ogm::OnesGM{F}) where {F}
    OnesGM(ogm.N; numtype=Val(F), graph=copy(ogm.graph), metadata=copy(ogm.metadata))
end
function show(io::IO, mime::MIME"text/plain", ogm::OnesGM{F}) where {F}
    println(f"OnesGM{{{F}}}(N={ogm.N})")
    print(f"graph: ")
    show(io, mime, ogm.graph)
    if !isnothing(ogm.metadata)
        print(f"\nmetadata: ")
        show(io, mime, ogm.metadata)
    end
end
graph(ogm::OnesGM) = ogm.graph

numstates(N) = 2^N
numstates(ogm::OnesGM) = numstates(ogm.N)

itostate(i, N) = digits(i - 1; base=2, pad=N)
itostate(i, ogm::OnesGM) = itostate(i, ogm.N)

statetoi(state) = 1 + sum(state[i] * 2^(i - 1) for i in 1:length(state))

allstates(N) = itostate.(1:numstates(N), N)
allstates(ogm::OnesGM) = allstates(ogm.N)

function findall_subvects(vect::Vector, subvect::Vector)
    vect_len = length(vect)
    subvect_len = length(subvect)
    matches = []
    for i in 1:length(vect)
        if (@view vect[mod1.(i:i+subvect_len-1, vect_len)]) == subvect
            push!(matches, i)
        end
    end
    matches
end

function subin!(vect::Vector, i, subvect::Vector)
    for j in 1:length(subvect)
        vect[mod1(i + j - 1, length(vect))] = subvect[j]
    end
end

################################################################################
# Plotting
################################################################################
function p_named_layouts(ogm::OnesGM, layout_name, layout_args)
    try
        invoke(p_named_layouts, Tuple{supertype(typeof(ogm)),Any,Any}, ogm, layout_name, layout_args)
    catch ArgumentError
        def_roffsets = false
        axis_labels = nothing

        if layout_name == :tree
            dim = 2
            groups = [[] for _ in 1:ogm.N+1]
            for i in 1:numstates(ogm)
                push!(groups[1+count(x -> x == 1, itostate(i, ogm))], i)
            end
            for g in groups
                sort!(g; by=i -> itostate(i, ogm), rev=true)
            end
            layout = Vector{Any}(undef, numstates(ogm))
            for (gi, group) in enumerate(groups)
                mid = length(group) / 2
                for (i, reali) in enumerate(group)
                    layout[reali] = (i - mid, float(gi))
                end
            end

            maxs = maximum(layout)
            axis_labels = (f"N({maxs[1]})", "")
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

plotgm_kwarg_defaults(_::OnesGM) = (; fnlabels=:repr)

function plot_ogm_min(args...; ecutoff=1.1, amin=0.2, amax=1.0, kwargs...)
    plot_ogm(args...;
        layout=makeprefunclayout(Spring(), filter_edges!, ecutoff),
        felabels=false,
        c_colormap=make_linalpha_cmap(:dense; amin, amax),
        kwargs...
    )
end

################################################################################
# Adding edges and running
################################################################################
function add_edges_cyclestart!(ogm::OnesGM, cycle; weight=1.0)
    for i in 1:numstates(ogm)
        state = itostate(i, ogm)
        matchindices = findall_subvects(state, cycle[1])
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

function add_edges_cycleall!(ogm::OnesGM, cycle; kwargs...)
    cycle_len = length(cycle)
    for si in 1:cycle_len
        cycle_ = cycle[mod1.(si:si+cycle_len-1, cycle_len)]
        add_edges_cyclestart!(ogm, cycle_; kwargs...)
    end
end

get_base_cycle1() = [[0, 0], [1, 0], [1, 1], [0, 1]]
get_base_cycle2() = [[0, 0, 0], [0, 1, 0], [1, 1, 1], [1, 0, 1]]
get_base_cycle3() = [[0, 0], [1, 0], [0, 1], [1, 1]]

function get_bc_int(i)
    xx = [
        [[0, 0], [1, 0], [1, 1], [0, 1]],
        [[0, 0], [1, 0], [0, 1], [1, 1]],
        [[0, 0], [1, 0], [0, 1]],
        [[0, 0], [1, 0], [1, 1]]
    ]
    xx[i]
end

function make_v1(N; cycle=get_base_cycle1(), kwargs...)
    ogm = OnesGM(N; metadata=Dict("chash" => hash(cycle), "ctype" => "simple"))
    add_edges_cyclestart!(ogm, cycle; kwargs...)
    ogm
end

function make_v2(N; cycle=get_base_cycle1(), kwargs...)
    ogm = OnesGM(N; metadata=Dict("chash" => hash(cycle), "ctype" => "full"))
    add_edges_cycleall!(ogm, cycle; kwargs...)
    ogm
end

function make_v3(N)
    ogm = OnesGM(N; metadata=Dict("ctype" => "full"))
    add_edges_cycleall!(ogm, [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    add_edges_cycleall!(ogm, [[1, 0, 0], [1, 1, 0], [1, 0, 1]])
    add_edges_cycleall!(ogm, [[0, 1, 0], [1, 1, 0], [0, 1, 1]])
    add_edges_cycleall!(ogm, [[0, 0, 1], [1, 0, 1], [0, 1, 1]])
    add_edges_cycleall!(ogm, [[0, 1, 1], [1, 1, 1]])
    add_edges_cycleall!(ogm, [[1, 0, 1], [1, 1, 1]])
    add_edges_cycleall!(ogm, [[1, 1, 0], [1, 1, 1]])
    ogm
end

function makespringplots(dirname, ns=3, ne=6; kwargs...)
    for n in ns:ne
        ogm = make_v1(n; kwargs...)
        fap = plot_ogm_min(ogm; node_size_scale=6.0, nlabels_fontsize=6, flabels=true, fnlabels=:repr)
        savefig(f"ones/{dirname}/", "spring", ogm, fap.figure)
    end
end

function makespringplots1()
    makespringplots("c1", 3, 8)
end

function makespringplots2()
    for i in 1:4
        for n in 3:6
            ogm = make_v1(n; cycle=get_bc_int(i))
            fap = plot_ogm_min(ogm; node_size_scale=6.0, nlabels_fontsize=6, flabels=true, fnlabels=:repr)
            savefig(f"ones/c{i}/", "spring", ogm, fap.figure)
        end
    end
end

function makespringplots3()
    for i in 1:4
        for n in 3:6
            ogm = make_v1(n; cycle=get_bc_int(i))
            fap = plot_ogm_min(ogm; node_size_scale=6.0, nlabels_fontsize=6, flabels=true, fnlabels=:repr, layout=:tree)
            savefig(f"ones/tc{i}/", "spring", ogm, fap.figure)
        end
    end
end
