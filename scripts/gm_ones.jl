using DrWatson
@quickactivate "TopoStochSim"

using Revise

includet(srcdir("gm_ca.jl"))

################################################################################
# Base setup
################################################################################
struct OnesGM{F} <: AbstractGraphModel
    N::Int
    graph::SimpleWeightedDiGraph{Int,F}
    function OnesGM(N;
        numtype::Val{F}=Val(Float64),
        graph=nothing
    ) where {F}
        if isnothing(graph)
            graph = SimpleWeightedDiGraph{Int,F}(numstates(N))
        elseif nv(graph) != numstates(N)
            throw(ArgumentError("invalid graph passed, does not have the correct number of nodes"))
        end

        new{F}(N, graph)
    end
end
function graph(ogm::OnesGM)
    ogm.graph
end
function copy(ogm::OnesGM{F}) where {F}
    OnesGM(ogm.N; numtype=Val(F), graph=copy(ogm.graph))
end

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
# Particular setups
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

function make_v1(N; kwargs...)
    ogm = OnesGM(N)
    add_edges_cyclestart!(ogm, get_base_cycle1(); kwargs...)
    ogm
end

function make_v2(N; kwargs...)
    ogm = OnesGM(N)
    add_edges_cycleall!(ogm, get_base_cycle1(); kwargs...)
    ogm
end

################################################################################
# Plotting
################################################################################
function p_do_layout(ogm::OnesGM, layout, roffset_devs)
    if isnothing(layout)
        layout = :tree
    end
    if isnothing(roffset_devs)
        roffset_devs = :auto
    end

    # this may or may not be filled
    axis_labels = nothing

    if isa(layout, NetworkLayout.AbstractLayout)
        layout = layout(p_get_adjmat_symsafe(ogm))
        dim = length(layout[1])
    elseif isa(layout, Function)
        layout = layout(ogm)
        dim = length(layout[1])
    elseif isa(layout, AbstractArray) && (length(layout) == numstates(ogm))
        dim = length(layout[1])
    else
        if isa(layout, Symbol)
            layout = (layout,)
        end
        (type, layout_args...) = layout

        # These are mostly the case, override where not
        do_default_roffsets = false
        dim = 3

        if type == :tree
            groups = [[] for _ in 1:ogm.N+1]
            for i in 1:numstates(ogm)
                push!(groups[1+count(x -> x == 1, itostate(i, ogm))], i)
            end
            for g in groups
                sort!(g; by=i -> itostate(i, ogm), rev=true)
            end
            layout = Vector{Any}(undef, numstates(ogm))
            for (gi, group) in enumerate(groups)
                mid = div(length(group), 2, RoundUp)
                for (i, reali) in enumerate(group)
                    layout[reali] = (float(i - mid), float(gi))
                end
            end

            maxs = maximum(layout)
            axis_labels = (f"N({maxs[1]})", "")
            dim = 2
        else
            throw(ArgumentError(f"layout of {layout} is not recognised"))
        end

        if roffset_devs == :auto
            if do_default_roffsets
                roffset_devs = 0.01
            else
                roffset_devs = nothing
            end
        elseif roffset_devs == true
            roffset_devs = 0.01
        elseif roffset_devs == false
            roffset_devs = nothing
        end

        if !isnothing(roffset_devs)
            if length(roffset_devs) == 1
                roffset_devs = Tuple(roffset_devs for _ in 1:length(layout[1]))
            end
            for i in eachindex(layout)
                offsets = Tuple(randn() * dev for dev in roffset_devs)
                layout[i] = layout[i] .+ offsets
            end
        end
    end
    (; layout, dim, axis_labels)
end

function plot_ogm_!(ax, ogm::OnesGM{F},
    do_layout_rslt, args...;
    axparent=nothing,
    # these apply to all
    flabels=:auto,
    fnlabels=nothing,
    felabels=nothing,
    node_size_scale=10.0,
    interactive=:auto,
    symbolic=:auto,
    # these only to concrete (non symbolic) cas
    c_colormap=:dense,
    c_colorscale=:auto,
    c_colorrange=:auto,
    c_colorbar=:auto,
    kwargs...
) where {F}
    layout, dim, axis_labels = do_layout_rslt
    if symbolic == :auto
        symbolic = F == Num
    end

    auto_kwargs = Dict{Symbol,Any}()

    # This is for possible later elementwise updates
    auto_kwargs[:node_color] = [:snow3 for _ in 1:numstates(ogm)]
    auto_kwargs[:node_size] = [node_size_scale for _ in 1:numstates(ogm)]

    # Do colored transitions
    if !symbolic
        weights = weight.(edges(ogm.graph))
        auto_kwargs[:edge_color] = weights
        if c_colorrange == :auto
            max_weight = maximum(weights)
            min_weight = minimum(weights)
            if max_weight == min_weight
                min_weight -= 0.1
            end
            c_colorrange = (min_weight, max_weight)
        end
        if c_colorscale == :auto
            c_colorscale = x -> Makie.pseudolog10(x) * log(10)
        end
        edge_colormap_attrs = (;
            colormap=c_colormap,
            colorrange=c_colorrange,
            colorscale=c_colorscale
        )
        auto_kwargs[:arrow_attr] = auto_kwargs[:edge_attr] = edge_colormap_attrs
    end

    # Label nodes and or edges
    if flabels == :auto
        flabels = nv(ogm.graph) < 100
        if isnothing(fnlabels)
            fnlabels = :repr
        end
        if isnothing(felabels)
            felabels = :repr
        end
    end
    if flabels
        if fnlabels == :repr
            auto_kwargs[:nlabels] = repr.(allstates(ogm))
        end
        if felabels == :repr
            auto_kwargs[:elabels] = repr.(weight.(edges(ogm.graph)))
            auto_kwargs[:elabels_shift] = 0.4
        end
    end

    plot = graphplot!(ax, ogm.graph, args...;
        layout,
        auto_kwargs...,
        kwargs...
    )

    if interactive == :auto
        interactive = (dim == 2)
    end
    if interactive
        make_interactive(ax, plot)
    end

    if !isnothing(axis_labels)
        if dim == 3
            xlabel!(ax.scene, axis_labels[1])
            ylabel!(ax.scene, axis_labels[2])
            zlabel!(ax.scene, axis_labels[3])
        elseif dim == 2
            ax.xlabel[] = axis_labels[1]
            ax.ylabel[] = axis_labels[2]
        end
    end

    if !symbolic && !isnothing(axparent)
        if (c_colorbar == :auto) || c_colorbar
            Colorbar(axparent[1, 2]; colormap=c_colormap, colorrange=c_colorrange, scale=c_colorscale)
        end
    end

    plot
end

function plot_ogm!(maybeax, ogm::OnesGM, args...;
    map=nothing, layout=nothing, roffset_devs=nothing, returnax=false, kwargs...
)
    if !isnothing(map)
        ogm = map(ogm)
    end

    do_layout_rslt = p_do_layout(ogm, layout, roffset_devs)

    if isa(maybeax, Makie.AbstractAxis)
        ax = maybeax
    elseif isa(maybeax, GridPosition)
        ax = p_make_ca_ax(do_layout_rslt.dim, maybeax)
    elseif isa(maybeax, GridLayout)
        ax = p_make_ca_ax(do_layout_rslt.dim, maybeax[1, 1])
    else
        throw(ArgumentError(f"maybeax type {typeof(maybeax)} is not recognised, must be either Makie.AbstractAxis or GridPosition"))
    end

    plot = plot_ogm_!(ax, ogm, do_layout_rslt, args...; kwargs...)

    if returnax
        ax, plot
    else
        plot
    end
end
function plot_ogm(ogm::OnesGM, args...;
    map=nothing, layout=nothing, roffset_devs=nothing, kwargs...
)
    if !isnothing(map)
        ogm = map(ogm)
    end

    do_layout_rslt = p_do_layout(ogm, layout, roffset_devs)

    fig = Figure()
    ax = p_make_ca_ax(do_layout_rslt.dim, fig[1, 1])

    plot = plot_ogm_!(ax, ogm, do_layout_rslt, args...; axparent=fig, kwargs...)

    Makie.FigureAxisPlot(fig, ax, plot)
end

function plotGM(ogm::OnesGM, args...; kwargs...)
    plot_ogm(ogm, args...; kwargs...)
end

function plot_ogm_min(args...; ecutoff=1.1, amin=0.2, amax=1.0, kwargs...)
    plot_ogm(args...;
        layout=makeprefunclayout(Spring(), filter_edges!, ecutoff),
        felabels=false,
        c_colormap=make_linalpha_cmap(:dense; amin, amax),
        kwargs...
    )
end

function makespringplots()
    for n in 3:8
        ogm = make_v1(n)
        fap = plot_ogm_min(ogm; node_size_scale=6., nlabels_fontsize=6, flabels=true, fnlabels=:repr)
        savefig("ones", savename("spring", ogm), fap.figure)
    end
end
