################################################################################
# Plotting util
################################################################################
struct FigureAxisAnything
    figure::Figure
    axis::Any
    obj::Any
end
display(faa::FigureAxisAnything) = display(faa.figure)

function make_linalpha_cmap(cmap; amin=0.0, amax=1.0)
    cmap = Makie.to_colormap(cmap)
    acmap = [coloralpha(color(c), a) for (c, a) in zip(cmap, LinRange(amin, amax, length(cmap)))]
    cgrad(acmap)
end

function makeprefunclayout(base_layout::NetworkLayout.AbstractLayout, f, args...; kwargs...)
    copyf = copyand(f, args...; kwargs...)
    function (ca)
        base_layout(p_get_adjmat_symsafe(copyf(ca)))
    end
end

# adding extras
# FIX: This is somewhat redundant, uses should be removed
function edgecolorbar(fig, plot)
    edgeploti = findfirst(typeof(plot_) <: Plot{GraphMakie.edgeplot} for plot_ in plot.plots)
    Colorbar(fig[1, 2], plot.plots[edgeploti])
end
edgecolorbar((fig, ax, plot)) = edgecolorbar(fig, plot)

function make_interactive(ax, plot)
    interactions_ = keys(interactions(ax))
    if :rectanglezoom in interactions_
        deregister_interaction!(ax, :rectanglezoom)
    end
    if :ndrag in interactions_
        deregister_interaction!(ax, :ndrag)
    end
    function node_drag_action(state, idx, event, axis)
        plot[:node_pos][][idx] = event.data
        plot[:node_pos][] = plot[:node_pos][]
    end
    ndrag = NodeDragHandler(node_drag_action)
    register_interaction!(ax, :ndrag, ndrag)
end
make_interactive((_, ax, plot)) = make_interactive(ax, plot)

################################################################################
# TODO: Remove: Old Main plotting
################################################################################
function p_get_adjmat_symsafe(graph::AbstractNumGraph)
    map(x -> if iszero(x)
            0.0
        else
            1.0
        end, adjacency_matrix(graph))
end
function p_get_adjmat_symsafe(graph::AbstractFloatGraph)
    adjacency_matrix(graph)
end
p_get_adjmat_symsafe(gm::AbstractGraphModel) = p_get_adjmat_symsafe(graph(gm))

################################################################################
# Main plotting
################################################################################
function p_safe_adjmat(graph::AbstractNumGraph)
    map(x -> if iszero(x)
            0.0
        else
            1.0
        end, adjacency_matrix(graph))
end
function p_safe_adjmat(graph::AbstractFloatGraph)
    adjacency_matrix(graph)
end
p_safe_adjmat(gm::AbstractGraphModel) = p_safe_adjmat(graph(gm))

function p_set_axis_labels!(ax, dim, axis_labels)
    if !isnothing(axis_labels)
        if dim == 2
            ax.xlabel[] = axis_labels[1]
            ax.ylabel[] = axis_labels[2]
        elseif dim == 3
            xlabel!(ax.scene, axis_labels[1])
            ylabel!(ax.scene, axis_labels[2])
            zlabel!(ax.scene, axis_labels[3])
        else
            throw(ArgumentError("the layout dimension was not 2 or 3"))
        end
    end
end
function p_make_ax(dim, place, axis_labels=nothing)
    ax = if dim == 2
        Axis(place)
    elseif dim == 3
        LScene(place)
    else
        throw(ArgumentError("the layout dimension was not 2 or 3"))
    end

    if !isnothing(axis_labels)
        p_set_axis_labels!(ax, dim, axis_labels)
    end

    ax
end

function p_named_layouts(gm::AbstractGraphModel, layout_name, layout_args)
    def_roffsets = false
    axis_labels = nothing

    if layout_name == :random2
        dim = 2
        layout = [rand(2) for _ in 1:nv(graph(gm))]
    else
        throw(ArgumentError(f"layout of name {layout_name} is not recognised"))
    end

    dim, layout, def_roffsets, axis_labels
end
function p_do_layout(gm::AbstractGraphModel, layout=nothing, roffset_devs=nothing)
    if isnothing(layout)
        layout = Spring()
    end
    if isnothing(roffset_devs)
        roffset_devs = :auto
    end

    def_roffsets = false
    axis_labels = nothing

    if isa(layout, NetworkLayout.AbstractLayout)
        layout = layout(p_safe_adjmat(gm))
        dim = length(layout[1])
    elseif isa(layout, Function)
        layout = layout(gm)
        dim = length(layout[1])
    elseif isa(layout, AbstractArray) && (length(layout) == nv(graph(gm)))
        dim = length(layout[1])
    else
        if isa(layout, Symbol)
            layout = (layout,)
        end
        (layout_name, layout_args...) = layout
        dim, layout, def_roffsets, axis_labels = p_named_layouts(gm, layout_name, layout_args)
    end

    if roffset_devs == :auto
        if def_roffsets
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

    (; dim, layout, axis_labels)
end

function p_ss_transform(transformtype, ss, dim)
    if transformtype == :lin
        ss
    elseif transformtype == :log
        [log(1 + x) for x in ss]
    elseif transformtype == :sqrt
        [sqrt(x) for x in ss]
    elseif transformtype == :linA
        [(x)^(1 / dim) for x in ss]
    elseif transformtype == :logA
        [(log(1 + x))^(1 / dim) for x in ss]
    elseif transformtype == :sqrtA
        [sqrt(x)^(1 / dim) for x in ss]
    else
        throw(ArgumentError(f"transformtype of {transformtype} is not recognised"))
    end
end

function plotgm_!(ax, gm::AbstractGraphModel{F},
    (; dim, layout), args...;
    axparent=nothing,       # needed for adding colorbars
    # these apply to all
    interactive=:auto,
    flabels=:auto,          # labeling nodes and edges
    fnlabels=:index,
    felabels=:repr,
    # the following can only be used for non-symbolics GMs
    # color mappings for edge weights
    e_color=:auto,
    e_colormap=:dense,
    e_colorscale=:pseudolog,
    e_colorrange=:auto,
    e_colorbar=:auto,
    # dealing with nodes and steadystate calcs
    n_size=30.0,
    n_color=:snow3,
    do_ss=:auto,
    ss_size=:sqrt,
    ss_color=true,
    ss_colormap=:viridis,
    ss_colorscale=identity,
    ss_colorrange=:auto,
    ss_colorbar=:auto,
    kwargs...
) where {F}
    auto_kwargs = Dict{Symbol,Any}()

    # Color edges by weight
    if e_color == :auto
        e_color = F <: AbstractFloat
    end
    if e_color
        weights = weight.(edges(gm.graph))
        auto_kwargs[:edge_color] = weights
        if e_colorrange == :auto
            max_weight = maximum(weights)
            min_weight = minimum(weights)
            if max_weight == min_weight
                min_weight -= 0.1
            end
            e_colorrange = (min_weight, max_weight)
        end
        if e_colorscale == :pseudolog
            e_colorscale = x -> Makie.pseudolog10(x) * log(10)
        end
        edge_colormap_attrs = (;
            colormap=e_colormap,
            colorrange=e_colorrange,
            colorscale=e_colorscale
        )
        auto_kwargs[:arrow_attr] = auto_kwargs[:edge_attr] = edge_colormap_attrs
    end

    # Dealing with nodes and steadystate calcs
    if do_ss == :auto
        do_ss = F <: AbstractFloat
    end
    if do_ss
        ss = supersteadystate(gm)
        if !(isnothing(ss_size) || (ss_size == false))
            auto_kwargs[:node_size] = n_size .* p_ss_transform(ss_size, ss, dim)
        end
        if ss_color
            auto_kwargs[:node_color] = ss
            if ss_colorrange == :auto
                max_weight = maximum(ss)
                min_weight = minimum(ss)
                if max_weight == min_weight
                    min_weight -= 0.1
                end
                ss_colorrange = (min_weight, max_weight)
            end
            auto_kwargs[:node_attr] = (;
                colormap=ss_colormap,
                colorrange=ss_colorrange,
                colorscale=ss_colorscale
            )
        end
    end
    if !haskey(auto_kwargs, :node_size)
        auto_kwargs[:node_size] = [n_size for _ in 1:numstates(gm)]
    end
    if !haskey(auto_kwargs, :node_color)
        auto_kwargs[:node_color] = [n_color for _ in 1:numstates(gm)]
    end

    # Label nodes and or edges
    if flabels == :auto
        flabels = numstates(gm) < 100
    end
    if flabels
        if fnlabels == :index
            auto_kwargs[:nlabels] = repr.(1:numstates(gm))
        elseif fnlabels == :repr
            auto_kwargs[:nlabels] = repr.(allstates(gm))
        elseif !(isnothing(fnlabels) || (fnlabels == false))
            throw(ArgumentError(f"fnlabels of {fnlabels} is not recognised"))
        end
        if felabels == :repr
            auto_kwargs[:elabels] = repr.(weight.(edges(graph(gm))))
            auto_kwargs[:elabels_shift] = 0.4
        elseif !(isnothing(felabels) || (felabels == false))
            throw(ArgumentError(f"felabels of {felabels} is not recognised"))
        end
    end

    plot = graphplot!(ax, graph(gm), args...;
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

    if !isnothing(axparent)
        if e_color && ((e_colorbar == :auto) || e_colorbar)
            Colorbar(axparent[1, 2]; colormap=e_colormap, colorrange=e_colorrange, scale=e_colorscale)
        end
        if ss_color && ((ss_colorbar == :auto) || ss_colorbar)
            Colorbar(axparent[1, 3]; colormap=ss_colormap, colorrange=ss_colorrange, scale=ss_colorscale)
        end
    end

    plot
end

p_kwarg_defaults(_::AbstractGraphModel) = (;)
function plotgm(gm::AbstractGraphModel, args...;
    layout=nothing, roffset_devs=nothing, kwargs...
)
    do_layout_rslt = p_do_layout(gm, layout, roffset_devs)

    fig = Figure()
    ax = p_make_ax(do_layout_rslt.dim, fig[1, 1], do_layout_rslt.axis_labels)

    plot = plotgm_!(ax, gm, do_layout_rslt, args...; axparent=fig, p_kwarg_defaults(gm)..., kwargs...)

    Makie.FigureAxisPlot(fig, ax, plot)
end


################################################################################
# Util
################################################################################
function savefig(subdir, prefix, gm::AbstractGraphModel, fig)
    savefig(subdir,
        savename(prefix, gm;
            allowedtypes=[DrWatson.default_allowed(gm)...; Dict],
            expand=["metadata"]),
        fig)
end
