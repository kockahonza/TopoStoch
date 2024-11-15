################################################################################
# Main plotting
################################################################################
struct FigureAxisAnything
    figure::Figure
    axis::Any
    obj::Any
end
display(faa::FigureAxisAnything) = display(faa.figure)

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

# axis setup
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

# layout setup
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

# fancy plotting
function plotgm_!(ax, gm::AbstractGraphModel{F},
    (; dim, layout), args...;
    axparent=nothing,       # needed for adding colorbars
    # these apply to all
    interactive=:auto,
    flabels=:auto,          # labeling nodes and edges
    fnlabels=:index,
    felabels=:repr,
    # the following can only be used for non-symbolics GMs
    ss=:auto,            # calculate steady state?
    # dealing with nodes
    n_size=30.0,
    n_color=:snow3,
    n_ss_size=:sqrt,
    n_ss_color=true,
    n_ss_colormap=:viridis,
    n_ss_colorscale=identity,
    n_ss_colorrange=:auto,
    n_ss_colorbar=:auto,
    # color mappings for edge weights
    e_color=:auto,
    e_colormap=:dense,
    e_colorscale=:auto,
    e_colorrange=:auto,
    e_colorbar=:auto,
    kwargs...
) where {F}
    auto_kwargs = Dict{Symbol,Any}()

    # Dealing with nodes and steadystate calcs
    if ss == :auto
        if F <: AbstractFloat
            ss = supersteadystate(gm)
        else
            ss = nothing
        end
    end
    if !isnothing(ss)
        # do node sizes which can only be based on ss
        if !(isnothing(n_ss_size) || (n_ss_size == false))
            if n_ss_size == :lin
                auto_kwargs[:node_size] = [n_size * x for x in ss]
            elseif n_ss_size == :log
                auto_kwargs[:node_size] = [n_size * log(1 + x) for x in ss]
            elseif n_ss_size == :sqrt
                auto_kwargs[:node_size] = [n_size * sqrt(x) for x in ss]
            elseif n_ss_size == :linA
                auto_kwargs[:node_size] = [n_size * (x)^(1 / dim) for x in ss]
            elseif n_ss_size == :logA
                auto_kwargs[:node_size] = [n_size * (log(1 + x))^(1 / dim) for x in ss]
            elseif n_ss_size == :sqrtA
                auto_kwargs[:node_size] = [n_size * sqrt(x)^(1 / dim) for x in ss]
            else
                throw(ArgumentError(f"n_ss_size of {n_ss_size} is not recognised"))
            end
        end
        # do node color mapping which can also only be based on ss
        if n_ss_color
            auto_kwargs[:node_color] = ss
            if n_ss_colorrange == :auto
                max_weight = maximum(ss)
                min_weight = minimum(ss)
                if max_weight == min_weight
                    min_weight -= 0.1
                end
                n_ss_colorrange = (min_weight, max_weight)
            end
            auto_kwargs[:node_attr] = (;
                colormap=n_ss_colormap,
                colorrange=n_ss_colorrange,
                colorscale=n_ss_colorscale
            )
        end
    end
    # Make sure code sizes and colors are set to some array for later modifications
    if !haskey(auto_kwargs, :node_size)
        auto_kwargs[:node_size] = [n_size for _ in 1:numstates(gm)]
    end
    if !haskey(auto_kwargs, :node_color)
        auto_kwargs[:node_color] = [n_color for _ in 1:numstates(gm)]
    end

    # Color edges
    if e_color == :auto
        if F <: AbstractFloat
            e_color = :weight
        else
            e_color = nothing
        end
    elseif e_color == false
        e_color = nothing
    end
    if !isnothing(e_color)
        auto_kwargs[:edge_color] = colors = if e_color == :weight
            if e_colorscale == :auto
                e_colorscale = :pseudolog
            end
            weight.(edges(graph(gm)))
        elseif e_color == :current
            if isnothing(ss)
                throw(ArgumentError(f"cannot use e_color of {e_color} without ss"))
            end
            if e_colorscale == :auto
                e_colorscale = :lin
            end
            [e.weight * ss[e.src] for e in edges(graph(gm))]
        else
            throw(ArgumentError(f"e_color of {e_color} is not recognised"))
        end

        if e_colorrange == :auto
            max_weight = maximum(colors)
            min_weight = minimum(colors)
            if max_weight == min_weight
                min_weight -= 0.1
            end
            e_colorrange = (min_weight, max_weight)
        end
        if e_colorscale == :pseudolog
            e_colorscale = x -> Makie.pseudolog10(x) * log(10)
        elseif e_colorscale == :lin
            e_colorscale = identity
        else
            throw(ArgumentError(f"e_colorscale of {e_colorscale} is not recognised"))
        end
        edge_colormap_attrs = (;
            colormap=e_colormap,
            colorrange=e_colorrange,
            colorscale=e_colorscale
        )
        auto_kwargs[:arrow_attr] = auto_kwargs[:edge_attr] = edge_colormap_attrs
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
        if !isnothing(e_color) && ((e_colorbar == :auto) || e_colorbar)
            Colorbar(axparent[1, 2]; colormap=e_colormap, colorrange=e_colorrange, scale=e_colorscale)
        end
        if !isnothing(ss) & n_ss_color && ((n_ss_colorbar == :auto) || n_ss_colorbar)
            Colorbar(axparent[1, 3]; colormap=n_ss_colormap, colorrange=n_ss_colorrange, scale=n_ss_colorscale)
        end
    end

    plot
end

plotgm_kwarg_defaults(_::AbstractGraphModel) = (;)

function plotgm!(maybeax, gm::AbstractGraphModel, args...;
    layout=nothing, roffset_devs=nothing, returnax=false, kwargs...
)
    do_layout_rslt = p_do_layout(gm, layout, roffset_devs)

    if isa(maybeax, Makie.AbstractAxis)
        ax = maybeax
    elseif isa(maybeax, GridPosition)
        ax = p_make_ax(do_layout_rslt.dim, maybeax, do_layout_rslt.axis_labels)
    elseif isa(maybeax, GridLayout)
        ax = p_make_ax(do_layout_rslt.dim, maybeax[1, 1], do_layout_rslt.axis_labels)
    else
        throw(ArgumentError(f"maybeax type {typeof(maybeax)} is not recognised, must be Makie.AbstractAxis, GridPosition or GridLayout"))
    end

    plot = plotgm_!(ax, gm, do_layout_rslt, args...; plotgm_kwarg_defaults(gm)..., kwargs...)

    if returnax
        ax, plot
    else
        plot
    end
end
function plotgm(gm::AbstractGraphModel, args...;
    layout=nothing, roffset_devs=nothing, kwargs...
)
    do_layout_rslt = p_do_layout(gm, layout, roffset_devs)

    fig = Figure()
    ax = p_make_ax(do_layout_rslt.dim, fig[1, 1], do_layout_rslt.axis_labels)

    plot = plotgm_!(ax, gm, do_layout_rslt, args...; axparent=fig, plotgm_kwarg_defaults(gm)..., kwargs...)

    Makie.FigureAxisPlot(fig, ax, plot)
end

################################################################################
# Plotting util
################################################################################
function make_linalpha_cmap(cmap; amin=0.0, amax=1.0)
    cmap = Makie.to_colormap(cmap)
    acmap = [coloralpha(color(c), a) for (c, a) in zip(cmap, LinRange(amin, amax, length(cmap)))]
    cgrad(acmap)
end

function makeprefunclayout(base_layout::NetworkLayout.AbstractLayout, f, args...; kwargs...)
    copyf = copyand(f, args...; kwargs...)
    function (ca)
        base_layout(p_safe_adjmat(copyf(ca)))
    end
end

# adding extras
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
# Saving
################################################################################
function savefig(subdir, prefix, gm::AbstractGraphModel, fig)
    savefig(subdir,
        savename(prefix, gm;
            allowedtypes=[DrWatson.default_allowed(gm)...; Dict],
            expand=["metadata"]),
        fig)
end

function save_transmat(gm::AbstractGraphModel; name=savename(f"tm_{string(typeof(gm).name.name)}", gm, "table"), short=false)
    file = open(datadir(name), "w")

    row_names = repr.(1:numstates(gm)) .* " / " .* repr.(allstates(gm))
    tm = transmat(gm; mat=true)
    if short
        header = [""; repr.(1:numstates(gm))]
        tm = map(x -> if iszero(x)
                0
            else
                1
            end, tm)
    else
        header = [""; row_names]
    end
    modified_matrix = [row_names tm]

    pretty_table(file, modified_matrix; header)
    close(file)
end

function save_ca_obj(gm::AbstractGraphModel; kwargs...)
    for (name, value) in kwargs
        save_object(datadir(savename(String(name), gm, "jld2")), value)
    end
end

