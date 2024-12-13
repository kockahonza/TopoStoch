################################################################################
# Plotting basics
################################################################################
struct FigureAxisAnything
    figure::Figure
    axis::Any
    obj::Any
end
display(faa::FigureAxisAnything) = display(faa.figure)
export FigureAxisAnything

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
function p_set_axis_labels!(ax::Makie.AbstractAxis, axis_labels=[])
    if isa(ax, Axis)
        ax.xlabel[] = axis_labels[1]
        ax.ylabel[] = axis_labels[2]
    elseif isa(ax, LScene)
        xlabel!(ax.scene, axis_labels[1])
        ylabel!(ax.scene, axis_labels[2])
        zlabel!(ax.scene, axis_labels[3])
    elseif isa(ax, Axis3)
        ax.xlabel[] = axis_labels[1]
        ax.ylabel[] = axis_labels[2]
        ax.zlabel[] = axis_labels[3]
    else
        throw(ArgumentError(f"unrecognised ax type {typeof(ax)}"))
    end
end

function p_make_ax(dim, maybeax, axis_labels=nothing; useaxis3=false)
    if isa(maybeax, Makie.AbstractAxis)
        # If it is already an ax then just do checks and return it
        if isa(maybeax, Axis) && (dim != 2)
            @warn "in p_make_ax got a prepared Axis but dim is not 2"
        end
        return maybeax
    elseif isa(maybeax, GridPosition)
        place = maybeax
    elseif isa(maybeax, GridLayout) || isa(maybeax, Figure)
        place = maybeax[1, 1]
    else
        throw(ArgumentError(f"maybeax type {typeof(maybeax)} is not recognised, must be Makie.AbstractAxis, GridPosition or GridLayout"))
    end

    # Creates the ax if not already existing
    ax = if dim == 2
        Axis(place)
    elseif dim == 3
        if !useaxis3
            LScene(place)
        else
            Axis3(place)
        end
    else
        throw(ArgumentError("the layout dimension was not 2 or 3"))
    end

    if !isnothing(axis_labels)
        p_set_axis_labels!(ax, axis_labels)
    end

    ax
end
export p_make_ax

function p_make_interactive(ax, plot)
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
p_make_interactive((_, ax, plot)) = p_make_interactive(ax, plot)

################################################################################
# Main plotgm stuff
################################################################################
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
export p_named_layouts, p_do_layout

# fancy plotting
function plotgm_!(ax, gm::AbstractGraphModel{F}, (; dim, layout);
    axparent=nothing,       # needed for adding colorbars
    # these apply to all
    interactive=:auto,
    flabels=:auto,          # labeling nodes and edges
    fnlabels=:index,
    felabels=:auto,
    makwargs=x -> (), # mutating function to apply to auto_kwargs just before plotting
    # the following can only be used for non-symbolics GMs
    ss=:auto,            # calculate steady state?
    ss_curgraph=nothing,
    # dealing with nodes
    n_size=15.0,
    n_color=:snow3,
    n_ss_size=:sqrt,
    n_ss_color=true,
    n_ss_colormap=:Oranges,
    n_ss_colorscale=identity,
    n_ss_colorrange=:auto,
    n_ss_colorbar=:auto,
    n_ss_colorbar_label=:auto,
    # color mappings for edge weights
    e_color=:auto,
    e_colormap=:BuPu,
    e_colorscale=:auto,
    e_colorrange=:auto,
    e_colorbar=:auto,
    e_colorbar_label=:auto,
    colorrange_threshold=1e-8,
    # other possible kwargs to be passed down
    zerothreshold=nothing,
    kwargs...
) where {F}
    auto_kwargs = Dict{Symbol,Any}()

    # Dealing with nodes and steadystate calcs
    if ss == :auto
        if F <: AbstractFloat
            @info "automatically calculating steady state as graph is concrete"
            ss = supersteadystate(gm)
        else
            ss = nothing
        end
    end
    if !isnothing(ss)
        # do node sizes which can only be based on ss
        if !(isnothing(n_ss_size) || (n_ss_size == false))
            n_size *= 4.0
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
                if (max_weight - min_weight) < colorrange_threshold
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

    # The graph that will eventually be plotted, this can be modified!
    plotgraph = graph(gm)

    # Color edges
    if e_color == :auto
        if F <: AbstractFloat
            e_color = :rates
        else
            e_color = nothing
        end
    elseif e_color == false
        e_color = nothing
    end
    if !isnothing(e_color)
        auto_kwargs[:edge_color] = colors = if e_color == :rates
            if e_colorscale == :auto
                e_colorscale = :pseudolog
            end
            if e_colorbar_label == :auto
                e_colorbar_label = "Transition rates"
            end
            weight.(edges(plotgraph))
        elseif e_color == :currents
            if isnothing(ss)
                throw(ArgumentError(f"cannot use e_color of {e_color} without ss"))
            end
            if e_colorbar_label == :auto
                e_colorbar_label = "Probability currents"
            end
            [e.weight * ss[e.src] for e in edges(plotgraph)]
        elseif e_color == :dcurrents
            if isnothing(ss)
                throw(ArgumentError(f"cannot use e_color of {e_color} without ss"))
            end
            if e_colorbar_label == :auto
                e_colorbar_label = "Net probability currents"
            end
            plotgraph = if isnothing(ss_curgraph)
                if isnothing(zerothreshold)
                    make_current_graph(gm, ss)
                else
                    make_current_graph(gm, ss; zerothreshold)
                end
            else
                ss_curgraph
            end
            weight.(edges(plotgraph))
        else
            throw(ArgumentError(f"e_color of {e_color} is not recognised"))
        end

        if e_colorrange == :auto
            max_weight = maximum(colors)
            min_weight = minimum(colors)
            if (max_weight - min_weight) < colorrange_threshold
                min_weight -= 0.1
            end
            e_colorrange = (min_weight, max_weight)
        end
        if e_colorscale == :auto
            e_colorscale = identity
        elseif e_colorscale == :pseudolog
            e_colorscale = x -> Makie.pseudolog10(x) * log(10)
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
        elseif isa(fnlabels, Function)
            auto_kwargs[:nlabels] = fnlabels.(allstates(gm))
        elseif !(isnothing(fnlabels) || (fnlabels == false))
            throw(ArgumentError(f"fnlabels of {fnlabels} is not recognised"))
        end
        if felabels == :auto
            if !(F <: AbstractFloat)
                felabels = :repr
            elseif !isnothing(e_color)
                felabels = :e_color
            end
        end
        if felabels == :e_color
            auto_kwargs[:elabels] = roundrepr.(auto_kwargs[:edge_color])
            auto_kwargs[:elabels_shift] = 0.4
        elseif felabels == :repr
            auto_kwargs[:elabels] = repr.(weight.(edges(plotgraph)))
            auto_kwargs[:elabels_shift] = 0.4
        elseif !(isnothing(felabels) || (felabels == false))
            throw(ArgumentError(f"felabels of {felabels} is not recognised"))
        end
    end

    makwargs(auto_kwargs)

    plot = graphplot!(ax, plotgraph;
        layout,
        auto_kwargs...,
        kwargs...
    )

    if interactive == :auto
        interactive = (dim == 2)
    end
    if interactive
        p_make_interactive(ax, plot)
    end

    if !isnothing(axparent)
        if !isnothing(e_color) && ((e_colorbar == :auto) || e_colorbar)
            if e_colorbar_label == :auto
                e_colorbar_label = string(e_color)
            end
            Colorbar(axparent[1, 2];
                colormap=e_colormap,
                colorrange=e_colorrange,
                scale=e_colorscale,
                label=e_colorbar_label
            )
        end
        if !isnothing(ss) && n_ss_color && ((n_ss_colorbar == :auto) || n_ss_colorbar)
            if n_ss_colorbar_label == :auto
                n_ss_colorbar_label = "Occupation probability"
            end
            Colorbar(axparent[1, 3];
                colormap=n_ss_colormap,
                colorrange=n_ss_colorrange,
                scale=n_ss_colorscale,
                label=n_ss_colorbar_label
            )
        end
    end

    plot
end

plotgm_kwarg_defaults(_::AbstractGraphModel) = (;)
export plotgm_kwarg_defaults

"""
Fancy function to plot a graph model with as many convenience features as possible
"""
function plotgm!(maybeax, gm::AbstractGraphModel;
    layout=nothing, roffset_devs=nothing, returnax=false, kwargs...
)
    do_layout_rslt = p_do_layout(gm, layout, roffset_devs)
    ax = p_make_ax(do_layout_rslt.dim, maybeax, do_layout_rslt.axis_labels)

    plot = plotgm_!(ax, gm, do_layout_rslt;
        plotgm_kwarg_defaults(gm)...,
        kwargs...
    )

    if returnax
        ax, plot
    else
        plot
    end
end
@doc (@doc plotgm!)
function plotgm(gm::AbstractGraphModel;
    layout=nothing, roffset_devs=nothing, kwargs...
)
    fig = Figure()
    do_layout_rslt = p_do_layout(gm, layout, roffset_devs)
    ax = p_make_ax(do_layout_rslt.dim, fig[1, 1], do_layout_rslt.axis_labels)

    plot = plotgm_!(ax, gm, do_layout_rslt;
        axparent=fig,
        plotgm_kwarg_defaults(gm)...,
        kwargs...
    )

    Makie.FigureAxisPlot(fig, ax, plot)
end
export plotgm!, plotgm

################################################################################
# plotcgm - plotting the net currents of some state
################################################################################
function p_do_clayout(
    gm::AbstractGraphModel,
    cgraph::SimpleWeightedDiGraph{<:Integer,<:AbstractFloat},
    layout=nothing, roffset_devs=nothing
)
    if isnothing(layout)
        layout = Spectral(dim=2)
    end

    if isa(layout, NetworkLayout.AbstractLayout)
        layout = layout(cgraph)
    elseif isa(layout, Function)
        layout = layout(cgraph)
    end

    p_do_layout(gm, layout, roffset_devs)
end
export p_do_clayout

"""
Plots the net probability currents at some particular state of a graph model.
Uses plotgm!_ and tries to keep as much functionality as possible.
"""
function plotcgm!(maybeax, gm::AbstractGraphModel, state=supersteadystate(gm);
    layout=nothing, roffset_devs=nothing, returnax=false, kwargs...
)
    curgraph = make_current_graph(gm, state)

    do_layout_rslt = p_do_clayout(gm, curgraph, layout, roffset_devs)
    ax = p_make_ax(do_layout_rslt.dim, maybeax, do_layout_rslt.axis_labels)

    plot = plotgm_!(ax, gm, do_layout_rslt;
        ss=state,
        ss_curgraph=curgraph,
        plotgm_kwarg_defaults(gm)...,
        kwargs...
    )

    if returnax
        ax, plot
    else
        plot
    end
end
@doc (@doc plotcgm!)
function plotcgm(gm::AbstractGraphModel{<:AbstractFloat}, state=supersteadystate(gm);
    layout=nothing, roffset_devs=nothing, kwargs...
)
    curgraph = make_current_graph(gm, state)

    fig = Figure()
    do_layout_rslt = p_do_clayout(gm, curgraph, layout, roffset_devs)
    ax = p_make_ax(do_layout_rslt.dim, fig[1, 1], do_layout_rslt.axis_labels)

    plot = plotgm_!(ax, gm, do_layout_rslt;
        axparent=fig,
        ss=state,
        ss_curgraph=curgraph,
        plotgm_kwarg_defaults(gm)...,
        kwargs...
    )

    Makie.FigureAxisPlot(fig, ax, plot)
end
export plotcgm!, plotcgm

"""
Plot pertrubation state, uses plotcgm
"""
function plot_ps(gm::AbstractGraphModel{<:AbstractFloat}, pstate; kwargs...)
    if maximum(imag, pstate) > 1e-8
        @warn "plot_ps is intended for fully real states but got non-zero imaginary components"
    end
    pstate = real(pstate)
    rss = pstate / maximum(abs, pstate)
    plotcgm(gm, pstate; n_ss_colormap=:redsblues, n_ss_colorscale=identity, n_ss_colorrange=(-1.0, 1.0), kwargs...)
end
export plot_ps

"""
Interactive visualization of a simulation, does a plotgm plot, highlights the
current node and adds controls.
"""
function plotgm_sim(gm::AbstractGraphModel, times, states; kwargs...)
    if length(times) != length(states)
        throw(ArgumentError("`times` and `states` must have the same length"))
    end
    fig = Figure()
    controls = GridLayout(fig[1, 1], 2, 5)
    plotplace = GridLayout(fig[2, 1])

    # setup the progress bar
    sg = SliderGrid(controls[2, :],
        (label="index", range=1:length(times), startvalue=1)
    )
    slider_i = sg.sliders[1]
    obs_i = slider_i.value

    # main controls
    but_prev = Button(controls[1, 1]; label="prev")
    but_next = Button(controls[1, 3]; label="next")
    on(but_prev.clicks) do _
        if obs_i[] > 1
            set_close_to!(slider_i, obs_i[] - 1)
        end
    end
    on(but_next.clicks) do _
        if obs_i[] < length(times)
            set_close_to!(slider_i, obs_i[] + 1)
        end
    end

    running = false
    running_task_channel = nothing
    but_playpause = Button(controls[1, 2]; label="play")
    speed_bar = SliderGrid(controls[1, 4],
        (label="speed", range=0.0:0.01:10.0, startvalue=1.0)
    )
    on(but_playpause.clicks) do _
        running = !running
        if !isnothing(running_task_channel)
            if running
                @warn "killing old runner despite the fact that it should not be running anymore!"
            end
            @info "removing runner task"
            put!(running_task_channel, true)
            running_task_channel = nothing
        end
        if running
            @info "adding runner task"
            running_task_channel = Channel(c -> begin
                exit = false
                while !exit
                    if isready(c) && (take!(c) == true)
                        exit = true
                    end
                    if obs_i[] >= length(times)
                        @info "killing runner as reached end of time series"
                        exit = true
                    end
                    delta_time = times[obs_i[]+1] - times[obs_i[]]
                    sleep(speed_bar.sliders[1].value[] * delta_time)
                    set_close_to!(slider_i, obs_i[] + 1)
                end
            end)
            but_playpause.label = "pause"
        else
            but_playpause.label = "play"
        end
    end

    # shows the time
    Label(controls[1, end], lift(i -> f"t={times[i]:.2f}", obs_i))

    # Make plot itself
    ax, plot = plotgm!(plotplace, gm; returnax=true, axparent=plotplace, kwargs...)

    # highlight
    cur_state_pos = lift(obs_i, plot.node_pos) do i, npos
        Point2f(npos[states[i]])
    end
    scatter!(ax, cur_state_pos; marker='○', markersize=40, color=:green)

    FigureAxisAnything(fig, ax, [plot, but_playpause])
end
function plotgm_sim(gm::AbstractGraphModel, h5data::HDF5.H5DataStore; kwargs...)
    plotgm_sim(gm, h5data["time"], h5data["path"]; kwargs...)
end
export plotgm_sim

################################################################################
# Other plotting util
################################################################################
function make_linalpha_cmap(cmap; amin=0.0, amax=1.0)
    cmap = Makie.to_colormap(cmap)
    acmap = [coloralpha(color(c), a) for (c, a) in zip(cmap, LinRange(amin, amax, length(cmap)))]
    cgrad(acmap)
end
export make_linalpha_cmap

function makeprefunclayout(base_layout::NetworkLayout.AbstractLayout, f, args...; kwargs...)
    copyf = copyand(f, args...; kwargs...)
    function (ca)
        base_layout(p_safe_adjmat(copyf(ca)))
    end
end
export makeprefunclayout

# More interactivity
function prep_sliders(variables; ranges=nothing, startvalues=nothing)
    slider_specs = []
    for i in eachindex(variables)
        range_ = isnothing(ranges) ? (0.0:0.1:10.0) : ranges[i]
        startvalue_ = isnothing(startvalues) ? 1.0 : startvalues[i]
        push!(slider_specs, (;
            label=string(variables[i]),
            range=range_,
            startvalue=startvalue_
        ))
    end

    fig = Figure()
    sg = SliderGrid(fig[1, 1], slider_specs...)

    fig, [x.value for x in sg.sliders], sg
end
export prep_sliders

################################################################################
# Saving
################################################################################
function savefig(subdir, basename, fig; kwargs...)
    mkpath(plotsdir(subdir))
    save(plotsdir(subdir, basename * ".png"), fig;
        px_per_unit=2.0,
        kwargs...
    )
end
savefig(basename, fig; kwargs...) = savefig("", basename, fig; kwargs...)
function savefig(subdir, prefix, gm::AbstractGraphModel, fig; kwargs...)
    savefig(subdir,
        savename(prefix, gm;
            allowedtypes=[DrWatson.default_allowed(gm)...; Dict],
            expand=["metadata"]
        ), fig; kwargs...)
end
savefig(prefix, gm::AbstractGraphModel, fig; kwargs...) = savefig("", prefix, gm, fig; kwargs...)
export savefig

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
export save_transmat

"""
Saves any julia object in a jld2 file in the data dir with a fancy name.
"""
function save_object_smart(basename::AbstractString, obj, assoc_objs...; subdir="", separator="__")
    filepath = datadir(subdir, basename * separator * join(savename.([obj; collect(assoc_objs)]), separator) * ".jld2")
    save_object(filepath, obj)
end
export save_object_smart
