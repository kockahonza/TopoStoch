################################################################################
# Complex plotting ComplexAllosteryGM objects
################################################################################
function p_do_layout(ca::ComplexAllosteryGM, layout, roffset_devs)
    if isnothing(layout)
        layout = :NR
    end
    if isnothing(roffset_devs)
        roffset_devs = :auto
    end

    # this may or may not be filled
    axis_labels = nothing

    if isa(layout, NetworkLayout.AbstractLayout)
        layout = layout(p_get_adjmat_symsafe(ca))
        dim = length(layout[1])
    elseif isa(layout, Function)
        layout = layout(ca)
        dim = length(layout[1])
    elseif isa(layout, AbstractArray) && (length(layout) == numstates(ca))
        dim = length(layout[1])
    else
        if isa(layout, Symbol)
            layout = (layout,)
        end
        (type, layout_args...) = layout

        # These are mostly the case, override where not
        do_default_roffsets = true
        dim = 3

        if type == :partlyc3d
            sub_layout, z_scale = layout_args
            layout = [
                (x, y, z * z_scale) for ((x, y), z) in zip(sub_layout(p_get_adjmat_symsafe(ca)), calc_numligands.(allstates(ca)))
            ]
            do_default_roffsets = false
        elseif type == :NEB
            layout = [
                (calc_numligands(st), calc_energy(st, ca), calc_numboundaries(st, ca)) for st in allstates(ca)
            ]
            maxs = maximum(layout)
            axis_labels = (f"N({maxs[1]})", f"E({maxs[2]})", f"boundaries({maxs[3]})")
        elseif type == :NRB
            layout = [
                Point{3,Float64}(calc_numligands(st), calc_numofconf(st), calc_numboundaries(st, ca)) for st in allstates(ca)
            ]
            maxs = maximum(layout)
            axis_labels = (f"N({maxs[1]})", f"#tense({maxs[2]})", f"#boundaries({maxs[3]})")
        elseif type == :NB
            layout = [
                (Float64(calc_numligands(st)), Float64(calc_numboundaries(st, ca))) for st in allstates(ca)
            ]
            maxs = maximum(layout)
            axis_labels = (f"N({maxs[1]})", f"boundaries({maxs[2]})")
            dim = 2
        elseif type == :NR
            layout = [
                (Float64(calc_numligands(st)), Float64(calc_numofconf(st))) for st in allstates(ca)
            ]
            maxs = maximum(layout)
            axis_labels = (f"N({maxs[1]})", f"Number of tense({maxs[2]})")
            dim = 2
        elseif type == :N
            layout = [
                (Float64(calc_numligands(st)), 0.0) for st in allstates(ca)
            ]
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

# This is the main logic function
function plot_ca_!(ax, ca::ComplexAllosteryGM{S,F},
    do_layout_rslt, args...;
    axparent=nothing,
    # these apply to all
    flabels=:auto,
    fnlabels=nothing,
    felabels=nothing,
    node_size_scale=20.0,
    interactive=:auto,
    symbolic=:auto,
    # these only to concrete (non symbolic) cas
    c_colormap=:dense,
    c_colorscale=:auto,
    c_colorrange=:auto,
    c_colorbar=:auto,
    kwargs...
) where {S,F}
    layout, dim, axis_labels = do_layout_rslt
    if symbolic == :auto
        symbolic = F == Num
    end

    auto_kwargs = Dict{Symbol,Any}()

    # This is for possible later elementwise updates
    auto_kwargs[:node_color] = [:snow3 for _ in 1:numstates(ca)]
    auto_kwargs[:node_size] = [node_size_scale for _ in 1:numstates(ca)]

    # Do colored transitions
    if !symbolic
        auto_kwargs[:edge_color] = weight.(edges(ca.graph))
        if c_colorrange == :auto
            max_weight = maximum(weight.(edges(ca.graph)))
            if iszero(max_weight)
                max_weight = 1.0
            end
            c_colorrange = (0.0, max_weight)
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
        flabels = nv(ca.graph) < 100
        if isnothing(fnlabels)
            fnlabels = :repr
        end
        if isnothing(felabels)
            felabels = :repr
        end
    end
    if flabels
        if fnlabels == :repr
            auto_kwargs[:nlabels] = repr.(allstates(ca))
        elseif fnlabels == :E
            auto_kwargs[:nlabels] = repr.(calc_energy.(allstates(ca), ca))
        end
        if felabels == :repr
            auto_kwargs[:elabels] = repr.(weight.(edges(ca.graph)))
            auto_kwargs[:elabels_shift] = 0.4
        end
    end

    plot = graphplot!(ax, ca.graph, args...;
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

# Simple callers
function plot_ca!(maybeax, ca::ComplexAllosteryGM, args...;
    map=nothing, layout=nothing, roffset_devs=nothing, returnax=false, kwargs...
)
    if !isnothing(map)
        ca = map(ca)
    end

    do_layout_rslt = p_do_layout(ca, layout, roffset_devs)

    if isa(maybeax, Makie.AbstractAxis)
        ax = maybeax
    elseif isa(maybeax, GridPosition)
        ax = p_make_ax(do_layout_rslt.dim, maybeax)
    elseif isa(maybeax, GridLayout)
        ax = p_make_ax(do_layout_rslt.dim, maybeax[1, 1])
    else
        throw(ArgumentError(f"maybeax type {typeof(maybeax)} is not recognised, must be either Makie.AbstractAxis or GridPosition"))
    end

    plot = plot_ca_!(ax, ca, do_layout_rslt, args...; kwargs...)

    if returnax
        ax, plot
    else
        plot
    end
end
function plot_ca(ca::ComplexAllosteryGM, args...;
    map=nothing, layout=nothing, roffset_devs=nothing, kwargs...
)
    if !isnothing(map)
        ca = map(ca)
    end

    do_layout_rslt = p_do_layout(ca, layout, roffset_devs)

    fig = Figure()
    ax = p_make_ax(do_layout_rslt.dim, fig[1, 1])

    plot = plot_ca_!(ax, ca, do_layout_rslt, args...; axparent=fig, kwargs...)

    Makie.FigureAxisPlot(fig, ax, plot)
end
# alias used by simGM
# function plotgm(ca::ComplexAllosteryGM, args...; kwargs...)
#     plot_ca(ca, args...; kwargs...)
# end

"""
Very simple and not so efficient caller that can use a graph observable
"""
function plot_ca!(maybeax, cao::Observable, args...;
    returnax=false, remove_colorbars=true, refresh_layout_obs=nothing, kwargs...
)
    ax, plot = plot_ca!(maybeax, cao.val, args...; returnax=true, kwargs...)

    on(cao) do ca
        layout = plot.node_pos[]
        delete!(ax, plot)
        if remove_colorbars
            delete!.(filter(x -> isa(x, Colorbar), ax.parent.content))
        end
        plot = plot_ca!(ax, ca, args...; kwargs..., layout)
    end

    if !isnothing(refresh_layout_obs)
        on(refresh_layout_obs) do _
            delete!(ax, plot)
            if remove_colorbars
                delete!.(filter(x -> isa(x, Colorbar), ax.parent.content))
            end
            plot = plot_ca!(ax, cao.val, args...; kwargs...)
        end
    end

    if returnax
        ax, plot
    else
        plot
    end
end

################################################################################
# More interactivity
################################################################################
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

function int_plot_ca(ca::ComplexAllosteryGM, args...;
    variables=Num.(get_variables(ca)),
    doeig=true, noplot=false,
    ranges=nothing, startvalues=nothing, kwargs...
)
    # Setup sliders and the graph observable
    fig, sobs, sg = prep_sliders(variables; ranges, startvalues)
    ca_factory = make_factory(ca, variables)

    ca_obs = lift(sobs...) do (svals...)
        ca_factory(svals...)
    end

    if noplot
        return FigureAxisAnything(fig, nothing, ca_obs)
    end

    buttons = fig[1, 2] = GridLayout()
    colsize!(fig.layout, 2, Fixed(60))
    refresh_layout_btn = Button(buttons[1, 1], label="layout")

    plotspace = fig[2, :] = GridLayout()

    extra_kwargs = Dict{Symbol,Any}()
    if doeig
        extra_kwargs[:node_size] = lift(supersteadystate, ca_obs)
    end

    # Make a suitable ax the the graph plotting
    ax, plot = plot_ca!(plotspace, ca_obs, args...;
        returnax=true,
        axparent=plotspace,
        refresh_layout_obs=refresh_layout_btn.clicks,
        roffset_devs=false,
        extra_kwargs...,
        kwargs...
    )

    FigureAxisAnything(fig, ax, ca_obs)
end

function int_eigen(ca::ComplexAllosteryGM, args...;
    variables=Num.(get_variables(ca)),
    ranges=nothing, startvalues=nothing, kwargs...
)
    # Setup sliders and the graph observable
    fig, sobs, _ = prep_sliders(variables; ranges, startvalues)
    ca_factory = make_factory(ca, variables)

    ca_obs = lift(sobs...) do (svals...)
        ca_factory(svals...)
    end

    buttons = fig[1, 2] = GridLayout()
    colsize!(fig.layout, 2, Fixed(60))
    refresh_layout_btn = Button(buttons[1, 1], label="layout")

    plotspace = fig[2, :] = GridLayout()

    # Make a suitable ax the the graph plotting
    ax, plot = plot_ca!(plotspace, ca_obs, args...; returnax=true, axparent=plotspace, refresh_layout_obs=refresh_layout_btn.clicks, roffset_devs=false, kwargs...)

    FigureAxisAnything(fig, ax, ca_obs)
end

# other interactive plots
function eq_stats_plot(ca::ComplexAllosteryGM{S,Num}) where {S}
    fig = Figure()

    variables = Num.(get_variables(ca.energy_matrices))
    @variables μ, kT
    append!(variables, μ)

    sliders = [(label=repr(var), range=0.0:0.01:10.0, startvalue=1.0) for var in variables]
    sg = SliderGrid(fig[1, 1], sliders...)
    variable_observables = [x.value for x in sg.sliders]

    overall_stats = GridLayout()
    fig[1, 2] = overall_stats
    colsize!(fig.layout, 2, 30)

    partition_function_f = build_function(substitute(calc_partition_function(ca), (kT => 1)), variables; expression=Val(false))

    avg_numligands_f = build_function(substitute(calc_avg_numligands(ca), (kT => 1)), variables; expression=Val(false))
    avg_numligands_o = lift((args...) -> f"N={avg_numligands_f(args):.3g}", variable_observables...)
    Label(overall_stats[1, 1], avg_numligands_o)

    avg_energy_f = build_function(substitute(calc_avg_energy(ca), (kT => 1)), variables; expression=Val(false))
    avg_energy_o = lift((args...) -> f"E={avg_energy_f(args):.3g}", variable_observables...)
    Label(overall_stats[2, 1], avg_energy_o)

    # Do the barplot
    state_reordering = sortperm(calc_numligands.(allstates(ca)); alg=Base.Sort.DEFAULT_STABLE)
    state_labels = repr.(allstates(ca))[state_reordering]
    ax = Axis(fig[2, :];
        xticks=(1:numstates(ca), state_labels),
        xticklabelrotation=pi / 5
    )
    ylims!(ax, 0.0, 0.5)

    gfs = substitute.(calc_gibbs_factor.(allstates(ca), ca), (kT => 1.0))
    gfs_fs = build_function.(gfs, Ref(variables); expression=Val(false))

    probabilities = lift(variable_observables...) do (vars...)
        pfval = partition_function_f(vars)
        [gfs_fs[i](vars) / pfval for i in 1:numstates(ca)][state_reordering]
    end

    bp = barplot!(ax, 1:numstates(ca), probabilities)

    Makie.FigureAxisPlot(fig, ax, bp)
end

function eq_coopbinding_plot(ca::ComplexAllosteryGM{S,Num}) where {S}
    fig = Figure()

    energy_vars = collect(get_variables(ca.energy_matrices))
    @variables μ, εsol, kT, c

    slider_vars = [energy_vars; εsol]

    sliders = [(label=repr(var), range=-5.0:0.1:10.0, startvalue=1.0) for var in slider_vars]
    sg = SliderGrid(fig[1, 1], sliders...)
    slider_observables = [x.value for x in sg.sliders]

    ax = Axis(fig[2, :])

    prepped_expression = substitute(calc_avg_numligands(ca), Dict(
        kT => 1,
        μ => εsol + log(c)
    ))
    avg_numligands_f = build_function(prepped_expression, [energy_vars; εsol; c]; expression=Val(false))

    conc_range = LinRange(0, 200, 5000)

    avg_numligands = lift(slider_observables...) do (vars...)
        [avg_numligands_f((vars..., cval)) for cval in conc_range]
    end

    sc = scatter!(ax, conc_range, avg_numligands)

    Makie.FigureAxisPlot(fig, ax, sc)
end

################################################################################
# Saving data
################################################################################
function save_transmat(ca::ComplexAllosteryGM; name=savename("transmat", ca, "table"), short=false)
    file = open(datadir(name), "w")

    row_names = repr.(1:numstates(ca)) .* " / " .* repr.(allstates(ca))
    tm = transmat(ca; mat=true)
    if short
        header = [""; repr.(1:numstates(ca))]
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

function save_ca_obj(ca::ComplexAllosteryGM; kwargs...)
    for (name, value) in kwargs
        save_object(datadir(savename(String(name), ca, "jld2")), value)
    end
end

function save_evals(ca::ComplexAllosteryGM, evals; name=savename("evals", ca))
    file = open(datadir(name), "w")
    println(file, "evals:")
    for eval in evals
        println(file, eval)
    end
    close(file)
end

function save_gibbs_factors(ca::ComplexAllosteryGM; name=savename("adjmat", ca))
    file = open(datadir(savename("gibbsfs", ca, "table")), "w")
    data = [repr.(allstates(ca)) repr.(calc_gibbs_factor.(allstates(ca), ca))]
    pretty_table(file, data; header=["state", "gibbs factor"])
    close(file)
end

# Symbolics helper functions
function get_rewriter()
    rule_divexp = @rule ~y / exp(~x) => ~y * exp(-~x)
    rules_explog = [
        (@rule exp(log(~z)) => ~z),
        (@rule exp(~y * log(~z)) => ~z^~y),
        (@rule exp(~x + ~y * log(~z)) => exp(~x) * ~z^~y)
    ]

    Rewriters.Prewalk(Rewriters.Chain([
        rule_divexp,
        rules_explog...,
        SymbolicUtils.default_simplifier()
    ]))
end
