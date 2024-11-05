################################################################################
# Complex plotting ComplexAllosteryGM objects
################################################################################
function p_get_adjmat_symsafe(ca::ComplexAllosteryGM{S,F}) where {S,F}
    if F == Num
        map(x -> if iszero(x)
                0.0
            else
                1.0
            end, adjacency_matrix(ca.graph))
    else
        adjacency_matrix(ca.graph)
    end
end

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

function p_make_ca_ax(dim, place)
    if dim == 2
        Axis(place)
    elseif dim == 3
        LScene(place)
    else
        throw(ErrorException("the layout dimension was not 2 or 3"))
    end
end

# This is the main logic function
function plot_ca_!(ax, ca::ComplexAllosteryGM{S,F},
    do_layout_rslt, args...;
    axparent=nothing,
    # these apply to all
    fancy_labels=:auto,
    node_size_scale=10.0,
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

    auto_kwargs[:node_color] = [:snow3 for _ in 1:numstates(ca)]
    auto_kwargs[:node_size] = [node_size_scale for _ in 1:numstates(ca)]

    if !symbolic
        # color transitions smartly
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

    if fancy_labels == :auto
        fancy_labels = symbolic && (nv(ca.graph) < 50)
    end
    if fancy_labels
        auto_kwargs[:nlabels] = repr.(allstates(ca))
        auto_kwargs[:elabels] = repr.(weight.(edges(ca.graph)))
        auto_kwargs[:elabels_shift] = 0.4
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
    layout=nothing, roffset_devs=nothing, returnax=false, kwargs...
)
    do_layout_rslt = p_do_layout(ca, layout, roffset_devs)

    if isa(maybeax, Makie.AbstractAxis)
        ax = maybeax
    elseif isa(maybeax, GridPosition)
        ax = p_make_ca_ax(do_layout_rslt.dim, maybeax)
    elseif isa(maybeax, GridLayout)
        ax = p_make_ca_ax(do_layout_rslt.dim, maybeax[1, 1])
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
    layout=nothing, roffset_devs=nothing, kwargs...
)
    do_layout_rslt = p_do_layout(ca, layout, roffset_devs)

    fig = Figure()
    ax = p_make_ca_ax(do_layout_rslt.dim, fig[1, 1])

    plot = plot_ca_!(ax, ca, do_layout_rslt, args...; axparent=fig, kwargs...)

    Makie.FigureAxisPlot(fig, ax, plot)
end
# alias used by simGM
function plotGM(ca::ComplexAllosteryGM, args...; kwargs...)
    plot_ca(ca, args...; kwargs...)
end

"""
Very simple and not so efficient caller that can use a graph observable
"""
function plot_ca!(maybeax, cao::Observable, args...;
    returnax=false, remove_colorbars=true, kwargs...
)
    ax, plot = plot_ca!(maybeax, cao.val, args...; returnax=true, kwargs...)

    on(cao) do ca
        empty!(ax)
        if remove_colorbars
            delete!.(filter(x -> isa(x, Colorbar), ax.parent.content))
        end
        plot_ca!(ax, ca, args...; kwargs...)
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
function prep_sliders(variables; ranges=nothing)
    slider_specs = []
    for i in eachindex(variables)
        range_ = isnothing(ranges) ? (0.0:0.1:10.0) : ranges[i]
        push!(slider_specs, (; label=string(variables[i]), range=range_))
    end

    fig = Figure()
    sg = SliderGrid(fig[1, 1], slider_specs...)

    fig, [x.value for x in sg.sliders], sg
end

function int_plot_ca(ca::ComplexAllosteryGM, args...;
    variables=Num.(get_variables(ca)), ranges=nothing, kwargs...
)
    # Setup sliders and the graph observable
    fig, sobs, _ = prep_sliders(variables; ranges)
    ca_factory = make_factory(ca, variables)

    ca_obs = lift(sobs...) do (svals...)
        ca_factory(svals...)
    end

    plotspace = fig[2, 1] = GridLayout()

    # Make a suitable ax the the graph plotting
    ax, plot = plot_ca!(plotspace, ca_obs, args...; returnax=true, axparent=plotspace, kwargs...)

    Makie.FigureAxisPlot(fig, ax, plot)
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
    am = transmat(ca; mat=true)
    if short
        header = [""; repr.(1:numstates(ca))]
        am = map(x -> if iszero(x)
                0
            else
                1
            end, am)
    else
        header = [""; row_names]
    end
    modified_matrix = [row_names am]

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
