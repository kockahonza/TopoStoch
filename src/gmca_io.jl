# New plotting interface
function p_named_layouts(ca::ComplexAllosteryGM, layout_name, layout_args)
    try
        invoke(p_named_layouts, Tuple{supertype(typeof(ca)),Any,Any}, ca, layout_name, layout_args)
    catch ArgumentError
        def_roffsets = true
        dim = 3
        axis_labels = nothing

        if layout_name == :partlyc3d
            sub_layout, z_scale = layout_args
            layout = [
                (x, y, z * z_scale) for ((x, y), z) in zip(sub_layout(p_safe_adjmat(ca)), calc_numligands.(allstates(ca)))
            ]
            def_roffsets = false
        elseif layout_name == :NEB
            layout = [
                (calc_numligands(st), calc_energy(st, ca), calc_numboundaries(st, ca)) for st in allstates(ca)
            ]
            maxs = maximum(layout)
            axis_labels = (f"N({maxs[1]})", f"E({maxs[2]})", f"boundaries({maxs[3]})")
        elseif layout_name == :NRB
            layout = [
                Point{3,Float64}(calc_numligands(st), calc_numofconf(st), calc_numboundaries(st, ca)) for st in allstates(ca)
            ]
            maxs = maximum(layout)
            axis_labels = (f"N({maxs[1]})", f"#tense({maxs[2]})", f"#boundaries({maxs[3]})")
        elseif layout_name == :NB
            layout = [
                (Float64(calc_numligands(st)), Float64(calc_numboundaries(st, ca))) for st in allstates(ca)
            ]
            maxs = maximum(layout)
            axis_labels = (f"N({maxs[1]})", f"boundaries({maxs[2]})")
            dim = 2
        elseif layout_name == :NR
            layout = [
                (Float64(calc_numligands(st)), Float64(calc_numofconf(st))) for st in allstates(ca)
            ]
            maxs = maximum(layout)
            axis_labels = (f"N({maxs[1]})", f"Number of tense({maxs[2]})")
            dim = 2
        elseif layout_name == :N
            layout = [
                (Float64(calc_numligands(st)), 0.0) for st in allstates(ca)
            ]
            maxs = maximum(layout)
            axis_labels = (f"N({maxs[1]})", "")
            dim = 2
        else
            rethrow()
        end

        dim, layout, def_roffsets, axis_labels
    end
end
function p_do_layout(ca::ComplexAllosteryGM, layout=nothing, roffset_devs=nothing)
    if isnothing(layout)
        layout = :NR
    end
    invoke(p_do_layout, Tuple{AbstractGraphModel,Any,Any}, ca, layout, roffset_devs)
end

# FIX: Everything from now on should be replaced by new stuff in gm_io.jl,
# keeping for now for observables stuff and interactivity
################################################################################
# More interactivity
################################################################################
"""
Very simple and not so efficient caller that can use a graph observable
"""
function plot_cao!(maybeax, cao::Observable, args...;
    returnax=false, remove_colorbars=true, refresh_layout_obs=nothing, kwargs...
)
    ax, plot = plotgm!(maybeax, cao.val, args...; returnax=true, kwargs...)

    on(cao) do ca
        layout = plot.node_pos[]
        delete!(ax, plot)
        if remove_colorbars
            delete!.(filter(x -> isa(x, Colorbar), ax.parent.content))
        end
        plot = plotgm!(ax, ca, args...; kwargs..., layout)
    end

    if !isnothing(refresh_layout_obs)
        on(refresh_layout_obs) do _
            delete!(ax, plot)
            if remove_colorbars
                delete!.(filter(x -> isa(x, Colorbar), ax.parent.content))
            end
            plot = plotgm!(ax, cao.val, args...; kwargs...)
        end
    end

    if returnax
        ax, plot
    else
        plot
    end
end

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
    noplot=false, ranges=nothing, startvalues=nothing, kwargs...
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

    # Make a suitable ax the the graph plotting
    ax, plot = plot_cao!(plotspace, ca_obs, args...;
        returnax=true,
        axparent=plotspace,
        refresh_layout_obs=refresh_layout_btn.clicks,
        roffset_devs=false,
        kwargs...
    )

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
function save_gibbs_factors(ca::ComplexAllosteryGM; name=savename("adjmat", ca))
    file = open(datadir(savename("gibbsfs", ca, "table")), "w")
    data = [repr.(allstates(ca)) repr.(calc_gibbs_factor.(allstates(ca), ca))]
    pretty_table(file, data; header=["state", "gibbs factor"])
    close(file)
end