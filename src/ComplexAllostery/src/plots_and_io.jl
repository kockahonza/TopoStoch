################################################################################
# Main plotting - add methods to GraphModels functions
################################################################################
import GraphModels: p_named_layouts, p_do_layout, plotgm_kwarg_defaults
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
            axis_labels = (f"#ligands({maxs[1]})", f"E({maxs[2]})", f"boundaries({maxs[3]})")
        elseif layout_name == :NTB
            layout = [
                Point{3,Float64}(calc_numligands(st), calc_numofconf(st), calc_numboundaries(st, ca)) for st in allstates(ca)
            ]
            maxs = maximum(layout)
            axis_labels = (f"#ligands({maxs[1]})", f"#tense({maxs[2]})", f"#boundaries({maxs[3]})")
        elseif layout_name == :NB
            layout = [
                (Float64(calc_numligands(st)), Float64(calc_numboundaries(st, ca))) for st in allstates(ca)
            ]
            maxs = maximum(layout)
            axis_labels = (f"#ligands({maxs[1]})", f"#boundaries({maxs[2]})")
            dim = 2
        elseif layout_name == :NT
            layout = [
                (Float64(calc_numligands(st)), Float64(calc_numofconf(st))) for st in allstates(ca)
            ]
            maxs = maximum(layout)
            axis_labels = (f"#ligands({maxs[1]})", f"#tense({maxs[2]})")
            dim = 2
        elseif layout_name == :N
            layout = [
                (Float64(calc_numligands(st)), 0.0) for st in allstates(ca)
            ]
            maxs = maximum(layout)
            axis_labels = (f"#ligands({maxs[1]})", "")
            dim = 2
        elseif layout_name == :NTc
            if ca.C != 2
                @warn "named layout NTc is specific to C=2, this will likely not work well"
            end

            groups = Matrix{Vector{Int}}(undef, ca.N * ca.B + 1, ca.N + 1)
            for i in eachindex(groups)
                groups[i] = Int[]
            end
            for (i, st) in enumerate(allstates(ca))
                numlig = calc_numligands(st)
                numtense = calc_numofconf(st)
                push!(groups[numlig+1, numtense+1], i)
            end

            layout = Vector{Point}(undef, numstates(ca))
            R = (length(layout_args) >= 1) ? layout_args[1] : 0.3
            for i in axes(groups, 1)
                for j in axes(groups, 2)
                    numlig = i - 1
                    numtense = j - 1
                    group = groups[i, j]
                    Ngroup = length(group)
                    for (group_i, global_i) in enumerate(group)
                        theta = ((group_i - 1) / Ngroup) * 2 * pi
                        layout[global_i] = Point2f(numlig, numtense)
                        if Ngroup != 1
                            layout[global_i] += Point2f(R * cos(theta), R * sin(theta))
                        end
                    end
                end
            end
            maxs = maximum(layout)
            axis_labels = (f"#ligands({maxs[1]})", f"#tense({maxs[2]})")
            dim = 2
            def_roffsets = false
        else
            rethrow()
        end

        dim, layout, def_roffsets, axis_labels
    end
end
function p_do_layout(ca::ComplexAllosteryGM, layout=nothing, roffset_devs=nothing)
    if isnothing(layout)
        layout = :NTc
    end
    invoke(p_do_layout, Tuple{AbstractGraphModel,Any,Any}, ca, layout, roffset_devs)
end

function plotgm_kwarg_defaults(ca::ComplexAllosteryGM{S,Num}) where {S}
    fnlabels = if numstates(ca) < 10
        :repr
    else
        labelstate_NR
    end
    (; fnlabels)
end
function plotgm_kwarg_defaults(ca::ComplexAllosteryGM{S,<:AbstractFloat}) where {S}
    if numstates(ca) < 10
        fnlabels = :repr
        felabels = true
    else
        fnlabels = labelstate_NR
        felabels = false
    end
    (;
        flabels=true,
        fnlabels,
        felabels,
        n_size=15.0,
        n_ss_size=false,
        n_ss_colorscale=log10,
        e_color=:dcurrents
    )
end

################################################################################
# ComplexAllostery specific plots
################################################################################
function plotca_macro(
    ca::ComplexAllosteryGM{S,F},
    ss=supersteadystate(ca);
    layout=nothing,
    n_size=20.0,
    n_ss_size=false,
    n_ss_colormap=:Oranges,
    n_ss_colorscale=log10,
    n_ss_colorrange=:auto,
    n_ss_colorbar=true,
    e_color=:dcurrents,
    e_colormap=:BuPu,
    e_colorscale=identity,
    e_colorrange=:auto,
    e_colorbar=true,
    interactive=:auto,
    kwargs...
) where {S,F<:AbstractFloat}
    group_vect = [(calc_numligands(st), calc_numofconf(st)) for st in allstates(ca)]
    groups = unique(group_vect)

    graph_ = if e_color == :rates
        graph(ca)
    elseif e_color == :dcurrents
        make_current_graph(ca, ss)
    else
        throw(ArgumentError(f"unrecognized e_color of {e_color}"))
    end

    macrograph, macrostates = groupgraph(graph_, group_vect; groups)
    macrostate, _ = groupsum(ss, group_vect; groups)

    if isnothing(layout)
        layout = :NT
    end
    if isa(layout, NetworkLayout.AbstractLayout)
        layout = layout(macrograph)
    elseif isa(layout, AbstractArray) && (length(layout) == nv(macrograph))
    else
        if isa(layout, Symbol)
            layout = (layout,)
        end
        (layout_name, layout_args...) = layout
        if layout_name == :NT
            layout = macrostates
        else
            throw(ArgumentError(f"layout of name {layout_name} is not recognised"))
        end
    end

    auto_kwargs = Dict{Symbol,Any}()

    # node size
    auto_kwargs[:node_size] = if !(isnothing(n_ss_size) || (n_ss_size == false))
        n_size *= 4.0
        if n_ss_size == :lin
            [n_size * x for x in macrostate]
        elseif n_ss_size == :log
            [n_size * log(1 + x) for x in macrostate]
        elseif n_ss_size == :sqrt
            [n_size * sqrt(x) for x in macrostate]
        elseif n_ss_size == :linA
            [n_size * (x)^(1 / 2) for x in macrostate]
        elseif n_ss_size == :logA
            [n_size * (log(1 + x))^(1 / 2) for x in macrostate]
        elseif n_ss_size == :sqrtA
            [n_size * sqrt(x)^(1 / 2) for x in macrostate]
        else
            throw(ArgumentError(f"n_ss_size of {n_ss_size} is not recognised"))
        end
    else
        [n_size for _ in 1:nv(macrograph)]
    end
    # node colors
    auto_kwargs[:node_color] = macrostate
    if n_ss_colorrange == :auto
        max_weight = maximum(macrostate)
        min_weight = minimum(macrostate)
        if (max_weight - min_weight) < 1e-8
            min_weight -= 0.1
        end
        n_ss_colorrange = (min_weight, max_weight)
    end
    auto_kwargs[:node_attr] = (;
        colormap=n_ss_colormap,
        colorrange=n_ss_colorrange,
        colorscale=n_ss_colorscale
    )
    # edge colors
    auto_kwargs[:edge_color] = colors = weight.(edges(macrograph))
    if e_colorrange == :auto
        max_weight = maximum(colors)
        min_weight = minimum(colors)
        if (max_weight - min_weight) < 1e-8
            min_weight -= 0.1
        end
        e_colorrange = (min_weight, max_weight)
    end
    auto_kwargs[:edge_attr] = auto_kwargs[:arrow_attr] = (;
        colormap=e_colormap,
        colorrange=e_colorrange,
        colorscale=e_colorscale
    )

    fap = graphplot(macrograph;
        layout,
        auto_kwargs...
    )

    if e_colorbar
        Colorbar(fap.figure[1, 2];
            colormap=e_colormap,
            colorrange=e_colorrange,
            scale=e_colorscale,
            label="net probability current"
        )
    end
    if n_ss_colorbar
        Colorbar(fap.figure[1, 3];
            colormap=n_ss_colormap,
            colorrange=n_ss_colorrange,
            scale=n_ss_colorscale,
            label="probability"
        )
    end

    if interactive == :auto
        interactive = isa(fap.axis, Axis)
    end
    if interactive
        GraphModels.p_make_interactive(fap.axis, fap.plot)
    end

    fap
end
export plotca_macro


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
        plot = plotgm!(ax, ca, args...; kwargs..., layout, roffset_devs=false)
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

function int_plot_ca(ca::ComplexAllosteryGM, args...;
    variables=get_variables(ca),
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

################################################################################
# Other interactive plots -- mostly old stuff but works and interesting!
################################################################################
function eq_stats_plot(ca::ComplexAllosteryGM{S,Num}) where {S}
    fig = Figure()

    variables = get_variables(ca.energy_matrices)
    mu, kT = symvar_mu, symvar_kT
    append!(variables, mu)

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
export eq_stats_plot

function eq_coopbinding_plot(ca::ComplexAllosteryGM{S,Num}) where {S}
    fig = Figure()

    energy_vars = get_variables(ca.energy_matrices)
    mu, kT = symvar_mu, symvar_kT
    @variables εsol, c

    slider_vars = [energy_vars; εsol]

    sliders = [(label=repr(var), range=-5.0:0.1:10.0, startvalue=1.0) for var in slider_vars]
    sg = SliderGrid(fig[1, 1], sliders...)
    slider_observables = [x.value for x in sg.sliders]

    ax = Axis(fig[2, :])

    prepped_expression = substitute(calc_avg_numligands(ca), Dict(
        kT => 1,
        mu => εsol + log(c)
    ))
    avg_numligands_f = build_function(prepped_expression, [energy_vars; εsol; c]; expression=Val(false))

    conc_range = LinRange(0, 200, 5000)

    avg_numligands = lift(slider_observables...) do (vars...)
        [avg_numligands_f((vars..., cval)) for cval in conc_range]
    end

    sc = scatter!(ax, conc_range, avg_numligands)

    Makie.FigureAxisPlot(fig, ax, sc)
end
export eq_coopbinding_plot

################################################################################
# Saving data
################################################################################
function save_gibbs_factors(ca::ComplexAllosteryGM; name=savename("adjmat", ca))
    file = open(datadir(savename("gibbsfs", ca, "table")), "w")
    data = [repr.(allstates(ca)) repr.(calc_gibbs_factor.(allstates(ca), ca))]
    pretty_table(file, data; header=["state", "gibbs factor"])
    close(file)
end
export save_gibbs_factors
