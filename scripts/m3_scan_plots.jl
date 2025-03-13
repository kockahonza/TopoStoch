using DrWatson
@quickactivate "TopoStochSim"

using JLD2
using GLMakie
using Printf
using MappedArrays

using GraphModels # not needed, just here to help the LSP
using ComplexAllostery
using ComplexAllostery.Model3

function scansummary(f::JLD2.JLDFile)
    range_names = f["range_names"]
    for rn in range_names
        @printf "%-10s - %3d\n" rn length(f[rn])
    end
    println(f["observable_names"])
end
scansummary(fn::AbstractString) = scansummary(jldopen(fn))

function _handle_observable(observable_names, obs)
    if isa(obs, Tuple)
        if length(obs) == 2
            obs, func = obs
            label = nothing
        elseif length(obs) == 3
            obs, func, label = obs
        else
            throw(ArgumentError(""))
        end
    else
        func = nothing
        label = nothing
    end
    if !isa(obs, Integer)
        maybe_obs = findfirst(x -> x == obs, observable_names)
        if !isnothing(maybe_obs)
            obs = maybe_obs
        else
            throw(ArgumentError(@sprintf "%s is an unrecognised observable name" string(obs)))
        end
    end
    obs = observable_names[obs]
    if isnothing(label)
        label = obs
    end
    (obs, label, func)
end

function _handle_rangei(range_names, i)
    if isa(i, AbstractString)
        maybe_xi = findfirst(x -> x == i, range_names)
        if !isnothing(maybe_xi)
            i = maybe_xi
        else
            throw(ArgumentError(@sprintf "cannot find range of name %s" i))
        end
    elseif !isa(i, Integer)
        throw(ArgumentError(@sprintf "invalid range identifier %s" string(i)))
    end
    if (i < 1) || (i > length(range_names))
        throw(ArgumentError(@sprintf "range index %d out of range" i))
    end
    i
end

function _centres_to_edges_clipped(range)
    edges = Vector{Float64}(undef, length(range) + 1)
    edges[1] = range[1]
    for i in 1:(length(range)-1)
        edges[i+1] = (range[i+1] + range[i]) / 2
    end
    edges[end] = range[end]
    edges
end

function scanheatmaps(f, xi, yi, observables;
    nrows=nothing, xscale=identity, yscale=identity, kwargs...
)
    # prep and handle inputs
    range_names = f["range_names"]
    observable_names = f["observable_names"]
    xi = _handle_rangei(range_names, xi)
    yi = _handle_rangei(range_names, yi)
    if xi == yi
        throw(ArgumentError("The two range variables cannot be the same"))
    end
    observables = [_handle_observable(observable_names, obs) for obs in observables]

    # make fig and more prep
    fig = Figure()
    nplots = length(observables)
    if isnothing(nrows)
        (nrows, ncols) = ntogridsize1(nplots)
    else
        ncols = ceil(Int, nplots / nrows)
    end

    # make sliders
    varis = filter(x -> !(x in (xi, yi)), 1:length(range_names))
    slider_specs = []
    for vari in varis
        range_name = range_names[vari]
        push!(slider_specs, (label=range_name, range=f[range_name]))
    end
    sg = SliderGrid(fig[1, 1:ncols*2], slider_specs...)
    var_obs = [x.selected_index for x in sg.sliders]

    # this does the complex bit
    function make_data_obs(data)
        lift(var_obs...) do (varvals...)
            indices = Any[(:) for _ in 1:length(range_names)]
            for (vari, varval) in zip(varis, varvals)
                indices[vari] = varval
            end
            if xi < yi
                data[indices...]
            else
                transpose(data[indices...])
            end
        end
    end

    # Do the plots
    xs = _centres_to_edges_clipped(f[range_names[xi]])
    ys = _centres_to_edges_clipped(f[range_names[yi]])
    cis = CartesianIndices((ncols, nrows))
    axs = []
    hms = []
    for ((obs_name, obs_label, func), ci) in zip(observables, cis)
        ax = Axis(fig[1+ci[2], ci[1]*2-1]; xscale, yscale)
        push!(axs, ax)

        data_obs = if isnothing(func)
            make_data_obs(f[obs_name])
        else
            make_data_obs(mappedarray(func, f[obs_name]))
        end

        hm = heatmap!(ax, xs, ys, data_obs; kwargs...)
        push!(hms, hm)
        Colorbar(fig[1+ci[2], ci[1]*2], hm)

        ax.title = obs_label
        ax.xlabel = range_names[xi]
        ax.ylabel = range_names[yi]
    end

    return FigureAxisAnything(fig, axs, hms)
end

function scanmultiscatters(f, xi, observables;
    nrows=nothing, xscale=identity, yscale=identity, axis=(;), kwargs...
)
    # prep and handle inputs
    range_names = f["range_names"]
    observable_names = f["observable_names"]
    xi = _handle_rangei(range_names, xi)
    observables = [_handle_observable(observable_names, obs) for obs in observables]

    # make fig and more prep
    fig = Figure()
    nplots = length(observables)
    if isnothing(nrows)
        (nrows, ncols) = ntogridsize1(nplots)
    else
        ncols = ceil(Int, nplots / nrows)
    end

    # make sliders
    varis = filter(x -> x != xi, 1:length(range_names))
    slider_specs = []
    for vari in varis
        range_name = range_names[vari]
        push!(slider_specs, (label=range_name, range=f[range_name]))
    end
    sg = SliderGrid(fig[1, 1:ncols], slider_specs...)
    var_obs = [x.selected_index for x in sg.sliders]

    # this does the complex bit
    function make_data_obs(data)
        lift(var_obs...) do (varvals...)
            indices = Any[(:) for _ in 1:length(range_names)]
            for (vari, varval) in zip(varis, varvals)
                indices[vari] = varval
            end
            data[indices...]
        end
    end

    # Do the plots
    xs = f[range_names[xi]]
    cis = CartesianIndices((ncols, nrows))
    axs = []
    scs = []
    for ((obs_name, obs_label, func), ci) in zip(observables, cis)
        ax = Axis(fig[1+ci[2], ci[1]]; xscale, yscale, axis...)
        push!(axs, ax)

        data_obs = if isnothing(func)
            make_data_obs(f[obs_name])
        else
            make_data_obs(mappedarray(func, f[obs_name]))
        end

        sc = scatter!(ax, xs, data_obs; kwargs...)
        push!(scs, sc)

        ax.title = obs_label
        ax.xlabel = range_names[xi]
    end

    return FigureAxisAnything(fig, axs, scs)
end

function scansinglescatter(f, xi, observables;
    xscale=identity, yscale=identity, axis=(;), kwargs...
)
    # prep and handle inputs
    range_names = f["range_names"]
    observable_names = f["observable_names"]
    xi = _handle_rangei(range_names, xi)
    observables = [_handle_observable(observable_names, obs) for obs in observables]

    # make fig and more prep
    fig = Figure()

    # make sliders
    varis = filter(x -> x != xi, 1:length(range_names))
    slider_specs = []
    for vari in varis
        range_name = range_names[vari]
        push!(slider_specs, (label=range_name, range=f[range_name]))
    end
    sg = SliderGrid(fig[1, 1], slider_specs...)
    var_obs = [x.selected_index for x in sg.sliders]

    # this does the complex bit
    function make_data_obs(data)
        lift(var_obs...) do (varvals...)
            indices = Any[(:) for _ in 1:length(range_names)]
            for (vari, varval) in zip(varis, varvals)
                indices[vari] = varval
            end
            data[indices...]
        end
    end

    # Do the plots
    xs = f[range_names[xi]]
    ax = Axis(fig[2, 1]; xscale, yscale, axis...)
    ax.xlabel = range_names[xi]
    scs = []
    for (obs_name, obs_label, func) in observables
        data_obs = if isnothing(func)
            make_data_obs(f[obs_name])
        else
            make_data_obs(mappedarray(func, f[obs_name]))
        end

        sc = scatter!(ax, xs, data_obs; label=obs_name, kwargs...)
        push!(scs, sc)
    end
    axislegend(ax)

    return FigureAxisAnything(fig, ax, scs)
end

################################################################################
# Specific plots for cluster data
################################################################################
function br1p1()
    f = jldopen("cluster_env/runs/m3_bigscan1/bigrun1_N4.jld2")
    scanheatmaps(f, "cPs", "cRs", ["noligs", "allligs", "notense", "onboundary"]; yscale=log10, xscale=log10)
end

# FIX: There may be something wrong with either this data or plots
obs_rtlmat_bl = ("rtl_mats", x -> x[1, 1], "rtl_bl")
obs_rtlmat_tl = ("rtl_mats", x -> x[1, 5], "rtl_tl")
obs_rtlmat_tr = ("rtl_mats", x -> x[5, 5], "rtl_tr")
obs_rtlmat_br = ("rtl_mats", x -> x[5, 1], "rtl_br")
function br2p1()
    f = jldopen("cluster_env/runs/m3_bigscan2/bigrun2_N4.jld2")
    fap = scanheatmaps(f, "cPs", "cRs", ["rtl_means"]; yscale=log10, xscale=log10, colorrange=(0.0, 0.55), colormap=:Blues)

    fap
end
function br2p2(; kwargs...)
    f = jldopen("cluster_env/runs/m3_bigscan2/bigrun2_N4.jld2")
    fap = scanheatmaps(f, "cPs", "cRs", [obs_rtlmat_tl, obs_rtlmat_tr, obs_rtlmat_bl, obs_rtlmat_br]; yscale=log10, xscale=log10, colorrange=(-6.0, 6.0), colormap=:RdBu, colorscale=Makie.Symlog10(1e-9), kwargs...)

    fap
end

obs_rtlssmat_bl = ("rtlss_mats", x -> x[1, 1], "rtlss_bl")
obs_rtlssmat_tl = ("rtlss_mats", x -> x[1, 5], "rtlss_tl")
obs_rtlssmat_tr = ("rtlss_mats", x -> x[5, 5], "rtlss_tr")
obs_rtlssmat_br = ("rtlss_mats", x -> x[5, 1], "rtlss_br")
obs_rtlssmat_l = ("rtlss_mats", x -> mean(x[1, begin+1:end-1]), "rtlss_l")
obs_rtlssmat_t = ("rtlss_mats", x -> mean(x[begin+1:end-1, 5]), "rtlss_t")
obs_rtlssmat_b = ("rtlss_mats", x -> mean(x[begin+1:end-1, 1]), "rtlss_b")
obs_rtlssmat_r = ("rtlss_mats", x -> mean(x[5, begin+1:end-1]), "rtlss_r")
function br3p1()
    f = jldopen("cluster_env/runs/m3_bigscan2ss/bigrun2ss_N4.jld2")
    fap = scanheatmaps(f, "cPs", "cRs", ["rtlss_means"]; yscale=log10, xscale=log10, colormap=:Blues)

    fap
end
function br3p2(; kwargs...)
    f = jldopen("cluster_env/runs/m3_bigscan2ss/bigrun2ss_N4.jld2")
    fap = scanheatmaps(f, "cPs", "cRs", [obs_rtlssmat_tl, obs_rtlssmat_tr, obs_rtlssmat_bl, obs_rtlssmat_br]; yscale=log10, xscale=log10, colorrange=(-0.05, 0.05), colormap=:RdBu, kwargs...)
    # fap = scanheatmaps(f, "cPs", "cRs", [obs_rtlssmat_tl, obs_rtlssmat_tr, obs_rtlssmat_bl, obs_rtlssmat_br]; yscale=log10, xscale=log10, colorrange=(-100, 100), colormap=:RdBu, colorscale=Makie.Symlog10(1e-20), kwargs...)

    fap
end
function br3p3(; kwargs...)
    f = jldopen("cluster_env/runs/m3_bigscan2ss/bigrun2ss_N4.jld2")
    fap = scanheatmaps(f, "cPs", "cRs", [obs_rtlssmat_tl, obs_rtlssmat_t, obs_rtlssmat_tr, obs_rtlssmat_l, "rtlss_means", obs_rtlssmat_r, obs_rtlssmat_bl, obs_rtlssmat_b, obs_rtlssmat_br]; yscale=log10, xscale=log10, colorrange=(-0.05, 0.05), colormap=:RdBu, kwargs...)
    # fap = scanheatmaps(f, "cPs", "cRs", [obs_rtlssmat_tl, obs_rtlssmat_tr, obs_rtlssmat_bl, obs_rtlssmat_br]; yscale=log10, xscale=log10, colorrange=(-100, 100), colormap=:RdBu, colorscale=Makie.Symlog10(1e-20), kwargs...)

    fap
end

function ms1p1()
    f = jldopen("cluster_env/runs/m3_macroscan/macroscan_N4.jld2")
    scansinglescatter(f, "cRs", ["numlig", "numtense", "numboundaries"]; xscale=log10, axis=(; limits=(nothing, (0, 3))))
end
