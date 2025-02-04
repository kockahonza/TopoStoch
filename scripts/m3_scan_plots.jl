using DrWatson
@quickactivate "TopoStochSim"

using JLD2
using GLMakie
using Printf

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
    if !isa(obs, Integer)
        maybe_obs = findfirst(x -> x == obs, observable_names)
        if !isnothing(maybe_obs)
            obs = maybe_obs
        else
            throw(ArgumentError(@sprintf "spec %s has unrecognised observable name" string(spec)))
        end
    end
    obs = observable_names[obs]
    obs
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

function scanheatmaps(f, xi, yi, observables; xscale=identity, yscale=identity)
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
    (nrows, ncols) = ntogridsize1(nplots)

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
    for (obs_name, ci) in zip(observables, cis)
        ax = Axis(fig[1+ci[2], ci[1]*2-1]; xscale, yscale)
        push!(axs, ax)
        hm = heatmap!(ax, xs, ys, make_data_obs(f[obs_name]))
        Colorbar(fig[1+ci[2], ci[1]*2], hm)
        ax.title = obs_name
        ax.xlabel = range_names[xi]
        ax.ylabel = range_names[yi]
    end

    return FigureAxisAnything(fig, axs, nothing)
end
