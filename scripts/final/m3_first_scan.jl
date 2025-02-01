using DrWatson
@quickactivate "TopoStochSim"

using GLMakie
using JLD2

using ComplexAllostery
using ComplexAllostery.Model3

################################################################################
# First scan test, keeps rs at 1 and changes eb, cP, cR
################################################################################
function scan1()
    ca = makeCAGM(4, vars_simplified())
    cca_fact = make_factory(
        ca,
        [symvars.sim_rs; symvars.energies[3];
            symvars.concentrations[1]; symvars.redconcentration]
    )

    noligs_is = Int[]
    allligs_is = Int[]
    notense_is = Int[]
    alltense_is = Int[]
    for (i, st) in enumerate(allstates(ca))
        if all(x -> x == 0, st.occupations)
            push!(noligs_is, i)
        end
        if all(x -> x == 1, st.occupations)
            push!(allligs_is, i)
        end
        if all(x -> x != 1, st.conformations)
            push!(notense_is, i)
        end
        if all(x -> x == 1, st.conformations)
            push!(alltense_is, i)
        end
    end
    onboundary_is = unique(vcat(noligs_is, allligs_is, notense_is, alltense_is))

    eb_range = LinRange(0, 5, 10)
    cP_range = LinRange(0, 1, 20)
    cR_range = LinRange(1, 20, 20)

    noligs_probs = Array{Float64,3}(undef, length(eb_range), length(cP_range), length(cR_range))
    allligs_probs = Array{Float64,3}(undef, length(eb_range), length(cP_range), length(cR_range))
    notense_probs = Array{Float64,3}(undef, length(eb_range), length(cP_range), length(cR_range))
    alltense_probs = Array{Float64,3}(undef, length(eb_range), length(cP_range), length(cR_range))
    boundary_probs = Array{Float64,3}(undef, length(eb_range), length(cP_range), length(cR_range))

    for (eb_i, eb) in enumerate(eb_range)
        for (cP_i, cP) in enumerate(cP_range)
            for (cR_i, cR) in enumerate(cR_range)
                cca = cca_fact(1, 1, 1, 0, eb, cP, cR)
                ss = supersteadystate(cca)
                noligs_probs[eb_i, cP_i, cR_i] = sum(ss[noligs_is])
                allligs_probs[eb_i, cP_i, cR_i] = sum(ss[allligs_is])
                notense_probs[eb_i, cP_i, cR_i] = sum(ss[notense_is])
                alltense_probs[eb_i, cP_i, cR_i] = sum(ss[alltense_is])
                boundary_probs[eb_i, cP_i, cR_i] = sum(ss[onboundary_is])
            end
        end
    end

    boundary_probs

    jldsave(datadir("scans/m3_scan1_1.jld2");
        eb_range, cP_range, cR_range,
        noligs_probs, allligs_probs, notense_probs, alltense_probs,
        boundary_probs
    )
end

function viz_scan1(xi, yi, f=jldopen(datadir("scans/m3_scan1_1.jld2")))
    fig = Figure()

    ranges = [f["eb_range"], f["cP_range"], f["cR_range"]]
    var_labels = ["eb", "cP", "cR"]

    varis = filter(x -> !(x in (xi, yi)), 1:length(ranges))
    slider_specs = []
    for vari in varis
        label = if !isnothing(var_labels)
            var_labels[vari]
        else
            string(vari)
        end
        push!(slider_specs, (label=label, range=ranges[vari]))
    end
    sg = SliderGrid(fig[1, 1:4], slider_specs...)
    var_obs = [x.selected_index for x in sg.sliders]

    function make_data_obs(data)
        lift(var_obs...) do (varvals...)
            indices = Any[(:) for _ in 1:length(ranges)]
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

    noligs_data = make_data_obs(f["noligs_probs"])
    allligs_data = make_data_obs(f["allligs_probs"])
    notense_data = make_data_obs(f["notense_probs"])
    alltense_data = make_data_obs(f["alltense_probs"])
    data = [noligs_data, allligs_data, notense_data, alltense_data]
    titles = ["noligs", "allligs", "notense", "alltense"]

    for (i, ci) in enumerate(CartesianIndices((1:2, 1:2)))
        ax = Axis(fig[1+ci[1], ci[2]*2-1])
        hm = heatmap!(ax, ranges[xi], ranges[yi], data[i])
        ax.title = titles[i]
        ax.xlabel = var_labels[xi]
        ax.ylabel = var_labels[yi]
        Colorbar(fig[1+ci[1], ci[2]*2], hm)
    end

    fig
end

################################################################################
"""
A somewhat general plotting thing similar to the above
"""
function scanheatmap(data, ranges, xi, yi; var_labels=nothing)
    fig = Figure()

    varis = filter(x -> !(x in (xi, yi)), 1:length(ranges))
    slider_specs = []
    for vari in varis
        label = if !isnothing(var_labels)
            var_labels[vari]
        else
            string(vari)
        end
        push!(slider_specs, (label=label, range=ranges[vari]))
    end
    sg = SliderGrid(fig[1, 1], slider_specs...)
    var_obs = [x.selected_index for x in sg.sliders]

    heatmapdata = lift(var_obs...) do (varvals...)
        indices = Any[(:) for _ in 1:length(ranges)]
        for (vari, varval) in zip(varis, varvals)
            indices[vari] = varval
        end
        if xi < yi
            data[indices...]
        else
            transpose(data[indices...])
        end
    end

    ax = Axis(fig[2, 1])
    hm = heatmap!(ax, ranges[xi], ranges[yi], heatmapdata)
    Colorbar(fig[2, 2], hm)

    fig
end
