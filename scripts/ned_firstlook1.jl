using DrWatson
@quickactivate "TopoStochSim"

using Printf
using DataFrames
using GLMakie

using GraphModels
using NonEqDigits

function gsummary1(g)
    # num of completely isolated nodes
    isolated = []
    # num of nodes with only incoming edges
    deadends = []
    # num of nodes with only outgoing edges
    neverback = []
    for ni in 1:nv(g)
        if isempty(neighbors(g, ni))
            push!(isolated, ni)
        else
            if isempty(outneighbors(g, ni))
                push!(deadends, ni)
            end
            if isempty(inneighbors(g, ni))
                push!(neverback, ni)
            end
        end
    end
    isolated, deadends, neverback
end
gsummary1(gm::AbstractGraphModel) = gsummary1(graph(gm))

function makedfbycodes(N=2:6, codes=ca_ucodes_bydeg();
    add_classes=true,
)
    if !isa(N, AbstractVector)
        N = [N]
    end

    df = DataFrame(
        code_i=Int[],
        code=Int[],
        N=Int[],
        numscc=Int[],
        numac=Int[],
        numnontrivac=Int[],
        numss=Int[]
    )
    for (code_i, code) in enumerate(codes)
        for n in N
            ned = make_ca_ned(n, code; show=false)
            numscc = length(strongly_connected_components(graph(ned)))
            acs = attracting_components(graph(ned))
            numac = length(acs)
            numnontrivac = count(x -> length(x) > 1, acs)
            numss = length(steadystates(ned))
            push!(df, (
                code_i,
                code,
                n,
                numscc,
                numac,
                numnontrivac,
                numss
            ))
        end
    end

    if add_classes
        cl_df = ca_load_classes_df()
        df.class = cl_df.class[df.code.+1]
        df.subclass_i = cl_df.subclass_i[df.code.+1]
        df.csc_i = cl_df.csc_i[df.code.+1]
    end

    df
end

function plotvsN_forcodes(var, df; colormap=:managua, roffset=0.1)
    fig = Figure()
    ax = Axis(fig[1, 1])
    cmap = cgrad(colormap)
    num_codes = maximum(df.code_i)
    for code in unique(df.code)
        subdf = df[df.code.==code, :]
        Ns = subdf[!, :N]
        vars = subdf[!, var]
        varoffsets = roffset .* rand()
        lines!(ax, Ns, vars .+ varoffsets; color=get(cmap, subdf[1, :code_i] / num_codes))
    end
    Colorbar(fig[1, 2]; colorrange=(0.0, num_codes), colormap=cmap)
    fig
end

function plotvsN_forallcodes(var, df; colormap=:managua, roffset=0.1)
    fig = Figure()
    ax = Axis(fig[1, 1])
    cmap = cgrad(colormap)
    for code in 0:255
        subdf = df[df.code.==code, :]
        Ns = subdf[!, :N]
        vars = subdf[!, var]
        Noffsets = roffset .* rand(length(Ns))
        varoffsets = roffset .* rand(length(Ns))
        lines!(ax, Ns .+ Noffsets, vars .+ varoffsets; color=get(cmap, code / 255.0))
    end
    Colorbar(fig[1, 2]; colorrange=(0.0, 255.0), colormap=cmap)
    fig
end

function ourclassesCN(N=2:6, codes=ca_ucodes_bydeg();
    add_ca_classes=false,
)
    if !isa(N, AbstractVector)
        N = [N]
    end

    df = DataFrame(
        # System specs
        code_i=Int[],
        code=Int[],
        N=Int[],
        # General properties
        numarrows=Int[],
        numscc=Int[],
        numac=Int[],
        numnontrivac=Int[],
        numss=Int[],
        onlynontrivaclen=Int[],
        # Custom classification properties
        class=Int[]
    )
    for (code_i, code) in enumerate(codes)
        for n in N
            ned = make_ca_ned(n, code; show=false)

            # Base stuff
            sccs = strongly_connected_components(graph(ned))
            acs = attracting_components(graph(ned))

            numarrows = ca_calc_numarrows(code)
            numscc = length(sccs)
            numac = length(acs)
            numnontrivac = count(x -> length(x) > 1, acs)
            numss = length(steadystates(ned))

            # classification
            cond = condensation(graph(ned), sccs)
            sing, nonsing = findsingularsccs(ned, sccs, cond)

            if numscc == 1 # fully connected
                class = 1
            elseif length(nonsing) == 0 # we have only singular points
                class = -1
            elseif length(nonsing) == 1 # almost fully connected
                class = 2
            elseif aresimplearrow(cond, nonsing) # single arrow + singular sccs
                class = 10
            else
                class = 0
            end

            push!(df, (
                code_i,
                code,
                n,
                numarrows,
                numscc,
                numac,
                numnontrivac,
                numss,
                get_onlynontrivaclen(ned),
                class
            ))
        end
    end

    if add_ca_classes
        cl_df = ca_load_classes_df()
        df.ca_class = cl_df.class[df.code.+1]
        df.ca_subclass_i = cl_df.subclass_i[df.code.+1]
        df.ca_csc_i = cl_df.csc_i[df.code.+1]
    end

    df
end

function plotdf!(ax, df, x, y, c=:code;
    gbi=:code_i,
    do3d=false,
    colormap=:managua,
    setlabels=true
)
    c_func = if isnothing(c)
        if do3d
            throw(ArgumentError("cannot use do3d without passing a color"))
        end

        colorinfo = nothing

        function (sdf)
            Cycled(1)
        end
    else
        cs = df[!, c]
        colorrange = (minimum(cs), maximum(cs))
        colorrange_span = colorrange[2] - colorrange[1]
        colorrange_base = colorrange[1]
        cmap = cgrad(colormap)

        colorinfo = (; colorrange, colormap)

        if do3d
            function (sdf)
                xx = sdf[!, c]
                xx, get(cmap, (xx .- colorrange_base) ./ colorrange_span)
            end
        else
            function (sdf)
                get(cmap, (sdf[!, c] .- colorrange_base) ./ colorrange_span)
            end
        end
    end

    gdf = groupby(df, gbi)
    for sdf in gdf
        numlines = size(sdf)[1]
        if do3d
            dim3vals, cvals = c_func(sdf)
            if numlines != 1
                lines!(ax, dim3vals, sdf[!, x], sdf[!, y]; color=cvals)
            else
                scatter!(ax, dim3vals, sdf[!, x], sdf[!, y]; color=cvals)
            end
        else
            if numlines != 1
                lines!(ax, sdf[!, x], sdf[!, y]; color=c_func(sdf))
            else
                scatter!(ax, sdf[!, x], sdf[!, y]; color=c_func(sdf))
            end
        end
    end

    if setlabels
        if do3d
            ax.xlabel = string(c)
            ax.ylabel = string(x)
            ax.zlabel = string(y)
        else
            ax.xlabel = string(x)
            ax.ylabel = string(y)
        end
    end

    colorinfo
end
function plotdf(args...; do3d=false, kwargs...)
    fig = Figure()
    if do3d
        ax = Axis3(fig[1, 1])
    else
        ax = Axis(fig[1, 1])
    end
    maybeci = plotdf!(ax, args...; do3d, kwargs...)

    if !isnothing(maybeci)
        Colorbar(fig[1, 2]; maybeci...)
    end

    FigureAxisAnything(fig, ax, nothing)
end

function makecondplot(code, N; kwargs...)
    ned = make_ca_ned(N, code)
    cond = condensation(graph(ned))
    graphplot(cond; kwargs...)
end
