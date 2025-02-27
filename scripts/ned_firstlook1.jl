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

function makedfbycodes(N=2:6, codes=ce_ucodes_bydeg();
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
            ned = make_ce_ned(n, code; show=false)
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
        cl_df = ce_load_classes_df()
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
