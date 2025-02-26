using DrWatson
@quickactivate "TopoStochSim"

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

function makeruledf(N=2:6)
    if !isa(N, AbstractVector)
        N = [N]
    end
    df = DataFrame(code=Int[], N=Int[], numscc=Int[], numac=Int[], numss=Int[])
    for code in 0:255
        for n in N
            ned = make_ce_ned(n, code; show=false)
            numscc = length(strongly_connected_components(graph(ned)))
            numac = length(attracting_components(graph(ned)))
            numss = length(steadystates(ned))
            push!(df, (code, n, numscc, numac, numss))
        end
    end
    df
end

function plotvsN_forallcodes(var, df=makeruledf(); colormap=:viridis, roffset=0.1)
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

function makeNdf(code, Ns=2:10)
    df = DataFrame(N=Int[], numscc=Int[], numac=Int[], numss=Int[])
    for N in Ns
        println(N)
        ned = make_ce_ned(N, code; show=false)
        numscc = length(strongly_connected_components(graph(ned)))
        numac = length(attracting_components(graph(ned)))
        numss = length(steadystates(ned))
        push!(df, (N, numscc, numac, numss))
    end
    df
end
