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

function makeruledf(N=3)
    if !isa(N, AbstractVector)
        N = [N]
    end
    df = DataFrame(numscc=Int[], numac=Int[], numss=Int[])
    for code in 1:255
        for n in N
            ned = make_ce_ned(n, code)
            numscc = length(strongly_connected_components(graph(ned)))
            numac = length(attracting_components(graph(ned)))
            numss = length(steadystates2(ned))
            push!(df, (numscc, numac, numss))
        end
    end
    df
end
