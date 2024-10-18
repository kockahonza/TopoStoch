using DrWatson
@quickactivate "TopoStochSim"

using Graphs, MetaGraphsNext, SimpleWeightedGraphs, GraphIO
using GLMakie, GraphMakie, NetworkLayout
using StatsBase
using StaticArrays, NamedArrays
using PyFormattedStrings

const gc_map = NamedArray{String}(4, 4, 4)
setnames!(gc_map, ["T", "C", "A", "G"], 1)
setnames!(gc_map, ["T", "C", "A", "G"], 2)
setnames!(gc_map, ["T", "C", "A", "G"], 3)

gc_map["G", "G", "G"] = "G"
gc_map["G", "G", "A"] = "G"
gc_map["G", "G", "C"] = "G"
gc_map["G", "G", "T"] = "G"
gc_map["G", "A", "G"] = "E"
gc_map["G", "A", "A"] = "E"
gc_map["G", "A", "C"] = "D"
gc_map["G", "A", "T"] = "D"
gc_map["G", "C", "G"] = "A"
gc_map["G", "C", "A"] = "A"
gc_map["G", "C", "C"] = "A"
gc_map["G", "C", "T"] = "A"
gc_map["G", "T", "G"] = "V"
gc_map["G", "T", "A"] = "V"
gc_map["G", "T", "C"] = "V"
gc_map["G", "T", "T"] = "V"
gc_map["A", "G", "G"] = "R"
gc_map["A", "G", "A"] = "R"
gc_map["A", "G", "C"] = "S"
gc_map["A", "G", "T"] = "S"
gc_map["A", "A", "G"] = "K"
gc_map["A", "A", "A"] = "K"
gc_map["A", "A", "C"] = "N"
gc_map["A", "A", "T"] = "N"
gc_map["A", "C", "G"] = "T"
gc_map["A", "C", "A"] = "T"
gc_map["A", "C", "C"] = "T"
gc_map["A", "C", "T"] = "T"
gc_map["A", "T", "G"] = "M"
gc_map["A", "T", "A"] = "I"
gc_map["A", "T", "C"] = "I"
gc_map["A", "T", "T"] = "I"
gc_map["C", "G", "G"] = "R"
gc_map["C", "G", "A"] = "R"
gc_map["C", "G", "C"] = "R"
gc_map["C", "G", "T"] = "R"
gc_map["C", "A", "G"] = "Q"
gc_map["C", "A", "A"] = "Q"
gc_map["C", "A", "C"] = "H"
gc_map["C", "A", "T"] = "H"
gc_map["C", "C", "G"] = "P"
gc_map["C", "C", "A"] = "P"
gc_map["C", "C", "C"] = "P"
gc_map["C", "C", "T"] = "P"
gc_map["C", "T", "G"] = "L"
gc_map["C", "T", "A"] = "L"
gc_map["C", "T", "C"] = "L"
gc_map["C", "T", "T"] = "L"
gc_map["T", "G", "G"] = "W"
gc_map["T", "G", "A"] = "*"
gc_map["T", "G", "C"] = "C"
gc_map["T", "G", "T"] = "C"
gc_map["T", "A", "G"] = "*"
gc_map["T", "A", "A"] = "*"
gc_map["T", "A", "C"] = "Y"
gc_map["T", "A", "T"] = "Y"
gc_map["T", "C", "G"] = "S"
gc_map["T", "C", "A"] = "S"
gc_map["T", "C", "C"] = "S"
gc_map["T", "C", "T"] = "S"
gc_map["T", "T", "G"] = "L"
gc_map["T", "T", "A"] = "L"
gc_map["T", "T", "C"] = "F"
gc_map["T", "T", "T"] = "F"

const bases = ["G", "A", "C", "T"]
const codes = unique(gc_map)
codetoi(code) = findfirst(x -> x == code, codes)
dnatocode(dna...) = gc_map[dna...]
dnatoi(dna...) = codetoi(dnatocode(dna...))

function makegraph()
    graph = SimpleWeightedDiGraph(length(codes))
    self_mutations = fill(0, length(codes))

    for d1 in bases
        for d2 in bases
            for d3 in bases
                mutations = []
                for m in bases
                    if m != d1
                        push!(mutations, (m, d2, d3))
                    end
                    if m != d2
                        push!(mutations, (d1, m, d3))
                    end
                    if m != d3
                        push!(mutations, (d1, d2, m))
                    end
                end
                basecode = dnatocode(d1, d2, d3)
                for mutation in mutations
                    mutationcode = dnatocode(mutation...)
                    if basecode == mutationcode
                        self_mutations[codetoi(basecode)] += 1
                    else
                        basei = codetoi(basecode)
                        mutationi = codetoi(mutationcode)
                        add_edge!(graph, basei, mutationi, get_weight(graph, basei, mutationi) + 1.0)
                    end
                end
            end
        end
    end

    graph, self_mutations
end

function plot(graph, self_mutations, args...; enable_drag=false, hide_below=nothing, kwargs...)
    if !isnothing(hide_below)
        graph = copy(graph)
        for e in edges(graph)
            if e.weight < hide_below
                rem_edge!(graph, e)
            end
        end
    end

    edge_color = [e.weight for e in edges(graph)]

    nodelabels = [f"{code}({sm})" for (code, sm) in zip(codes, self_mutations)]

    f, a, p = graphplot(graph, args...; nlabels=nodelabels, node_size=15, edge_color, kwargs...)
    # Colorbar(f[1, 2], p.plots[1].plots[1])

    if enable_drag
        deregister_interaction!(a, :rectanglezoom)
        function node_drag_action(state, idx, event, axis)
            p[:node_pos][][idx] = event.data
            p[:node_pos][] = p[:node_pos][]
        end
        ndrag = NodeDragHandler(node_drag_action)
        register_interaction!(a, :ndrag, ndrag)
    end

    Makie.FigureAxisPlot(f, a, p)
end

function main()
    g, sm = makegraph()
    plot(g, sm; enable_drag=true, layout=Shell(), hide_below=1.5, curve_distance_usage=false)
end

main()
