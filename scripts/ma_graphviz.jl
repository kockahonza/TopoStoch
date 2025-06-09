using NonEqDigits
using GraphvizDotLang: digraph, subgraph, node, edge

function kaka(ned::NonEqDigitsGM{Loop,2,L};
    force_order=true,
    double_arrows=true,
    kwargs...
) where {L}
    g = digraph(;
        rankdir="LR",
        ranksep="1",
        kwargs...
    )

    node_clusters = [subgraph(g, (@sprintf "cluster %s" string(i))) for i in 0:L]
    cluster_is = [[] for _ in 0:L]

    for i in 1:numstates(ned)
        st = itostate(i, ned)
        num_ones = count(x -> x == 1, st)
        cluster_index = 1 + num_ones
        node_clusters[cluster_index] |> node(string(i); label=join(st))
        push!(cluster_is[cluster_index], i)
    end
    if force_order
        for ci in 1:(L-1)
            for i1 in cluster_is[ci]
                for i2 in cluster_is[ci+1]
                    g |> edge(string(i1), string(i2);
                        style="invis",
                        # weight="100"
                    )
                end
            end
        end
    end

    if double_arrows
        ng = graph(ned)
        for i in 1:numstates(ned)
            for j in (i+1):numstates(ned)
                itoj = !iszero(get_weight(ng, i, j))
                jtoi = !iszero(get_weight(ng, j, i))
                if itoj && jtoi
                    g |> edge(string(i), string(j);
                        dir="both",
                        # weight="10"
                    )
                elseif itoj
                    g |> edge(
                        string(i),
                        string(j);
                        # weight="10"
                    )
                elseif jtoi
                    g |> edge(
                        string(j),
                        string(i);
                        # weight="10"
                    )
                end
            end
        end
    else
        for e in edges(graph(ned))
            g |> edge(
                string(e.src),
                string(e.dst);
                # weight="10"
            )
        end
    end

    g
end
