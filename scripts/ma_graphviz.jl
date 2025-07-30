using DrWatson;
@quickactivate "TopoStochSim";

using GraphModels
using NonEqDigits
using GraphvizDotLang
using GraphvizDotLang: digraph, subgraph, node, edge, attr
using ColorSchemes

using PythonCall
pd = pyimport("pydot")

function plot_0s_v_1s_gv(ned::NonEqDigitsGM{Loop,2,L};
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
        for ci in 1:(length(cluster_is)-1)
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

function molaut_gv(ma::MolAut;
    cluster=false,
    force_layout=cluster,
    simple_nodes=false,
    dim_nonacs=nothing,
    highlight_acs=nothing,
    edge_colormap=nothing,
    node_colors=nothing,
    kwargs...
)
    mg = ma.mg

    g = digraph(;
        rankdir="LR",
        # ranksep="1",
        kwargs...
    )

    if simple_nodes
        g |> attr(:node;
            shape="point",
            fixedsize="true",
            width="0.1",
        )
    end

    # prep node styling/coloring etc
    if highlight_acs == true
        highlight_acs = "crimson"
    end
    if dim_nonacs == true
        dim_nonacs = "gray"
    end
    if !isnothing(highlight_acs) || !isnothing(dim_nonacs)
        acs = attracting_components(mg)
        vs_in_acs = reduce(vcat, acs)
    end
    function node_style(l)
        rslt = Dict{Symbol,Any}()

        if !isnothing(node_colors)
            rslt[:fillcolor] = "#" * hex(node_colors[l])
            rslt[:style] = "filled"
            if !isnothing(highlight_acs)
                if code_for(mg, l) in vs_in_acs
                    rslt[:color] = highlight_acs
                    rslt[:penwidth] = "4"
                end
            end
        else
            if !isnothing(highlight_acs)
                if code_for(mg, l) in vs_in_acs
                    rslt[:color] = highlight_acs
                end
            end
        end

        if !isnothing(dim_nonacs)
            if !(code_for(mg, l) in vs_in_acs)
                rslt[:color] = dim_nonacs
            end
        end

        rslt
    end

    # add all nodes to graph
    if cluster
        clusters = [subgraph(g, "cluster_$(i)_ones") for i in 0:ma.L]
        cluster_labels = [[] for _ in 1:length(clusters)]

        for l in labels(mg)
            num1s = count(x -> x == 1, l)
            c = clusters[1+num1s]
            c |> node(mg[l];
                node_style(l)...
            )
            push!(cluster_labels[1+num1s], l)
        end

        if force_layout
            for ci in 1:(length(cluster_labels)-1)
                for s in cluster_labels[ci]
                    for d in cluster_labels[ci+1]
                        g |> edge(mg[s], mg[d];
                            style="invis",
                            weight="100"
                        )
                    end
                end
            end
        end
    else
        for l in labels(mg)
            g |> node(mg[l];
                node_style(l)...
            )
        end
    end

    if !isnothing(edge_colormap)
        if edge_colormap isa Symbol
            edge_colormap = colorschemes[edge_colormap]
        end
        edge_vals = unique(Graphs.weights(mg))
        mn, mx = extrema(edge_vals)
        if isapprox(mn, mx)
            @warn "skipping coloring as all edged have the same value"
            cfunc = nothing
        else
            cfunc = cval -> "#" * hex(get(edge_colormap, (cval - mn) / (mx - mn)))
        end
    else
        cfunc = nothing
    end

    # prep edge styling
    function edge_style(s, d)
        rslt = Dict{Symbol,Any}()

        if !isnothing(cfunc)
            rslt[:color] = cfunc(mg[s, d])
        end

        if !isnothing(highlight_acs)
            if (code_for(mg, s) in vs_in_acs) && (code_for(mg, d) in vs_in_acs)
                rslt[:color] = highlight_acs
            end
        end
        if !isnothing(dim_nonacs)
            if !(code_for(mg, s) in vs_in_acs) || !(code_for(mg, d) in vs_in_acs)
                rslt[:color] = dim_nonacs
            end
        end

        rslt
    end

    # add edges to graph
    for (s, d) in edge_labels(mg)
        g |> edge(mg[s], mg[d];
            edge_style(s, d)...
        )
    end

    g
end

function parse_dot_positions(fname; size=nothing)
    g = pd.graph_from_dot_file(fname)[0]

    # make a list of all nodes that we later filter through
    nodes = collect(g.get_nodes())
    for sg in g.get_subgraphs()
        for n in sg.get_nodes()
            push!(nodes, n)
        end
    end

    rslt = Dict()

    i = 0
    for n in nodes
        name = pyconvert(String, n.get_name())
        if occursin(r"^[01]{6}$", name)
            i += 1

            l = [x == '0' ? 0 : (x == '1' ? 1 : throw(ErrorException("could not parse name"))) for x in name]

            posstring = pyconvert(String, n.get_pos())[2:end-1]
            pos = parse.(Float64, split(posstring, ","))

            rslt[l] = pos
        end
    end

    if !isnothing(size)
        if i != size
            throw(ArgumentError("did not find the expected number of nodes"))
        end
    end

    rslt
end

"""
Output of this can be used as a layout for plotgm
"""
function get_gv_layout_positions(ma::MolAut;
    normalize_x=false,
    normalize_y=false,
    kwargs...,
)
    g = molaut_gv(ma; kwargs...)

    fname = tempname() * ".dot"
    GraphvizDotLang.save(g, fname; format="dot")

    posdict = parse_dot_positions(fname; size=2^ma.L)

    poslist = []
    for v in 1:nv(ma.mg)
        l = label_for(ma.mg, v)
        push!(poslist, posdict[l])
    end

    if normalize_x != false
        if normalize_x == true
            normalize_x = (0, ma.L)
        end

        xs = getindex.(poslist, 1)
        xmin, xmax = extrema(xs)

        normxs = normalize_x[1] .+ ((xs .- xmin) ./ (xmax - xmin)) .* (normalize_x[2] - normalize_x[1])
        for i in 1:length(poslist)
            poslist[i][1] = normxs[i]
        end
    end

    if normalize_y != false
        if normalize_y == true
            normalize_y = (0, ma.L)
        end

        ys = getindex.(poslist, 2)
        xmin, xmax = extrema(ys)

        normxs = normalize_y[1] .+ ((ys .- xmin) ./ (xmax - xmin)) .* (normalize_y[2] - normalize_y[1])
        for i in 1:length(poslist)
            poslist[i][2] = normxs[i]
        end
    end


    poslist
end
