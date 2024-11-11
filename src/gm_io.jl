struct FigureAxisAnything
    figure::Figure
    axis::Any
    obj::Any
end
display(faa::FigureAxisAnything) = display(faa.figure)

function make_linalpha_cmap(cmap; amin=0.0, amax=1.0)
    cmap = Makie.to_colormap(cmap)
    acmap = [coloralpha(color(c), a) for (c, a) in zip(cmap, LinRange(amin, amax, length(cmap)))]
    cgrad(acmap)
end

function p_get_adjmat_symsafe(graph::AbstractSimpleWeightedGraph{T,F}) where {T,F}
    if F == Num
        map(x -> if iszero(x)
                0.0
            else
                1.0
            end, adjacency_matrix(graph))
    else
        adjacency_matrix(graph)
    end
end
p_get_adjmat_symsafe(gm::AbstractGraphModel) = p_get_adjmat_symsafe(graph(gm))

function makeprefunclayout(base_layout::NetworkLayout.AbstractLayout, f, args...; kwargs...)
    copyf = copyand(f, args...; kwargs...)
    function (ca)
        base_layout(p_get_adjmat_symsafe(copyf(ca)))
    end
end

function p_make_ca_ax(dim, place)
    if dim == 2
        Axis(place)
    elseif dim == 3
        LScene(place)
    else
        throw(ErrorException("the layout dimension was not 2 or 3"))
    end
end

