################################################################################
# Getting the transition matrices in a general way
################################################################################
"""
The standard Markov chain/Master equation transition matrix where the
ij th element corresponds to a transition/rate from j to i.
"""
function transmat(graph::SimpleWeightedDiGraph; mat=false)
    tm = transpose(adjacency_matrix(graph))
    mat ? Matrix(tm) : sparse(tm)
end
function transmat(gm::AbstractGraphModel; mat=false)
    transmat(graph(gm); mat)
end
export transmat

"""
Makes the modified transition matrix that can be used to get the time
derivative as matrix multiplication only. See latex_notes/MasterEq this
is W. Follows the same convention as `transmat`.
"""
function etransmat!(tm::AbstractMatrix)
    n = size(tm)[1]
    for i in 1:n
        tm[i, i] -= sum(tm[k, i] for k in 1:n)
    end
end
@doc (@doc etransmat!)
function etransmat(tm::AbstractMatrix)
    ctm = copy(tm)
    etransmat!(ctm)
    ctm
end
function etransmat(graph::SimpleWeightedDiGraph, args...; kwargs...)
    tm = transmat(graph, args...; kwargs...)
    etransmat!(tm)
    tm
end
function etransmat(gm::AbstractGraphModel, args...; kwargs...)
    etransmat(graph(gm), args...; kwargs...)
end
export etransmat!, etransmat

"""Makes a general symbolic version of the etransmat"""
function etransmat_safe(gm::AbstractGraphModel, args...; kwargs...)
    tm = transmat(gm, args...; kwargs...)
    setm, rt = safemat_general(tm)
    etransmat!(setm)
    setm, rt
end
export etransmat_safe

################################################################################
# Dealing with etransmat and its eigensystem
################################################################################
function steadystates(gm::AbstractGraphModel{<:AbstractFloat};
    threshold=1e-8,
    checks=true,
    checkimag=checks,
    checkss=checks,
    checkothers=checks,
    returnothers::Val{Returnothers}=Val(false),
    sortothers=true
) where {Returnothers}
    esys = eigen!(etransmat(gm; mat=true))
    n = length(esys.values)

    steadystates = Vector{Vector{Float64}}(undef, 0)
    if Returnothers
        nonss = Eigensystem(numstates(gm); numtype=Val(ComplexF64))
    end
    for i in 1:n
        eval = esys.values[i]
        evec = (@view esys.vectors[:, i])

        if abs(eval) < threshold
            ss = evec / sum(evec)
            if checkimag
                maximag = maximum(x -> abs(imag(x)), ss)
                if maximag > threshold
                    @warn f"when aligning steady state eigenvector encountered imag part of {maximag} which is over {threshold}"
                end
            end
            rss = real(ss)
            if checkss
                firstviolation = findfirst(x -> (x < 0.0) || (x > 1.0), rss)
                if !isnothing(firstviolation)
                    @warn f"getting steady states with a component of {rss[firstviolation]} which is not between 0 and 1"
                end
            end
            push!(steadystates, rss)
        else
            if checkothers
                if abs(sum(evec)) > threshold
                    @warn f"found a non-steady eigenstate the components of which do not sum to 0"
                end
            end
            if Returnothers
                push!(nonss, (eval, evec))
            end
        end
    end

    if !Returnothers
        steadystates
    else
        if sortothers
            order = sortperm(real.(nonss.evals); rev=true)
            nonss = subes(nonss, order)
        end
        steadystates, nonss
    end
end
fsteadystates(args...; kwargs...) = steadystates(args...; returnothers=Val(true), kwargs...)
export steadystates, fsteadystates

supersteadystate(args...; kwargs...) = sum(steadystates(args...; kwargs..., returnothers=Val(false)))
export supersteadystate

function nonss_analysis(gm::AbstractGraphModel{<:AbstractFloat}; kwargs...)
    steadystates_, es = steadystates(gm; kwargs..., returnothers=Val(true), sortothers=true)
end
export nonss_analysis

"""
Returns list of lists each of which is a group of indices of those eigenvalues/states
which are degenerate with respect to each other.
"""
function degenerate_groups(es::Eigensystem; dthreshold=nothing, onlynontriv=true)
    if isnothing(dthreshold)
        dthreshold = 1e-6
    end
    evals = es.evals
    graph = SimpleGraph(length(es))
    for i in 1:length(evals)
        for j in i+1:length(evals)
            if abs(evals[i] - evals[j]) < dthreshold
                add_edge!(graph, i, j)
            end
        end
    end
    if edges == []
        nothing
    else
        comp = connected_components(graph)
        if onlynontriv
            filter!(x -> length(x) > 1, comp)
        end
        comp
    end
end
export degenerate_groups

function degenerate_map_(es::Eigensystem, groups)
    map = Vector{Union{Nothing,Int}}(undef, length(es))
    map .= nothing
    for (group_i, group) in enumerate(groups)
        for i in group
            map[i] = group_i
        end
    end
    map
end
function degenerate_map(es::Eigensystem, groups=nothing; returngroups=false, kwargs...)
    if isnothing(groups)
        groups = degenerate_groups(es; kwargs...)
    end
    map = degenerate_map_(es, groups)
    if !returngroups
        map
    else
        groups, map
    end
end
export degenerate_map

"""
Calculates the current matrix (antisymmetric by definition) from the transition
matrix (works with either transmat or etransmat, as diagonal elements cancel)
in the physics convention where J_ij corresponds to the net probability current
flowing from j to i.
"""
function calc_currents(etm::AbstractMatrix, ss::AbstractVector)
    cmat = similar(etm)
    for i in axes(etm, 1)
        for j in axes(etm, 2)
            cmat[i, j] = etm[i, j] * ss[j] - etm[j, i]ss[i]
        end
    end
    cmat
end
export calc_currents

function make_current_graph(gm::AbstractGraphModel{F},
    state::AbstractVector=supersteadystate(gm);
    zerothreshold=eps(F)
) where {F<:AbstractFloat}
    cmat = calc_currents(transmat(gm), state)
    cadjmat = spzeros(size(cmat))
    for i in axes(cadjmat, 1)
        for j in axes(cadjmat, 2)
            x = cmat[j, i]
            if (x > 0) && !isapprox(x, 0.0; atol=zerothreshold)
                cadjmat[i, j] = x
            end
        end
    end
    SimpleWeightedDiGraph(cadjmat)
end
export make_current_graph

################################################################################
# Stuff that changes the graph like editing/filtering edges
################################################################################
function groupgraph(graph::SimpleWeightedDiGraph, group_vect::AbstractVector;
    groups=unique(group_vect)
)
    if length(group_vect) != nv(graph)
        throw(ArgumentError("the passed group information does not have the correct length"))
    end
    if !allunique(groups)
        throw(ArgumentError("the passed `groups` list is not valid as it has repeats"))
    end
    itogroup = [findfirst(x -> x == i, groups) for i in group_vect]

    grouped_graph = SimpleWeightedDiGraph(length(groups))
    for e in edges(graph)
        gsrc = itogroup[src(e)]
        gdst = itogroup[dst(e)]
        if !isnothing(gsrc) && !isnothing(gdst)
            inc_edge!(grouped_graph, itogroup[src(e)], itogroup[dst(e)], weight(e))
        end
    end

    grouped_graph, groups
end
groupgraph(gm::AbstractGraphModel, ga; kwargs...) = groupgraph(graph(gm), ga; kwargs...)
groupgraph(gm::AbstractGraphModel, f::Function; kwargs...) = groupgraph(gm, f.(allstates(gm)); kwargs...)
export groupgraph

function groupsum(v::AbstractVector, group_vect::AbstractVector;
    groups=unique(group_vect)
)
    if length(group_vect) != length(v)
        throw(ArgumentError("the passed group information does not have the correct length"))
    end
    if !allunique(groups)
        throw(ArgumentError("the passed `groups` list is not valid as it has repeats"))
    end
    itogroup = [findfirst(x -> x == i, groups) for i in group_vect]

    grouped_v = fill(zero(eltype(v)), length(groups))
    for (i, x) in enumerate(v)
        grouped_v[itogroup[i]] += x
    end

    grouped_v, groups
end
export groupsum

function filter_edges!(graph::AbstractFloatGraph, threshold)
    for e in edges(graph)
        if weight(e) < threshold
            rem_edge!(graph, e)
        end
    end
end
function filter_edges!(gm::AbstractGraphModel, args...)
    filter_edges!(graph(gm), args...)
end
filter_edges(gm, args...) = copyand(filter_edges!, args...)(gm)
export filter_edges!, filter_edges

function keep_best_only!(graph::AbstractFloatGraph, threshold=1e-10)
    for vertex in 1:nv(graph)
        neightbors_ = copy(neighbors(graph, vertex))

        max_rate = 0.0
        for n in neightbors_
            nrate = get_weight(graph, vertex, n)
            if nrate > max_rate
                max_rate = nrate
            end
        end

        to_remove = []
        for n in neightbors_
            nrate = get_weight(graph, vertex, n)
            if (max_rate - nrate) > threshold
                push!(to_remove, n)
            end
        end
        for n in to_remove
            rem_edge!(graph, vertex, n)
        end
    end
end
function keep_best_only!(gm::AbstractGraphModel, args...)
    keep_best_only!(graph(gm), args...)
end
keep_best_only(gm, args...) = copyand(keep_best_only!, args...)(gm)
export keep_best_only!, keep_best_only

# TODO: FIX: This is not finished!!
function find_cycles(graph, next_vertex_func=next_vertex_choose_max)
    cycles = []
    cycle_map = Vector{Union{Nothing,Int}}(nothing, nv(graph))
    for vertex in 1:nv(graph)
        if isnothing(cycle_map[vertex])
            cycle = []
            true_cycle_start_i = findfirst(x -> x == vertex, cycle)
            while isnothing(true_cycle_start_i)
                push!(cycle, vertex)
                vertex = next_vertex_func(graph, vertex)
                true_cycle_start_i = findfirst(x -> x == vertex, cycle)
            end
            for v in cycle
                cycle_map[v] = length(cycles) + 1
            end
            cycle = cycle[true_cycle_start_i:end]
            push!(cycles, cycle)
        end
    end
    (cycles, cycle_map)
end
