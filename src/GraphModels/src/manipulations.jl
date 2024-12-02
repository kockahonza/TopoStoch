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
# Doing stuff with the etransmat
################################################################################
function steadystates(gm::AbstractGraphModel{<:AbstractFloat};
    threshold=1e-8,
    checkimag=true,
    checkothers=true,
    returnothers::Val{Returnothers}=Val(false)
) where {Returnothers}
    esys = eigen!(etransmat(gm; mat=true))
    n = length(esys.values)

    steadystates = Vector{Vector{Float64}}(undef, 0)
    if Returnothers
        otherevals = Vector{ComplexF64}(undef, 0)
        otherevecs = Vector{Vector{ComplexF64}}(undef, 0)
    end
    for i in 1:n
        eval = esys.values[i]
        evec = (@view esys.vectors[:, i])

        if abs(eval) < threshold
            ss = evec / sum(evec)
            if checkimag
                maximag = maximum(x -> abs(imag(x)), ss)
                if maximag > threshold
                    @warn f"when aligning ss eigenvector encountered imag part of {maximag} over {threshold}"
                end
            end
            push!(steadystates, real(ss))
        else
            if checkothers
                if abs(sum(evec)) > threshold
                    @warn f"found a non-physical eigenstate the components of which do not sum to 0"
                end
            end
            if Returnothers
                push!(otherevals, eval)
                push!(otherevecs, evec)
            end
        end
    end

    if !Returnothers
        steadystates
    else
        steadystates, (; evals=otherevals, evecs=otherevecs)
    end
end
export steadystates

supersteadystate(args...; kwargs...) = sum(steadystates(args...; kwargs..., returnothers=Val(false)))
export supersteadystate

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

function make_current_graph(gm::AbstractGraphModel, state::AbstractVector)
    cmat = calc_currents(transmat(gm), state)
    cadjmat = spzeros(size(cmat))
    for i in axes(cadjmat, 1)
        for j in axes(cadjmat, 2)
            x = cmat[j, i]
            if x > 0
                cadjmat[i, j] = x
            end
        end
    end
    SimpleWeightedDiGraph(cadjmat)
end
export make_current_graph

################################################################################
# Editing/filtering edges in concrete GraphModel graphs
################################################################################
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
