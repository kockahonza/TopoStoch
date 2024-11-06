using DrWatson
@quickactivate "TopoStochSim"

include(srcdir("gm.jl"))

using Graphs, MetaGraphsNext, SimpleWeightedGraphs
using GLMakie, GraphMakie, NetworkLayout
using StatsBase
using StaticArrays
using PyFormattedStrings

################################################################################
# Simple four state mode
################################################################################
struct SimpleFourStateGM <: AbstractStatefulGraphModel{2,Int}
    Nx::Int
    Ny::Int
    γe::Float64
    γi::Float64
    graph::SimpleWeightedDiGraph{Int,Float64}
    states::Vector{State{2,Int}}
    allowed
    function SimpleFourStateGM(Nx, Ny; γe=0.95, γi=0.05)
        states = []
        for x in 0:(Nx-1)
            for y in 0:(Ny-1)
                for substate in 1:4
                    push!(states, State(SA[x, y], substate))
                end
            end
        end

        allowed((; extstate, substate)::State{2,Int}) = (0 <= extstate[1] < Nx) && (0 <= extstate[2] < Ny) && (1 <= substate <= 4)

        graph = SimpleWeightedDiGraph(length(states))

        # Build self so that we can use statetoi
        self = new(Nx, Ny, γe, γi, graph, states, allowed)

        # Add the transitions
        for x in 0:(Nx-1)
            for y in 0:(Ny-1)
                # The favoured ones
                # A->B
                dest = State(SA[x, y+1], 2)
                if allowed(dest)
                    add_edge!(graph, statetoi(self, State(SA[x, y], 1)), statetoi(self, dest), γe)
                end
                # B->C
                dest = State(SA[x+1, y], 3)
                if allowed(dest)
                    add_edge!(graph, statetoi(self, State(SA[x, y], 2)), statetoi(self, dest), γe)
                end
                # C->D
                dest = State(SA[x, y-1], 4)
                if allowed(dest)
                    add_edge!(graph, statetoi(self, State(SA[x, y], 3)), statetoi(self, dest), γe)
                end
                # D->A
                dest = State(SA[x-1, y], 1)
                if allowed(dest)
                    add_edge!(graph, statetoi(self, State(SA[x, y], 4)), statetoi(self, dest), γe)
                end

                # The unfavoured
                # A->D
                dest = State(SA[x, y], 4)
                if allowed(dest)
                    add_edge!(graph, statetoi(self, State(SA[x, y], 1)), statetoi(self, dest), γi)
                end
                # D->C
                dest = State(SA[x, y], 3)
                if allowed(dest)
                    add_edge!(graph, statetoi(self, State(SA[x, y], 4)), statetoi(self, dest), γi)
                end
                # C->B
                dest = State(SA[x, y], 2)
                if allowed(dest)
                    add_edge!(graph, statetoi(self, State(SA[x, y], 3)), statetoi(self, dest), γi)
                end
                # B->A
                dest = State(SA[x, y], 1)
                if allowed(dest)
                    add_edge!(graph, statetoi(self, State(SA[x, y], 2)), statetoi(self, dest), γi)
                end
            end
        end

        self
    end
end
function graph(sfs::SimpleFourStateGM)
    sfs.graph
end
function states(sfs::SimpleFourStateGM)
    sfs.states
end

function sfs_subtype_to_letter(subtype)
    'A' + subtype - 1
end

function plotGM(sfs::SimpleFourStateGM, args...; ecolorbar=true, kwargs...)
    states_ = states(sfs)
    offsets = SA[Point(1, 1), Point(1, -1), Point(-1, -1), Point(-1, 1)] .* 0.15
    layout = [Point(s.extstate[1], s.extstate[2]) + offsets[s.substate] for s in states_]

    ilabels = [sfs_subtype_to_letter(s.substate) for s in states_]
    node_color = [:snow3 for _ in 1:length(states_)]
    node_size = [15 for _ in 1:length(states_)]

    edge_color = [e.weight for e in edges(graph(sfs))]
    edge_attr = arrow_attr = (; colormap=:blues)

    fap = graphplot(sfs.graph, args...; layout, ilabels, node_color, node_size, edge_color, edge_attr, arrow_attr, kwargs...)

    if ecolorbar
        edgecolorbar(fap)
    end

    fap
end

################################################################################
# 3D counting allostery model
################################################################################
struct Allostery3D1 <: AbstractStatefulGraphModel{3,Int}
    N::Int
    graph::SimpleWeightedDiGraph{Int,Float64}
    states::Vector{State{3,Int}}
    allowed
    function Allostery3D1(N; weights=@SVector [1.0 for _ in 1:8])
        if weights == :random
            weights = @SVector [rand() for _ in 1:8]
        end

        states = []
        for T in 0:N
            for OT in 0:(T)
                for OR in 0:(N-T)
                    push!(states, State(SA[T, OR, OT], 1))
                end
            end
        end
        states

        function allowed((; extstate, substate)::State{3,Int})
            T, OR, OT = extstate
            return (T <= N) && (0 <= OT <= T) && (0 <= OR <= (N - T)) && (substate == 1)
        end

        graph = SimpleWeightedDiGraph(length(states))

        # Build self so that we can use statetoi
        self = new(N, graph, states, allowed)

        # Add (de)oxygenation edges
        for i in eachindex(states)
            (T, OR, OT) = states[i].extstate
            oxR = State(SA[T, OR+1, OT], 1)
            oxT = State(SA[T, OR, OT+1], 1)
            if allowed(oxR)
                di = statetoi(self, oxR)
                add_edge!(graph, i, di, weights[1])
                add_edge!(graph, di, i, weights[2])
            end
            if allowed(oxT)
                di = statetoi(self, oxT)
                add_edge!(graph, i, di, weights[3])
                add_edge!(graph, di, i, weights[4])
            end
        end

        # Add edges for R-T and T-R
        for i in eachindex(states)
            (T, OR, OT) = states[i].extstate
            # R->T
            if T < N
                if OR < N - T   # at least one free R
                    dest = State(SA[T+1, OR, OT], 1)
                    di = statetoi(self, dest)
                    add_edge!(graph, i, di, weights[5])
                end
                if OR > 0       # at least one oxy R
                    dest = State(SA[T+1, OR-1, OT+1], 1)
                    di = statetoi(self, dest)
                    add_edge!(graph, i, di, weights[6])
                end
            end
            # T->R
            if T > 0
                if OT < T
                    dest = State(SA[T-1, OR, OT], 1)
                    di = statetoi(self, dest)
                    add_edge!(graph, i, di, weights[7])
                end
                if OT > 0
                    dest = State(SA[T-1, OR+1, OT-1], 1)
                    di = statetoi(self, dest)
                    add_edge!(graph, i, di, weights[8])
                end
            end
        end

        self
    end
end
function graph(agm::Allostery3D1)
    agm.graph
end
function states(agm::Allostery3D1)
    agm.states
end

function plotGM(agm::Allostery3D1, args...; ecolorbar=true, kwargs...)
    states_ = states(agm)
    layout = [Point(s.extstate[2], s.extstate[3], s.extstate[1]) for s in states_]

    node_color = [:snow4 for _ in 1:length(states_)]
    node_size = [15 for _ in 1:length(states_)]

    edge_color = [e.weight for e in edges(graph(agm))]

    fap = graphplot(agm.graph, args...; layout, node_color, node_size, edge_color, kwargs...)

    if ecolorbar
        edgecolorbar(fap)
    end

    fap
end

################################################################################
# Complex network allostery model
################################################################################
struct CAState{T}
    conformations::Vector{T}
    occupations::Vector{T}
end

struct ComplexAllostery <: AbstractGraphModel
    N::Int # Number of monomers
    C::Int # Number of possible conformational states per monomer
    B::Int # Number of ligand binding sites per monomer
end

################################################################################
# Sample run functions for the future when I don't remember how to use these
################################################################################
function run_sample_sfs()
    sfs = SimpleFourStateGM(5, 5; γi=0.05, γe=1.95)
    simGM(sfs, 10000000; delay=0.01)
end

function run_sample_agm()
    agm = Allostery3D1(6; weights=:random)
    simGM(agm, 1000; delay=0.05)
    # simGM(agm, 1000; delay=0.05, layout=Spring(dim=3))
end
