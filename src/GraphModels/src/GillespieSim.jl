"""
Gillespie algorithm that should actually be fast and parallelizeable, save
output in jld2 files which includes the whole trajectory or observables,
optionally does some plotting
"""
module GillespieSim
using ..GraphModels

mutable struct Gillespie
    # The model info in a fast, accesible manner
    graph::AbstractGraph{Int}
    # Saving data
    # Current state
    cur_time::Float64
    cur_state::Int
end
# export Gillespie

function new_gillespie_run(
    gm::AbstractGraphModel{<:AbstractFloat}
)
    Gillespie(graph(gm), 0, rand(1:numstates(gm)))
end
export new_gillespie_run

"""
Returns the time until next transition and the destination
"""
function next_step(gs::Gillespie)
    (rand(), rand(neighbors(gs.graph, gs.cur_state)))
end
export next_step

function run_gillespie!(gs::Gillespie;
    exit_time=nothing,
    exit_steps=nothing
)
    sim_steps = 0
    sim_time = 0.0
    exit = false
    while !exit
        # Make a step
        (delta_time, next_state) = next_step(gs)
        gs.cur_time += delta_time
        gs.cur_state = next_state
        sim_steps += 1
        sim_time += delta_time

        # Do saving/callbacks
        # @show gs.cur_time, gs.cur_state

        # Check exit conditions
        if !isnothing(exit_time)
            if sim_time >= exit_time
                # @info "Exited due to max simulation time met"
                exit = true
            end
        end
        if !isnothing(exit_steps)
            if sim_steps >= exit_steps
                # @info "Exited due to max simulation number of steps met"
                exit = true
            end
        end
    end
end
export run_gillespie!

end
