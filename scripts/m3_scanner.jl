using DrWatson
@quickactivate "TopoStochSim"

using JLD2

using ComplexAllostery
using ComplexAllostery.Model3

using LinearAlgebra
using Base.Threads

################################################################################
# General and fast, threaded scan function for the simplified model 3
################################################################################
function general_m3_scan(N, observables, save_fname=nothing,
    obs_names=[f"obs_{i}" for i in 1:length(observables)];
    r1s=[1.0],
    r2s=[1.0],
    r3s=[1.0],
    alphas=[0.0],
    ebs=[0.0],
    cPs=[1.0],
    cRs=[1.0],
    do_ss=true,
    save_metadata=(;)
)
    ca = makeCAGM(N, vars_simplified())
    cca_fact::Function = make_factory(
        ca,
        [symvars.sim_rs; symvars.energies[3];
            symvars.concentrations[1]; symvars.redconcentration]
    )

    cis = CartesianIndices(length.((r1s, r2s, r3s, alphas, ebs, cPs, cRs)))

    data = [fill(0.0, cis.indices) for _ in observables]

    # By default this is 11 but it seems much much better to parallelize here
    # instead of having a parallel eigensolver
    prev_blas_num_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)

    @threads for ci in cis
        r1_i, r2_i, r3_i, alpha_i, eb_i, cP_i, cR_i = ci.I
        cca = cca_fact(r1s[r1_i], r2s[r2_i], r3s[r3_i], alphas[alpha_i], ebs[eb_i], cPs[cP_i], cRs[cR_i])

        for (darray, observable) in zip(data, observables)
            if do_ss
                ss = supersteadystate(cca)
                darray[ci] = observable(cca, ss)
            else
                darray[ci] = observable(cca)
            end
        end
    end

    BLAS.set_num_threads(prev_blas_num_threads)

    if !isnothing(save_fname)
        obs_kwargs = Symbol.(obs_names) .=> data
        jldsave(save_fname; ca, r1s, r2s, r3s, alphas, ebs, cPs, cRs,
            observable_names=obs_names, obs_kwargs..., save_metadata...
        )
    end

    data
end

"""
Constructs an observable for `general_m3_scan` that is just the ss probability
sum over sum nodes.
"""
struct ObsTotalAreaProbability
    area_indices::Vector{Int}
end
function (tap::ObsTotalAreaProbability)(cca, ss)
    sum(@view ss[tap.area_indices])
end

# Some really basic observables
function make_allligs_obs(N, C=2, B=1)
    ai = Int[]
    for (i, st) in enumerate(ComplexAllostery.allstates_complexallostery(N, C, B))
        if all(x -> x == B, st.occupations)
            push!(ai, i)
        end
    end
    ObsTotalAreaProbability(ai)
end
function make_noligs_obs(N, C=2, B=1)
    ai = Int[]
    for (i, st) in enumerate(ComplexAllostery.allstates_complexallostery(N, C, B))
        if all(x -> x == 0, st.occupations)
            push!(ai, i)
        end
    end
    ObsTotalAreaProbability(ai)
end
function make_alltense_obs(N, C=2, B=1)
    ai = Int[]
    for (i, st) in enumerate(ComplexAllostery.allstates_complexallostery(N, C, B))
        if all(x -> x == 1, st.conformations)
            push!(ai, i)
        end
    end
    ObsTotalAreaProbability(ai)
end
function make_notense_obs(N, C=2, B=1)
    ai = Int[]
    for (i, st) in enumerate(ComplexAllostery.allstates_complexallostery(N, C, B))
        if all(x -> x != 1, st.conformations)
            push!(ai, i)
        end
    end
    ObsTotalAreaProbability(ai)
end
function make_onboundary_obs(N, C=2, B=1)
    ai = Int[]
    for (i, st) in enumerate(ComplexAllostery.allstates_complexallostery(N, C, B))
        if all(x -> x == B, st.occupations)
            push!(ai, i)
        elseif all(x -> x == 0, st.occupations)
            push!(ai, i)
        elseif all(x -> x == 1, st.conformations)
            push!(ai, i)
        elseif all(x -> x != 1, st.conformations)
            push!(ai, i)
        end
    end
    ObsTotalAreaProbability(ai)
end
