using DrWatson
@quickactivate "cluster_env"

include("../../scripts/m3_scanner.jl")

function bigrun1(N=4; testrun=false)
    observables = [make_allligs_obs(N), make_noligs_obs(N), make_alltense_obs(N), make_notense_obs(N), make_onboundary_obs(N)]
    obs_names = ["allligs", "noligs", "alltense", "notense", "onboundary"]

    fname = f"./bigrun1_N{N}.jld2"
    if ispath(fname)
        @error f"Path {fname} already exists, exiting so that it is not overriden"
    end
    if !isdir(dirname(fname))
        mkpath(dirname(fname))
    end

    if !testrun
        general_m3_scan(N, observables, fname, obs_names;
            r1s=10 .^ LinRange(-5, 3, 15),
            r2s=10 .^ LinRange(-5, 3, 15),
            r3s=10 .^ LinRange(-5, 3, 15),
            alphas=[0, 0.5],
            ebs=LinRange(0, 6, 5),
            cPs=10 .^ LinRange(-5, 1, 15),
            cRs=10 .^ LinRange(-1, 2, 15)
        )
    else
        general_m3_scan(N, observables, fname, obs_names;
            r1s=10 .^ LinRange(-5, 3, 15),
            alphas=[0, 0.5],
            ebs=LinRange(0, 6, 5),
        )
    end
end
