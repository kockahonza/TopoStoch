using DrWatson
@quickactivate "TopoStochSim"

using Printf
using DataFrames
using GLMakie

using GraphModels
using NonEqDigits

const ca_fs_codes = [77, 232, 150, 105, 178, 23, 204, 51]
const ca_ntfs_codes = ca_fs_codes[begin:end-2]

function obs_numacs(gm)
    length(attracting_components(graph(gm)))
end

function obs_numntsccs(gm)
    length(filter(x -> length(x) != 1, strongly_connected_components(graph(gm))))
end

function obs_numssccs(gm)
    x, _ = findsingularsccs(gm, strongly_connected_components(graph(gm)))
    length(x)
end
function obs_numnssccs(gm)
    _, x = findsingularsccs(gm, strongly_connected_components(graph(gm)))
    length(x)
end

function obs_numntacs(gm)
    length(filter(x -> length(x) != 1, attracting_components(graph(gm))))
end

function fsplot1(codes=ca_fs_codes)
    plot_eachcode_vsN2(codes, 3:8, [
        "num sccs" => x -> length(strongly_connected_components(graph(x))),
        "num singular sccs" => obs_numssccs,
        "num non-singular sccs" => obs_numnssccs,
        "num nt acs" => obs_numntacs,
    ])
end
