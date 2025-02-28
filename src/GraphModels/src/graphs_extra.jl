"""
Returns true if there are no edges going from any vertex in `g` but not in `vs`
to any vertex in `vs`. This algorithm is O(|vs|^2*|g|)
"""
function areunreachable(g, vs)
    for v in vs
        for inn in inneighbors(g, v)
            if !(inn in vs)
                return false
            end
        end
    end
    return true
end
"""
Returns true if there are no edges going from any vertex `vs` to any edge in
`g` but not in `vs`. This algorithm is O(|vs|^2*|g|)
"""
function aredeadend(g, vs)
    for v in vs
        for outn in outneighbors(g, v)
            if !(outn in vs)
                return false
            end
        end
    end
    return true
end
export areunreachable, aredeadend

"""
Tests whether the vertices `vs` form a simple arrow in the graph `g`. This means
that the vertices are all connected by a single directed path (in the same
direction) with no loops and no edges to any other vertices.
"""
function aresimplearrow(g, vs)
    if length(vs) < 2
        @debug "calling aresimplearrow with a vertex vector of length 0 or 1, it's not clear what this should do, returning true"
        return true
    end

    # find the only vertex with no incoming edges
    starti = nothing
    for (vi, v) in enumerate(vs)
        ins = inneighbors(g, v)
        if length(ins) == 0
            if isnothing(starti)
                starti = vi
            else
                return false
            end
        elseif length(ins) != 1 # and also make sure there are no other incoming edges
            return false
        end
    end
    if isnothing(starti)
        return false
    end

    # keep track of which vertices we traversed already
    visited = falses(length(vs))

    # start traversing this maybe arrow
    visited[starti] = true
    v = vs[starti]
    while true
        ons = outneighbors(g, v)
        if length(ons) == 0 # we reached the end of this simple path
            break
        elseif length(ons) != 1 # if there's a branch it's not an arrow
            return false
        end
        # check if we've visited this edge before, if so it's not an arrow
        newvi = findfirst(v2 -> v2 == ons[1], vs)
        if visited[newvi]
            return false
        end
        # move to the next vertex
        visited[newvi] = true
        v = ons[1]
    end

    # finally this is an arrow if we've visited all passed vertices
    return all(visited)
end
export aresimplearrow
