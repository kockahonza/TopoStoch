timestamp() = Dates.format(Dates.now(), "yymmdd_HMS")
export timestamp

function roundrepr(x)
    f"{x:.3g}"
end
roundrepr(x::Num) = repr(x)
export roundrepr

function inc_edge!(g, u, v, w)
    add_edge!(g, u, v, get_weight(g, u, v) + w)
end
export inc_edge!

function add_edge_ifnotself!(g, v, w)
    if v != w
        add_edge!(g, v, w)
    end
end
export add_edge_ifnotself!

# NOTE: This would definitely be better as a macro but idk how to write them
"""
Takes a mutating function along with any args and kwargs and returns a
closure that will return a copy of the passed object to which the passed
function was applied along with args and kwargs.
"""
function copyand(f!, args...; kwargs...)
    function (obj)
        cobj = copy(obj)
        f!(cobj, args...; kwargs...)
        cobj
    end
end
export copyand

function make_grid(n;
    aspect_ratio=1.0, prioritize=:none,
    min_rows=1, max_rows=typemax(Int),
    min_cols=1, max_cols=typemax(Int)
)
    # Calculate initial number of columns based on aspect ratio
    cols = ceil(Int, sqrt(n * aspect_ratio))
    rows = ceil(Int, n / cols)

    # Adjust based on prioritization
    if prioritize == :rows
        rows = min(max(min_rows, rows), max_rows)
        cols = ceil(Int, n / rows)
    elseif prioritize == :cols
        cols = min(max(min_cols, cols), max_cols)
        rows = ceil(Int, n / cols)
    end

    # Ensure the grid fits within the specified limits
    rows = min(max(min_rows, rows), max_rows)
    cols = min(max(min_cols, cols), max_cols)

    return rows, cols
end
export make_grid
