import Base: repr, string

mutable struct Code
    first::Char
    last::Char
    code::Vector{Char}
    function Code(n::Int=1, first='a', last='z')
        c = new(first, last, [first])
        for _ in 2:n
            increment!(c)
        end
        c
    end
end

string(c::Code) = string(c.code...)
repr(c::Code) = string(c)

function increment!(c::Code)
    last = c.code[end]
    if last == c.last
        push!(c.code, c.first)
    else
        c.code[end] = last + 1
    end
    c
end
