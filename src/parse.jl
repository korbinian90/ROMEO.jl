# The array syntax the command line tools accept:
#
#     [1.5,3.0]   [1.5, 3.0]   1.5 3.0   [1.5 3.0]   array, in any spacing
#     1:3         3.5:3.5:14                         range, with optional step
#     :                                              everything
#     3                                              a single number
#
# `nargs = '+'` hands the value over as several strings when the user separates
# them with spaces, which is why the input is a collection and gets joined.

"""
    parse_array(str)

Parse the array syntax the command line tools accept: an array literal
(`"[1.5,3.0]"`, brackets optional, comma or space separated), a range
(`"1:3"`, `"3.5:3.5:14"`), a bare `":"`, or a single number.

Returns `Colon()` for `":"`, the number itself for a single value, and a
`Vector{Int}` or `Vector{Float64}` otherwise - never a `Matrix` and never a lazy
range, so the result reads the same in a provenance record as it behaves in the
code. Throws an `ArgumentError` naming the offending text if it is not one of
those forms.
"""
function parse_array(str)
    s = String(strip(join(str, " ")))
    isempty(s) && throw(ArgumentError("expected a number, array or range, got nothing"))
    s == ":" && return Colon()
    # One layer of brackets is optional: "[1,2,3]" and "1,2,3" are the same.
    s = String(strip(strip(s), ['[', ']']))
    parts, isrange = _tokens(s)
    # All-integer input stays integer: echo indices are used to index with.
    ints = _parse_ints(parts)
    if ints !== nothing
        isrange && return _collect_range(ints, s)
        return length(ints) == 1 ? ints[1] : ints
    end
    floats = _parse_floats(parts, s)
    isrange && return _collect_range(floats, s)
    return length(floats) == 1 ? floats[1] : floats
end

# The parts of an array, separated by commas or spaces, or of a range,
# separated by colons. Written out so that no dynamic call is left.
function _tokens(s::String)
    isrange = any(==(':'), s)
    parts = String[]
    current = IOBuffer()
    for c in s
        if c == ':' || (!isrange && (c == ',' || isspace(c)))
            push!(parts, String(take!(current)))
        else
            write(current, c)
        end
    end
    push!(parts, String(take!(current)))
    if isrange
        parts = [String(strip(p)) for p in parts]
    else
        parts = filter(!isempty, parts)
    end
    return parts, isrange
end

function _parse_ints(parts)
    ints = Int[]
    for p in parts
        v = tryparse(Int, p)
        v === nothing && return nothing
        push!(ints, v)
    end
    return ints
end

function _parse_floats(parts, whole)
    floats = Float64[]
    for p in parts
        v = tryparse(Float64, p)
        v === nothing && throw(ArgumentError("\"$p\" in \"$whole\" is not a number"))
        push!(floats, v)
    end
    return floats
end

function _collect_range(bounds::Vector, whole)
    length(bounds) == 2 && return collect(bounds[1]:bounds[2])
    length(bounds) == 3 && return collect(bounds[1]:bounds[2]:bounds[3])
    throw(ArgumentError("\"$whole\" is not a range: expected start:stop or start:step:stop"))
end

"""
    parse_weight_flags(str)

Parse ROMEO's `--weights` bit string (`"1010"`, up to six flags) into a
`BitVector` of length 6. Returns `nothing` when `str` is not a bit string, which
is how the caller tells a set of flags from a named weighting like `"romeo3"`.
"""
function parse_weight_flags(str)
    (isempty(str) || length(str) > 6) && return nothing
    all(c -> c == '0' || c == '1', str) || return nothing
    flags = falses(6)
    for (i, c) in enumerate(str)
        flags[i] = c == '1'
    end
    return flags
end
