# Parsing the array syntax the command line tools accept, without `eval`.
#
# Every tool used to call `eval(Meta.parse(...))` on whatever the user typed
# after -t, -e or -s. That runs arbitrary Julia from the command line, it keeps
# the tools from ever being statically compiled, and it leaves the raw strings
# in the provenance record where the resolved numbers belong. One parser here,
# in the package with no dependencies, replaces all of it.
#
# The accepted syntax is what the tools' help text has always promised:
#
#     [1.5,3.0]   [1.5, 3.0]   1.5 3.0   [1.5 3.0]   array, in any spacing
#     1:3         3.5:3.5:14                         range, with optional step
#     :                                              everything
#     3                                              a single number
#
# `nargs = '+'` hands the value over as several strings when the user separates
# them with spaces, which is why the input is a collection and gets joined
# before parsing.

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
    s = strip(join(str, " "))
    isempty(s) && throw(ArgumentError("expected a number, array or range, got nothing"))
    s == ":" && return Colon()
    # One layer of brackets is optional: "[1,2,3]" and "1,2,3" are the same.
    s = strip(strip(s), ['[', ']'])
    if occursin(':', s)
        parts = strip.(split(s, ':'))
        bounds = _parse_number.(parts, s)
        length(bounds) == 2 && return _collect_range(bounds[1]:bounds[2])
        length(bounds) == 3 && return _collect_range(bounds[1]:bounds[2]:bounds[3])
        throw(ArgumentError("\"$s\" is not a range: expected start:stop or start:step:stop"))
    end
    parts = filter(!isempty, split(s, r"[,\s]+"))
    nums = _parse_number.(parts, s)
    length(nums) == 1 && return nums[1]
    return _narrow(nums)
end

function _parse_number(part, whole)
    v = tryparse(Int, part)
    v === nothing || return v
    v = tryparse(Float64, part)
    v === nothing || return v
    throw(ArgumentError("\"$part\" in \"$whole\" is not a number"))
end

# An all-integer input stays integer, as it did when this went through Julia's
# own parser: echo indices are used to index with, and [1.0, 2.0] would not.
_narrow(nums) = all(n -> n isa Int, nums) ? Vector{Int}(nums) : Vector{Float64}(nums)
_collect_range(r) = _narrow(collect(r))

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
