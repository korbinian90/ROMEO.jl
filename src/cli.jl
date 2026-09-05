"""
    CLI

A small command line parser for the tools of this family. It has no
dependencies and no dynamic dispatch, so that a tool can be compiled with
`juliac --trim`. It covers what the tools use of ArgParse: long and short
names, flags, options with one value, with several values (`nargs='+'`) and
with an optional value (`nargs='?'`), `--name=value`, combined short flags
(`-qQ`), `--help` and `--version`.

    spec = CLI.Spec("romeo", version, [CLI.Option("--phase", "-p", "The phase image"), ...])
    values = CLI.parse(spec, args)   # Dict of the given options, or nothing after --help / --version
"""
module CLI

"""
    Option(long, short, help; nargs=:one)

An option with its `--long` name, `-s`hort name (or `""`) and help text.
`nargs` is `:none` for a flag, `:one` for a value, `:many` for one or more
values and `:optional` for a value that may be left out.
"""
struct Option
    long::String
    short::String
    help::String
    nargs::Symbol
end
Option(long, short, help; nargs=:one) = Option(long, short, help, nargs)

struct Spec
    program::String
    version::String
    options::Vector{Option}
end

"The parsed command line: the values of every option that was given, by name without dashes."
const Values = Dict{String,Vector{String}}

name(o::Option) = o.long[3:end]
metavar(o::Option) = uppercase(name(o))

function find_long(spec::Spec, long::AbstractString)
    for o in spec.options
        o.long == long && return o
    end
    return nothing
end

function find_short(spec::Spec, short::AbstractString)
    for o in spec.options
        o.short == short && return o
    end
    return nothing
end

# a token starting with a dash is an option name, unless it is a negative number
is_option(tok::AbstractString) = length(tok) > 1 && tok[1] == '-' && tryparse(Float64, tok) === nothing

"""
    parse(spec, args) -> Values or nothing

Parses `args`. Prints the help or version and returns `nothing` when asked for
them, throws an `ArgumentError` for a wrong command line.
"""
function parse(spec::Spec, args::AbstractVector{<:AbstractString})
    values = Values()
    i = firstindex(args)
    while i <= lastindex(args)
        tok = String(args[i])
        i += 1
        if tok == "--help" || tok == "-h"
            print_help(spec)
            return nothing
        elseif tok == "--version"
            print_stdout(spec.version * "\n")
            return nothing
        elseif startswith(tok, "--")
            long, value = split_assignment(tok)
            o = find_long(spec, long)
            o === nothing && throw(ArgumentError("unrecognized option $long"))
            if value === nothing
                i = consume!(values, o, args, i)
            else
                o.nargs == :none && throw(ArgumentError("option $long takes no value"))
                values[name(o)] = [value]
            end
        elseif is_option(tok)
            # -p value, -pvalue, and combined flags -qQ
            j = nextind(tok, 1)
            while j <= lastindex(tok)
                short = "-" * tok[j]
                o = find_short(spec, short)
                o === nothing && throw(ArgumentError("unrecognized option $short"))
                rest = tok[nextind(tok, j):end]
                if o.nargs == :none
                    values[name(o)] = String[]
                elseif isempty(rest)
                    i = consume!(values, o, args, i)
                else
                    values[name(o)] = [rest]
                    break
                end
                j = nextind(tok, j)
            end
        else
            throw(ArgumentError("too many arguments"))
        end
    end
    return values
end

# "--name=value" as ("--name", "value"), anything else as (tok, nothing)
function split_assignment(tok::String)
    for i in eachindex(tok)
        tok[i] == '=' && return tok[1:prevind(tok, i)], tok[nextind(tok, i):end]
    end
    return tok, nothing
end

# Takes the values of option `o` from args starting at i, returns the next i
function consume!(values::Values, o::Option, args, i)
    vals = String[]
    if o.nargs == :none
    elseif o.nargs == :one
        (i > lastindex(args) || is_option(args[i])) && throw(ArgumentError("option $(o.long) requires a value"))
        push!(vals, String(args[i]))
        i += 1
    else # :many, :optional
        while i <= lastindex(args) && !is_option(args[i])
            push!(vals, String(args[i]))
            i += 1
            o.nargs == :optional && break
        end
        o.nargs == :many && isempty(vals) && throw(ArgumentError("option $(o.long) requires at least one value"))
    end
    values[name(o)] = vals
    return i
end

# In a program compiled with juliac --trim, Base.stdout cannot be called through
# (its declared type is abstract) and logging is turned off. Both are recognised
# through the logger juliac replaces, and the raw streams are used instead.
static_binary() = Base.CoreLogging.current_logger_for_env(Base.CoreLogging.Error, :cli, @__MODULE__) === nothing
print_stdout(s::String) = static_binary() ? print(Core.stdout, s) : print(s)
print_stderr(s::String) = static_binary() ? print(Core.stderr, s) : print(stderr, s)

"Prints `msg` to standard output."
info(msg::String) = print_stdout(msg * "\n")

"Warns about `msg`: through the logging system, or on standard error in a static binary."
function warn(msg::String)
    if static_binary()
        print_stderr("Warning: " * msg * "\n")
    else
        @warn msg
    end
    return nothing
end

spaces(n::Int) = n <= 0 ? "" : String(fill(UInt8(' '), n))

# Words of `text` wrapped to `width` columns, each line prefixed with `indent`
function wrap(text::String, width::Int, indent::String)
    lines = String[]
    line = indent
    for word in split(text)
        if length(line) > length(indent) && length(line) + 1 + length(word) > width
            push!(lines, line)
            line = indent * word
        elseif length(line) > length(indent)
            line *= " " * word
        else
            line *= word
        end
    end
    push!(lines, line)
    return join(lines, "\n")
end

function usage_item(o::Option)
    n = isempty(o.short) ? o.long : o.short
    o.nargs == :none && return "[$n]"
    o.nargs == :one && return "[$n $(metavar(o))]"
    o.nargs == :many && return "[$n $(metavar(o))...]"
    return "[$n [$(metavar(o))]]"
end

function help_item(o::Option)
    names = isempty(o.short) ? o.long : o.short * ", " * o.long
    o.nargs == :none && return names
    o.nargs == :one && return names * " " * metavar(o)
    o.nargs == :many && return names * " " * metavar(o) * "..."
    return names * " [" * metavar(o) * "]"
end

function help(spec::Spec)
    io = IOBuffer()
    items = join((usage_item(o) for o in spec.options), " ")
    print(io, wrap("usage: $(spec.program) $items [--version] [-h]", 80, "  ")[3:end], "\n\n")
    print(io, "optional arguments:\n")
    for o in spec.options
        item = help_item(o)
        if length(item) <= 22
            print(io, "  ", item, spaces(22 - length(item)), wrap(o.help, 80, spaces(24))[25:end], "\n")
        else
            print(io, "  ", item, "\n", wrap(o.help, 80, spaces(26)), "\n")
        end
    end
    print(io, "  --version             show version information and exit\n")
    print(io, "  -h, --help            show this help message and exit\n")
    return String(take!(io))
end

print_help(spec::Spec) = print_stdout(help(spec))

"""
    format(x)

The text of a value for a settings record: numbers as they print, a vector of
numbers as `[a, b]`, a vector of strings joined by spaces, and `(not set)` for
an empty one. Does not use `show` of arrays, which cannot be compiled statically.
"""
format(x::Real) = string(x)
format(x::Bool) = x ? "true" : "false"
format(x::AbstractString) = isempty(x) ? "(not set)" : String(x)
format(::Nothing) = "(not set)"
format(v::AbstractVector{<:AbstractString}) = isempty(v) ? "(not set)" : join(v, " ")
format(v::AbstractVector{<:Real}) = isempty(v) ? "(not set)" : "[" * join((string(x) for x in v), ", ") * "]"

end # module
