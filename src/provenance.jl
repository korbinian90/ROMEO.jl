# Provenance records: what a run was told to do, and what has to be cited for the
# methods it actually ran.
#
# The writer lives here, in the package with no dependencies, so every package in
# the family can reach it without adding a version constraint pointing back up
# the dependency graph. Citation text follows ownership instead: each package
# registers the references for the methods it implements, at load time. ROMEO
# therefore knows about ROMEO and best-path unwrapping and nothing else;
# MriResearchTools adds MCPC-3D-S and the homogeneity correction, CLEARSWI adds
# CLEAR-SWI, and so on. That keeps one definition per reference while letting the
# dependency arrows all point the same way.

"""
    CITATIONS

Reference text for the methods in the toolbox, keyed by method. ROMEO seeds this
with its own; other packages add theirs through [`register_citation!`](@ref).
"""
const CITATIONS = Dict{Symbol,String}()

"""
    NOTICES

Non-citation facts a method carries, keyed like [`CITATIONS`](@ref). Written into
the citations file when that method ran, because that file is what someone reads
before publishing or before shipping a product.
"""
const NOTICES = Dict{Symbol,String}()

"""
    LABELS

Human-readable method name for each key, used as the heading above its reference
in the citations file. A reader should be able to tell which step of the run each
reference is for without recognising the paper.
"""
const LABELS = Dict{Symbol,String}()

"""
    register_citation!(key, text; notice=nothing, label=nothing)

Register the reference for a method, to be written by [`write_provenance`](@ref)
when that method is used. Call this from a package's `__init__` for the methods
that package implements, so the text lives with the code rather than in whichever
tool happens to print it.

`label` is the method name shown as a heading above the reference, so a reader
can see which step of the run it belongs to. Two keys may share a label, in which
case their references appear together under one heading - which is how a method
with more than one reference is expressed. Defaults to the key.

Re-registering the same key with identical text is a no-op; changing the text of
an existing key warns, because two packages disagreeing about a reference is a
bug worth hearing about.
"""
function register_citation!(key::Symbol, text::AbstractString; notice=nothing, label=nothing)
    text = _dedent(text)
    if haskey(CITATIONS, key) && CITATIONS[key] != text
        @warn "citation for :$key re-registered with different text; keeping the first" key
    else
        CITATIONS[key] = text
    end
    if notice !== nothing
        NOTICES[key] = _dedent(notice)
    end
    if label !== nothing
        LABELS[key] = String(label)
    end
    return key
end

_label(key) = get(LABELS, key, string(key))

# Reference text is written indented for readability in the source; strip that
# back out so the written file is left-aligned.
function _dedent(text)
    lines = split(text, '\n')
    join([first(lines); lstrip.(lines[2:end])], '\n')
end

function __init__()
    register_citation!(:romeo, """Dymerska, B., Eckstein, K., Bachrata, B., Siow, B., Trattnig, S., Shmueli, K., Robinson, S.D., 2020.
                                  Phase Unwrapping with a Rapid Opensource Minimum Spanning TreE AlgOrithm (ROMEO).
                                  Magnetic Resonance in Medicine.
                                  https://doi.org/10.1002/mrm.28563""";
                       label = "ROMEO Unwrapping")
    register_citation!(:bestpath, """Abdul-Rahman, H.S., Gdeisat, M.A., Burton, D.R., Lalor, M.J., Lilley, F., Moore, C.J., 2007.
                                     Fast and robust three-dimensional best path phase unwrapping algorithm.
                                     Applied Optics 46, 6623-6635.
                                     https://doi.org/10.1364/AO.46.006623""";
                       label = "Best-path Unwrapping")
    register_citation!(:julia, """Bezanson, J., Edelman, A., Karpinski, S., Shah, V.B., 2017.
                                  Julia: A fresh approach to numerical computing.
                                  SIAM Review 59, 65-98.
                                  https://doi.org/10.1137/141000671""";
                       label = "Julia Scientific Programming Language")
end

# ArgParse hands a multi-value option over as the raw strings the user typed, so
# print the text for those and keep the bracketed form for resolved numeric ones.
function _fmt(v::AbstractArray)
    isempty(v) && return "(not set)"
    all(x -> x isa AbstractString, v) && return join(v, " ")
    return string(collect(v))
end
_fmt(v::AbstractString) = isempty(v) ? "(not set)" : String(v)
_fmt(v) = string(v)

function _timestamp()
    # Deliberately not using Dates: this package has no dependencies and a
    # provenance record is not worth adding one for.
    t = round(Int, time())
    return string(Libc.strftime("%Y-%m-%dT%H:%M:%S", t), " (local), unix ", t)
end

_default_describe(path) = try abspath(String(path)) catch; String(path) end

"""
    write_provenance(dir, tool; version, args, settings, cite, optional=Symbol[],
                     inputs=(), packages=(), describe=abspath)

Write `settings_<tool>.txt` and `citations_<tool>.txt` into `dir`, recording what
was run and what has to be cited for it.

* `version` the application version string.
* `args` the raw command line arguments.
* `settings` any key-value collection of resolved settings. Written sorted, and
  array values are written out rather than skipped, so the echo times actually
  used are recoverable from the record.
* `cite` the methods that were actually used, as keys of [`CITATIONS`](@ref).
  Pass only what ran: a citation for a method the user did not use is as wrong as
  a missing one. Any [`NOTICES`](@ref) entry for those methods is included too.
* `optional` further methods to list under "Optional citations".
* `inputs` `name => path` pairs for the input files.
* `packages` modules whose versions did the work, recorded alongside the Julia
  version so a result can be traced to the code that produced it.
* `describe` how to render an input path. Defaults to the absolute path; callers
  that can read the file cheaply pass something that adds the dimensions.

# Examples
```julia-repl
julia> write_provenance("out", "romeo"; version="4.7.1", args=ARGS, settings,
                        cite=[:romeo, :aspire], inputs=["phase" => fn_phase],
                        packages=[ROMEO])
```
"""
function write_provenance(dir, tool; version, args, settings, cite,
                          optional=Symbol[], inputs=(), packages=(),
                          describe=_default_describe)
    dir = abspath(dir)
    mkpath(dir)
    _write_settings(dir, tool; version, args, settings, inputs, packages, describe)
    write_citations(dir, tool; cite, optional)
    return dir
end

"""
    package_version(m::Module)

The version of the package `m`, as a string, for a provenance record.

`pkgversion` reads a package's Project.toml through path metadata that the
module carries at runtime, and it returns `nothing` rather than throwing when
that metadata is not there, which is how a settings file ends up recording the
literal word "nothing". A sysimage built with `--strip-metadata` is one case;
a module defined outside a package is another. Packages in this family
therefore also record their own version in a `PKG_VERSION` constant, which is
evaluated while the package is precompiled and is part of the image after
that, so it survives both.
"""
function package_version(m::Module)
    if isdefined(m, :PKG_VERSION)
        return string(getglobal(m, :PKG_VERSION))
    end
    v = try pkgversion(m) catch; nothing end
    return isnothing(v) ? "unknown" : string(v)
end

function _write_settings(dir, tool; version, args, settings, inputs, packages, describe)
    open(joinpath(dir, "settings_$tool.txt"), "w") do io
        println(io, "# $tool $version")
        println(io, "# written: ", _timestamp())
        println(io, "# command: ", join(args, " "))
        println(io)

        println(io, "[versions]")
        println(io, "julia: ", VERSION)
        for m in packages
            println(io, nameof(m), ": ", package_version(m))
        end
        println(io)

        if !isempty(inputs)
            println(io, "[inputs]")
            for (name, path) in inputs
                isnothing(path) && continue
                println(io, name, ": ", describe(path))
            end
            println(io)
        end

        println(io, "[settings]")
        for key in sort(collect(keys(settings)); by=string)
            key == "header" && continue
            println(io, key, ": ", _fmt(settings[key]))
        end
    end
end

"""
    write_citations(dir, tool; cite, optional=Symbol[])

Write `citations_<tool>.txt` into `dir`, covering only the methods named in
`cite`. Keys are looked up in the citation registry, which each package fills in
for the methods it implements (see [`register_citation!`](@ref)); a key with no
registered citation is warned about rather than silently omitted, because a
missing reference is the failure this is meant to prevent. Any notice attached to
a used method is written below the references, and `optional` keys that are
registered but were not used are listed separately.

The registry itself is `ROMEO.CITATIONS`, `ROMEO.NOTICES` and `ROMEO.LABELS`.
None are exported: the names are too generic to put in every user's
namespace, and `register_citation!` plus this function are the intended
interface.
"""
function write_citations(dir, tool; cite, optional=Symbol[])
    known = filter(k -> haskey(CITATIONS, k), unique(cite))
    missing_keys = setdiff(unique(cite), known)
    if !isempty(missing_keys)
        @warn """No citation is registered for $(join(missing_keys, ", ")), so it is missing from the record.
                 Each package registers the references for the methods it implements, so this means the
                 owning package is either not loaded or too old to register it. Check its version.""" maxlog=1
    end
    open(joinpath(dir, "citations_$tool.txt"), "w") do io
        println(io, "# Citations for the methods this run actually used.")
        println(io, "# Methods that were available but not used are deliberately absent.")
        println(io)
        _write_labelled(io, known)

        notices = [k for k in known if haskey(NOTICES, k)]
        if !isempty(notices)
            println(io, "# Notices for the methods used:")
            println(io)
            for k in notices
                println(io, NOTICES[k])
                println(io)
            end
        end

        opt = filter(k -> haskey(CITATIONS, k) && k ∉ known, unique(optional))
        if !isempty(opt)
            println(io, "# Optional citations:")
            println(io)
            _write_labelled(io, opt)
        end
    end
end

# One heading per method, not per reference: consecutive keys sharing a label are
# two references for the same method, and belong under one heading.
function _write_labelled(io, ks)
    previous = nothing
    for k in ks
        label = _label(k)
        if label != previous
            println(io, "## ", label)
            previous = label
        end
        println(io, CITATIONS[k])
        println(io)
    end
end
