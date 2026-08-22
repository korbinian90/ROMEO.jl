@testset "Provenance" begin

# ROMEO carries the writer and the references for its own methods, and nothing
# else. Anything more would mean a package knowing about methods it does not
# implement.
@test haskey(ROMEO.CITATIONS, :romeo)
@test haskey(ROMEO.CITATIONS, :bestpath)
@test !haskey(ROMEO.CITATIONS, :clearswi)
@test !haskey(ROMEO.CITATIONS, :tgv)

settings = Dict{String,Any}(
    "phase" => "p.nii",
    "TEs" => [4.0, 8.0, 12.0],           # an array: must be recorded, not skipped
    "weights" => "romeo3",
    "header" => "should not be written", # deliberately excluded
)

dir = mktempdir()
write_provenance(dir, "tool"; version="9.9.9", args=["-p", "p.nii", "-t", "[4,8,12]"],
                 settings, cite=[:romeo], optional=[:julia],
                 inputs=["phase" => "p.nii"], packages=[ROMEO])

s = read(joinpath(dir, "settings_tool.txt"), String)
c = read(joinpath(dir, "citations_tool.txt"), String)

@test occursin("# tool 9.9.9", s)
@test occursin("-p p.nii -t [4,8,12]", s)
@test occursin("julia: $VERSION", s)
@test occursin("ROMEO: $(pkgversion(ROMEO))", s)
# Array settings used to be dropped, which lost the echo times - the one setting
# a result can least afford to be missing.
@test occursin("TEs: [4.0, 8.0, 12.0]", s)
@test !occursin("should not be written", s)

# Citations must cover what ran and nothing else.
@test occursin("Phase Unwrapping with a Rapid Opensource", c)
@test !occursin("ASPIRE", c)
@test occursin("# Optional citations:", c)
# Reference text is stored indented for readability; the file must not be.
@test !occursin("\n   Magnetic Resonance in Medicine", c)

# Registering is idempotent for identical text, and complains when two packages
# disagree about a reference rather than silently keeping one.
n = length(ROMEO.CITATIONS)
register_citation!(:romeo, ROMEO.CITATIONS[:romeo])
@test length(ROMEO.CITATIONS) == n
@test_logs (:warn,) register_citation!(:romeo, "something else entirely")
@test occursin("Rapid Opensource", ROMEO.CITATIONS[:romeo]) # first registration kept

# A method whose owning package is not loaded must be reported, not silently
# dropped: a missing citation is the failure this file exists to prevent.
dir2 = mktempdir()
@test_logs (:warn,) write_provenance(dir2, "t2"; version="1", args=String[],
                                     settings=Dict{String,Any}(), cite=[:not_a_real_method])

# A notice registered with a citation travels with it, exactly once.
register_citation!(:test_method, "Some Reference."; notice="A NOTICE.")
dir3 = mktempdir()
write_provenance(dir3, "t3"; version="1", args=String[], settings=Dict{String,Any}(),
                 cite=[:test_method, :test_method])
c3 = read(joinpath(dir3, "citations_t3.txt"), String)
@test count("Some Reference.", c3) == 1
@test occursin("A NOTICE.", c3)
delete!(ROMEO.CITATIONS, :test_method); delete!(ROMEO.NOTICES, :test_method)

end
