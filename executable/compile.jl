using PackageCompiler

fnmodule = "executable/sharedLibrary/UnwrappingC.jl"
builddir = "~/builddir"
builddir_win = "C:/builddir"
snoopfile = "executable/unwrapping/snoopfile.jl" #"bqunwrap/unwrapTest2.jl"
fncprogram = "executable/c_program.c"

build_executable("executable/unwrapping/Module.jl", "ROMEO"; builddir = builddir, snoopfile = snoopfile, cpu_target = "x86-64")
