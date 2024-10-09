var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = ROMEO","category":"page"},{"location":"#ROMEO","page":"Home","title":"ROMEO","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for ROMEO.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [ROMEO]","category":"page"},{"location":"#ROMEO.calculateweights-Tuple{Any}","page":"Home","title":"ROMEO.calculateweights","text":"calculateweights(wrapped; weights=:romeo, kwargs...)\n\nCalculates weights for all edges. size(weights) == [3, size(wrapped)...]\n\nOptional keyword arguments:\n\nweights: Options are [:romeo] | :romeo2 | :romeo3 | :bestpath.\nmag: Additional mag weights are used.\nmask: Unwrapping is only performed inside the mask.\nphase2: A second reference phase image (possibly with different echo time).  It is used for calculating the phasecoherence weight.\nTEs: The echo times of the phase and the phase2 images as a tuple (eg. (5, 10) or [5, 10]).\n\n\n\n\n\n","category":"method"},{"location":"#ROMEO.unwrap","page":"Home","title":"ROMEO.unwrap","text":"unwrap(wrapped::AbstractArray; keyargs...)\n\nROMEO unwrapping for 3D and 4D data.\n\nOptional keyword arguments:\n\nTEs: Required for 4D data. The echo times for multi-echo data. In the case of single-echo    application with phase and the phase2 as a tuple (eg. (5, 10) or [5, 10]).\nweights: Options are [:romeo] | :romeo2 | :romeo3 | :bestpath.\nmag: The magnitude is used to improve the unwrapping-path.\nmask: Unwrapping is only performed inside the mask.\nphase2: A second reference phase image (possibly with different echo time).   It is used for calculating the phasecoherence weight. This is automatically   done for 4D multi-echo input and therefore not required.\ncorrectglobal=false: If true corrects global n2π offsets.\nindividual=false: If true perform individual unwrapping of echos.   Type ?unwrap_individual for more information\ntemplate=1: echo that is spatially unwrapped (if individual is false)\nmaxseeds=1: higher values allow more seperate regions\nmerge_regions=false: spatially merge neighboring regions after unwrapping\ncorrect_regions=false: bring each regions median closest to 0 by adding n2π\nwrap_addition=0: [0;π], allows 'linear unwrapping', neighbors can have more   (π+wrap_addition) phase difference\ntemporal_uncertain_unwrapping=0.5: Useful for veins. It uses spatial   unwrapping on voxels that have high uncertainty values after temporal   unwrapping. A higher quality threshold re-unwraps more voxels spatially.\n\nExamples\n\njulia> using MriResearchTools\njulia> phase = readphase(\"phase_3echo.nii\")\njulia> unwrapped = unwrap(phase; TEs=[1,2,3])\njulia> savenii(unwrapped, \"unwrapped.nii\"; header=header(phase))\n\n\n\n\n\n","category":"function"},{"location":"#ROMEO.unwrap!","page":"Home","title":"ROMEO.unwrap!","text":"unwrap(wrapped::AbstractArray; keyargs...)\n\nROMEO unwrapping for 3D and 4D data.\n\nOptional keyword arguments:\n\nTEs: Required for 4D data. The echo times for multi-echo data. In the case of single-echo    application with phase and the phase2 as a tuple (eg. (5, 10) or [5, 10]).\nweights: Options are [:romeo] | :romeo2 | :romeo3 | :bestpath.\nmag: The magnitude is used to improve the unwrapping-path.\nmask: Unwrapping is only performed inside the mask.\nphase2: A second reference phase image (possibly with different echo time).   It is used for calculating the phasecoherence weight. This is automatically   done for 4D multi-echo input and therefore not required.\ncorrectglobal=false: If true corrects global n2π offsets.\nindividual=false: If true perform individual unwrapping of echos.   Type ?unwrap_individual for more information\ntemplate=1: echo that is spatially unwrapped (if individual is false)\nmaxseeds=1: higher values allow more seperate regions\nmerge_regions=false: spatially merge neighboring regions after unwrapping\ncorrect_regions=false: bring each regions median closest to 0 by adding n2π\nwrap_addition=0: [0;π], allows 'linear unwrapping', neighbors can have more   (π+wrap_addition) phase difference\ntemporal_uncertain_unwrapping=0.5: Useful for veins. It uses spatial   unwrapping on voxels that have high uncertainty values after temporal   unwrapping. A higher quality threshold re-unwraps more voxels spatially.\n\nExamples\n\njulia> using MriResearchTools\njulia> phase = readphase(\"phase_3echo.nii\")\njulia> unwrapped = unwrap(phase; TEs=[1,2,3])\njulia> savenii(unwrapped, \"unwrapped.nii\"; header=header(phase))\n\n\n\n\n\n","category":"function"},{"location":"#ROMEO.unwrap_individual","page":"Home","title":"ROMEO.unwrap_individual","text":"unwrap_individual(wrapped::AbstractArray{T,4}; TEs, keyargs...) where T\n\nPerforms individual unwrapping of the echoes instead of temporal unwrapping. Still uses multi-echo information to improve the quality map. This function is identical to unwrap with the flag individual=true. The syntax is identical to unwrap, but doesn't support the temporal_uncertain_unwrapping and template options: Base.Docs.DocStr(svec(\"    unwrap(wrapped::AbstractArray; keyargs...)\\n\\nROMEO unwrapping for 3D and 4D data.\\n\\n###  Optional keyword arguments:\\n\\n- TEs: Required for 4D data. The echo times for multi-echo data. In the case of single-echo \\n    application with phase and the phase2 as a tuple (eg. (5, 10) or [5, 10]).\\n- weights: Options are [:romeo] | :romeo2 | :romeo3 | :bestpath.\\n- mag: The magnitude is used to improve the unwrapping-path.\\n- mask: Unwrapping is only performed inside the mask.\\n- phase2: A second reference phase image (possibly with different echo time).\\n    It is used for calculating the phasecoherence weight. This is automatically\\n    done for 4D multi-echo input and therefore not required.\\n- correctglobal=false: If true corrects global n2π offsets.\\n- individual=false: If true perform individual unwrapping of echos.\\n    Type ?unwrap_individual for more information\\n- template=1: echo that is spatially unwrapped (if individual is false)\\n- maxseeds=1: higher values allow more seperate regions\\n- merge_regions=false: spatially merge neighboring regions after unwrapping\\n- correct_regions=false: bring each regions median closest to 0 by adding n2π\\n- wrap_addition=0: [0;π], allows 'linear unwrapping', neighbors can have more\\n    (π+wrapaddition) phase difference\\n- `temporaluncertainunwrapping=0.5: Useful for veins. It uses spatial\\n    unwrapping on voxels that have high uncertainty values after temporal\\n    unwrapping. A higher quality threshold re-unwraps more voxels spatially.\\n\\n# Examples\\n``julia-repl\\njulia> using MriResearchTools\\njulia> phase = readphase(\\\"phase3echo.nii\\\")\\njulia> unwrapped = unwrap(phase; TEs=[1,2,3])\\njulia> savenii(unwrapped, \\\"unwrapped.nii\\\"; header=header(phase))\\n\\n\"), unwrap(wrapped::AbstractArray; keyargs...)\n\n\nROMEO unwrapping for 3D and 4D data.\n\n### Optional keyword arguments:\n\n  * `TEs`: Required for 4D data. The echo times for multi-echo data. In the case of single-echo    application with phase and the phase2 as a tuple (eg. (5, 10) or [5, 10]).\n  * `weights`: Options are [`:romeo`] | `:romeo2` | `:romeo3` | `:bestpath`.\n  * `mag`: The magnitude is used to improve the unwrapping-path.\n  * `mask`: Unwrapping is only performed inside the mask.\n  * `phase2`: A second reference phase image (possibly with different echo time).   It is used for calculating the phasecoherence weight. This is automatically   done for 4D multi-echo input and therefore not required.\n  * `correctglobal=false`: If `true` corrects global n2π offsets.\n  * `individual=false`: If `true` perform individual unwrapping of echos.   Type `?unwrap_individual` for more information\n  * `template=1`: echo that is spatially unwrapped (if `individual` is `false`)\n  * `maxseeds=1`: higher values allow more seperate regions\n  * `merge_regions=false`: spatially merge neighboring regions after unwrapping\n  * `correct_regions=false`: bring each regions median closest to 0 by adding n2π\n  * `wrap_addition=0`: [0;π], allows 'linear unwrapping', neighbors can have more   (π+wrap_addition) phase difference\n  * `temporal_uncertain_unwrapping=0.5`: Useful for veins. It uses spatial   unwrapping on voxels that have high uncertainty values after temporal   unwrapping. A higher quality threshold re-unwraps more voxels spatially.\n\n# Examples\n\n\njulia-repl julia> using MriResearchTools julia> phase = readphase(\"phase_3echo.nii\") julia> unwrapped = unwrap(phase; TEs=[1,2,3]) julia> savenii(unwrapped, \"unwrapped.nii\"; header=header(phase)) ``` , Dict{Symbol, Any}(:typesig => Union{}, :module => ROMEO, :linenumber => 28, :binding => ROMEO.unwrap!, :path => \"/home/runner/work/ROMEO.jl/ROMEO.jl/src/unwrapping.jl\"))\n\n\n\n\n\n","category":"function"},{"location":"#ROMEO.unwrap_individual!","page":"Home","title":"ROMEO.unwrap_individual!","text":"unwrap_individual(wrapped::AbstractArray{T,4}; TEs, keyargs...) where T\n\nPerforms individual unwrapping of the echoes instead of temporal unwrapping. Still uses multi-echo information to improve the quality map. This function is identical to unwrap with the flag individual=true. The syntax is identical to unwrap, but doesn't support the temporal_uncertain_unwrapping and template options: Base.Docs.DocStr(svec(\"    unwrap(wrapped::AbstractArray; keyargs...)\\n\\nROMEO unwrapping for 3D and 4D data.\\n\\n###  Optional keyword arguments:\\n\\n- TEs: Required for 4D data. The echo times for multi-echo data. In the case of single-echo \\n    application with phase and the phase2 as a tuple (eg. (5, 10) or [5, 10]).\\n- weights: Options are [:romeo] | :romeo2 | :romeo3 | :bestpath.\\n- mag: The magnitude is used to improve the unwrapping-path.\\n- mask: Unwrapping is only performed inside the mask.\\n- phase2: A second reference phase image (possibly with different echo time).\\n    It is used for calculating the phasecoherence weight. This is automatically\\n    done for 4D multi-echo input and therefore not required.\\n- correctglobal=false: If true corrects global n2π offsets.\\n- individual=false: If true perform individual unwrapping of echos.\\n    Type ?unwrap_individual for more information\\n- template=1: echo that is spatially unwrapped (if individual is false)\\n- maxseeds=1: higher values allow more seperate regions\\n- merge_regions=false: spatially merge neighboring regions after unwrapping\\n- correct_regions=false: bring each regions median closest to 0 by adding n2π\\n- wrap_addition=0: [0;π], allows 'linear unwrapping', neighbors can have more\\n    (π+wrapaddition) phase difference\\n- `temporaluncertainunwrapping=0.5: Useful for veins. It uses spatial\\n    unwrapping on voxels that have high uncertainty values after temporal\\n    unwrapping. A higher quality threshold re-unwraps more voxels spatially.\\n\\n# Examples\\n``julia-repl\\njulia> using MriResearchTools\\njulia> phase = readphase(\\\"phase3echo.nii\\\")\\njulia> unwrapped = unwrap(phase; TEs=[1,2,3])\\njulia> savenii(unwrapped, \\\"unwrapped.nii\\\"; header=header(phase))\\n\\n\"), unwrap(wrapped::AbstractArray; keyargs...)\n\n\nROMEO unwrapping for 3D and 4D data.\n\n### Optional keyword arguments:\n\n  * `TEs`: Required for 4D data. The echo times for multi-echo data. In the case of single-echo    application with phase and the phase2 as a tuple (eg. (5, 10) or [5, 10]).\n  * `weights`: Options are [`:romeo`] | `:romeo2` | `:romeo3` | `:bestpath`.\n  * `mag`: The magnitude is used to improve the unwrapping-path.\n  * `mask`: Unwrapping is only performed inside the mask.\n  * `phase2`: A second reference phase image (possibly with different echo time).   It is used for calculating the phasecoherence weight. This is automatically   done for 4D multi-echo input and therefore not required.\n  * `correctglobal=false`: If `true` corrects global n2π offsets.\n  * `individual=false`: If `true` perform individual unwrapping of echos.   Type `?unwrap_individual` for more information\n  * `template=1`: echo that is spatially unwrapped (if `individual` is `false`)\n  * `maxseeds=1`: higher values allow more seperate regions\n  * `merge_regions=false`: spatially merge neighboring regions after unwrapping\n  * `correct_regions=false`: bring each regions median closest to 0 by adding n2π\n  * `wrap_addition=0`: [0;π], allows 'linear unwrapping', neighbors can have more   (π+wrap_addition) phase difference\n  * `temporal_uncertain_unwrapping=0.5`: Useful for veins. It uses spatial   unwrapping on voxels that have high uncertainty values after temporal   unwrapping. A higher quality threshold re-unwraps more voxels spatially.\n\n# Examples\n\n\njulia-repl julia> using MriResearchTools julia> phase = readphase(\"phase_3echo.nii\") julia> unwrapped = unwrap(phase; TEs=[1,2,3]) julia> savenii(unwrapped, \"unwrapped.nii\"; header=header(phase)) ``` , Dict{Symbol, Any}(:typesig => Union{}, :module => ROMEO, :linenumber => 28, :binding => ROMEO.unwrap!, :path => \"/home/runner/work/ROMEO.jl/ROMEO.jl/src/unwrapping.jl\"))\n\n\n\n\n\n","category":"function"},{"location":"#ROMEO.voxelquality-Tuple{Any}","page":"Home","title":"ROMEO.voxelquality","text":"voxelquality(phase::AbstractArray; keyargs...)\n\nCalculates a quality for each voxel. The voxel quality can be used to create a mask. The quality range is [0;1]\n\nExamples\n\njulia> qmap = voxelquality(phase_3echo; TEs=[1,2,3]);\njulia> mask = robustmask(qmap);\n\nTakes the same inputs as romeo/unwrap: Base.Docs.DocStr(svec(\"    unwrap(wrapped::AbstractArray; keyargs...)\\n\\nROMEO unwrapping for 3D and 4D data.\\n\\n###  Optional keyword arguments:\\n\\n- TEs: Required for 4D data. The echo times for multi-echo data. In the case of single-echo \\n    application with phase and the phase2 as a tuple (eg. (5, 10) or [5, 10]).\\n- weights: Options are [:romeo] | :romeo2 | :romeo3 | :bestpath.\\n- mag: The magnitude is used to improve the unwrapping-path.\\n- mask: Unwrapping is only performed inside the mask.\\n- phase2: A second reference phase image (possibly with different echo time).\\n    It is used for calculating the phasecoherence weight. This is automatically\\n    done for 4D multi-echo input and therefore not required.\\n- correctglobal=false: If true corrects global n2π offsets.\\n- individual=false: If true perform individual unwrapping of echos.\\n    Type ?unwrap_individual for more information\\n- template=1: echo that is spatially unwrapped (if individual is false)\\n- maxseeds=1: higher values allow more seperate regions\\n- merge_regions=false: spatially merge neighboring regions after unwrapping\\n- correct_regions=false: bring each regions median closest to 0 by adding n2π\\n- wrap_addition=0: [0;π], allows 'linear unwrapping', neighbors can have more\\n    (π+wrapaddition) phase difference\\n- `temporaluncertainunwrapping=0.5: Useful for veins. It uses spatial\\n    unwrapping on voxels that have high uncertainty values after temporal\\n    unwrapping. A higher quality threshold re-unwraps more voxels spatially.\\n\\n# Examples\\n``julia-repl\\njulia> using MriResearchTools\\njulia> phase = readphase(\\\"phase3echo.nii\\\")\\njulia> unwrapped = unwrap(phase; TEs=[1,2,3])\\njulia> savenii(unwrapped, \\\"unwrapped.nii\\\"; header=header(phase))\\n\\n\"), unwrap(wrapped::AbstractArray; keyargs...)\n\n\nROMEO unwrapping for 3D and 4D data.\n\n### Optional keyword arguments:\n\n  * `TEs`: Required for 4D data. The echo times for multi-echo data. In the case of single-echo    application with phase and the phase2 as a tuple (eg. (5, 10) or [5, 10]).\n  * `weights`: Options are [`:romeo`] | `:romeo2` | `:romeo3` | `:bestpath`.\n  * `mag`: The magnitude is used to improve the unwrapping-path.\n  * `mask`: Unwrapping is only performed inside the mask.\n  * `phase2`: A second reference phase image (possibly with different echo time).   It is used for calculating the phasecoherence weight. This is automatically   done for 4D multi-echo input and therefore not required.\n  * `correctglobal=false`: If `true` corrects global n2π offsets.\n  * `individual=false`: If `true` perform individual unwrapping of echos.   Type `?unwrap_individual` for more information\n  * `template=1`: echo that is spatially unwrapped (if `individual` is `false`)\n  * `maxseeds=1`: higher values allow more seperate regions\n  * `merge_regions=false`: spatially merge neighboring regions after unwrapping\n  * `correct_regions=false`: bring each regions median closest to 0 by adding n2π\n  * `wrap_addition=0`: [0;π], allows 'linear unwrapping', neighbors can have more   (π+wrap_addition) phase difference\n  * `temporal_uncertain_unwrapping=0.5`: Useful for veins. It uses spatial   unwrapping on voxels that have high uncertainty values after temporal   unwrapping. A higher quality threshold re-unwraps more voxels spatially.\n\n# Examples\n\n\njulia-repl julia> using MriResearchTools julia> phase = readphase(\"phase_3echo.nii\") julia> unwrapped = unwrap(phase; TEs=[1,2,3]) julia> savenii(unwrapped, \"unwrapped.nii\"; header=header(phase)) ``` , Dict{Symbol, Any}(:typesig => Union{}, :module => ROMEO, :linenumber => 28, :binding => ROMEO.unwrap!, :path => \"/home/runner/work/ROMEO.jl/ROMEO.jl/src/unwrapping.jl\"))\n\nSee also unwrap\n\n\n\n\n\n","category":"method"}]
}
