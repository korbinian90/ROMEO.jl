include("UnwrappingExecutable.jl")
# 7T EPI patients
d = "/media/barbara/hdd2/DATA/phase_unwrapping_comparison/7T/19560318BLBL_201501201000_analysis/DCwithMarques_SingleGEFM_mag_20/"
@time unwrapping_main([d*"PhaseCombined_withGE7.nii", "-o", d*"test_romeo2.nii", "-m", d*"MagCombSS_withGE7.nii", "-t", "ones(57)", "-k", d*"MagCombSS_withGE7_masked_mask.nii"])
