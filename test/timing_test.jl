using MRI

phasefile = raw"F:\MRI\scanner_nifti\Paper\SWI_paper_7T_volunteers\paper1\nifti\20\reform\Image.nii"
magfile = raw"F:\MRI\scanner_nifti\Paper\SWI_paper_7T_volunteers\paper1\nifti\19\reform\Image.nii"
phaseni = readphase(phasefile)
magni = readmag(magfile)

@time unwrap(phaseni; mag=magni)
