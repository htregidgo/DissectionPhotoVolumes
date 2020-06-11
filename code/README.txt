Code to process photos from UW

FUNCTIONS / SCRIPTS (please see further info inside the individual functions / scripts)

UWphoto_startup.m:
Checks that freesurfer has been initialised and sets up some environment 
variables.

UWphoto_script_SoftReconstructions.m:
A script to select source photos and call the soft reconstruction function.

ReconPhotoVolume_joint_multires.m: 
Function to reconstruct a photo volume using a probabilistic atlas as reference

ReconPhotoVolume_joint_hard_multires.m:
Function to reconstruct a photo volume using a hard binary segmentation as reference

makeMovieFramesHardMask.m: 
Prepares volumes and surfaces for making a registration movie using Freeview.
Works with the output of ReconPhotoVolume_joint_hard_multires.
It may be good in the future to write an equivalent script for ReconPhotoVolume_joint_multires.m?



OLD SEGMENTATION CODE

SegmentPhotosWithSAMSEG.m: 
Segments a photo reconstruction with a modified version of SAMSEG.

SegmentMRIsWithSAMSEG.m: 
Segments an ex vivo MRIs with a modified version of SAMSEG. 


DIRECTORIES

functions: 
A bunch of functions required by the different functions / scripts in this directory,
including the cost functions for the optimization in ReconPhotoVolume_joint*.m

data:
Some SAMSEG-compatible sharedGMMparameter files, and a binarized version of the SAMSEG
template including only the cerebrum, which is sometimes useful to initialize the segmentation
(by registering this template with other tool or even manually, and bypassing the affine
registration in the python SAMSEG code)






PerspectiveNumberingBatch.*.m: different versions of software to manually annotate slices

MakeAffineAtlases.m: makes affine atlases with / without cerebellum and brainstem


SegmentPhotosSAMSEGscript.m: script that uses the function above to segment all cases in a directory

SegmentMRIsWithSAMSEG.m: script to run SAMSEG on the MRIs

CompareVolumes.m: compares volumes from MRI and photos

%%%%%%%%%%%%%%%%%%%%%%%%%
Old stuff...

ReconPhotsVolume_average_onlyAtlas.m: version that uses average brain (register slices only to atlas)
ReconPhotoVolume_average.m: version that uses average brain (registers to nrighboring slices too). 
   Has a hard flag in line 178 to choose between intensities or only masks. It does something a bit crappy,
   whic is linearly registering to the average of 3 channels. Need to code this properly, and also, penalyze
   determinants so they don't go to to zero! One more thing: grouping affine transforms of slices
   within the same photo would be good!

ReconPhotoVolume.m: prototype of code to reconstruct volume from MRI + photographs
ReconPhotoVolume_onlyGray.m: version that uses only grayscale

ReconPhotoVolume_binary.m: version that uses only binary mask

ReconPhotoVolume_joint.m: uses average brain to reconstruct photos 
ReconPhotoVolume_joint_hard.m: tversion with hard labels

ReconPhotoScript.m: uses ReconPhotoVolume_binary to reconstruct all samples in a directory




