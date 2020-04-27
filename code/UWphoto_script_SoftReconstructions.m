%UWPHOTO_SCRIPT_SOFTRECONSTRUCTIONS Gui script to reconstruct a directory
%
%   This script should be run after using UWPHOTO_startup.m to ensure that
%   the correct paths have been added.
%
%   The script will open UI windows to ask for two directories and possibly 
%   one file. These are:
%       -Top level photo directory 
%       a directory containing all the cases to be reconstructed. The file
%       structure should be 
%           <topLevelDir>/<CASE_ID>/"<CASE_ID> MATLAB"
%       these matlab directories should contain matching .tiff files and
%       .mat files for the perspective corrected slice photos and slice
%       masks respectively. Optionally a file 
%           <topLevelDir>/<CASE_ID>/slice_order.mat 
%       may be included containing a vector "slice_order" to define the
%       correct order of slices where the masks have been incorrectly
%       numbered.
%
%       -Top level results directory
%       Destination for reconstructed files. A new filestructure will be
%       created in this location, with folders for individual cases.
%
%       -Reference File 
%       This is the location of the probabilistic atlas file used for
%       registration. By default this is in 
%           <UWphoto top dir>/prob_atlases/onlyCerebrum.smoothed.nii.gz
%
%
%   The code is currently set up for source photos with pixel sizes of
%   0.1mm and slice thickness of 4mm. The resulting volume have anisotropic
%   voxel dimensions of 0.5mm x 0.5mm x 4mm.
%
%   The source photo dimensions can be modified by examining PHOTO_RES and
%   SLICE_THICKNESS in the code. 
%
% Henry Tregidgo h.tregidgo@ucl.ac.uk



%% initialise required variables

% flag to force recalculation
forceFlag = false;

% freesurfer home
FREESURFER_HOME = getenv('FREESURFER_HOME');

if isempty(FREESURFER_HOME)
    
    FREESURFER_HOME='/home/henry/Documents/Brain/Freesurfer/Stable/freesurfer';
    
    if ~exist(FREESURFER_HOME,'dir')
        FREESURFER_HOME=uigetdir(pwd,'Freesurfer directory');
        
        if isequal(0,FREESURFER_HOME)
            error('freesurfer directory not set')
        end
    end
    
    setenv('FREESURFER_HOME',FREESURFER_HOME)
end


% resolutions and scheduling
PHOTO_RES = 0.1;
SLICE_THICKNESS = 4;
TARGET_RES = [4 2 1 0.5]; 
scheduleITs = [50 40 30; 25 20 15; 12 10 8; 6 5 5]; 


% define I/O

topPhotoDir = '';
if ~exist(topPhotoDir,'dir')
    topPhotoDir=uigetdir(pwd,'Top level photo directory');
    
    if isequal(0,topPhotoDir)
        error('top level photo directory not set')
    end
end

dlist_cases = dir(fullfile(topPhotoDir,'*/*MATLAB')); % list cases

topOutDir = '';
if ~exist(topOutDir,'dir')
    topOutDir=uigetdir(topPhotoDir,'Top level results directory');
    
    if isequal(0,topOutDir)
        error('top level results directory not set')
    end
end

inputREFERENCE = '';
if ~exist(inputREFERENCE,'file')
    PHOTO_RECON_HOME = getenv('PHOTO_RECON_HOME');
    
    if ~isempty(PHOTO_RECON_HOME)
        
        HemiType = questdlg('Which hemisphere is present in the images?',...
            'Hemisphere selection','Left','Right','Both','Both');
        
        switch HemiType
            case 'Right'
                inputREFERENCE = fullfile(PHOTO_RECON_HOME,'prob_atlases',...
                    'right_onlyCerebrum.smoothed.nii.gz');
            case 'Left'
                inputREFERENCE = fullfile(PHOTO_RECON_HOME,'prob_atlases',...
                    'left_onlyCerebrum.smoothed.nii.gz');
            otherwise
                inputREFERENCE = fullfile(PHOTO_RECON_HOME,'prob_atlases',...
                    'onlyCerebrum.smoothed.nii.gz');
        end
    end
    if ~exist(inputREFERENCE,'file')
        
        [inputREFERENCEname,inputREFERENCEpath]=uigetfile('*.nii.gz','Select Reference file');
        inputREFERENCE = fullfile(inputREFERENCEpath,inputREFERENCEname);
        
        if isequal(0,inputREFERENCE)
            error('soft reference not set')
        end
    end
end

viewpoint = questdlg('Is the anterior or posterior side of the slice showing?',...
    'Viewpoint selection','Anterior','Posterior','Posterior');


problems_flag = true(size(dlist_cases));

%% iterate through cases

for il = 1:(length(dlist_cases))
    
    [~,caseID,~] = fileparts(dlist_cases(il).folder);
    inputPhotoDir = fullfile(dlist_cases(il).folder,dlist_cases(il).name);
    
    outputDir = fullfile(topOutDir,caseID,'soft');
    
    if ~exist(outputDir,'dir')
        mkdir(outputDir)
    end
    
    outputVol = fullfile(outputDir,[caseID,'_soft.mgz']);
    outputVolMask = fullfile(outputDir,[caseID,'_soft_mask.mgz']);
    outputWarpedRef = fullfile(outputDir,[caseID,'_soft_regatlas.mgz']);
    outputMat = fullfile(outputDir,[caseID,'_soft_history.mat']);
    
    try
        if forceFlag || ~exist(outputVol,'file')
            ReconPhotoVolume_joint_multires(inputPhotoDir,inputREFERENCE,outputVol,...
                outputVolMask,outputWarpedRef,outputMat,PHOTO_RES,SLICE_THICKNESS,...
                TARGET_RES,scheduleITs,'APswitch',viewpoint)
            
        end
        
        problems_flag(il)=false;
    catch
        warning(['Problem with running case ',caseID])
    end
    
end

%% check for errors

if any(problems_flag)
    
    incomplete_volumes = {dlist_cases(problems_flag).name};
    
    
    fprintf('\n========================================\n')
    fprintf('reconstructions complete with errors:')
    for il=1:length(incomplete_volumes)
        fprintf(['\n',incomplete_volumes{il}])
    end
    fprintf('\n========================================\n')
    
else
    
    fprintf('\n========================================\n')
    fprintf('reconstructions complete with no errors!')
    fprintf('\n========================================\n')
    
end