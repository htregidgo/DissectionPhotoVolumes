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
TARGET_RES = [4 2 1 0.5]; %0.1];%[4 2];
scheduleITs = [50 40 30; 25 20 15; 12 10 8; 6 5 5]; %; 2 2 2];%[45,30,20;15,10,10];

PHOTO_RECON_HOME=getenv('PHOTO_RECON_HOME');


% define I/O 

topPhotoDir = fullfile(PHOTO_RECON_HOME,'Photo_data_updated');
if ~exist(topPhotoDir,'dir')
    topPhotoDir=uigetdir(pwd,'Top level photo directory');
    
    if isequal(0,topPhotoDir)
        error('top level photo directory not set')
    end
end

dlist_cases = dir(fullfile(topPhotoDir,'*/*MATLAB')); % list cases

topOutDir = fullfile(PHOTO_RECON_HOME,'Results');
if ~exist(topOutDir,'dir')
    topOutDir=uigetdir(pwd,'Top level results directory');
    
    if isequal(0,topOutDir)
        error('top level results directory not set')
    end
end

inputREFERENCE = fullfile(PHOTO_RECON_HOME,'prob_atlases',...
    'onlyCerebrum.smoothed.nii.gz');
if ~exist(inputREFERENCE,'file')
    [inputREFERENCEname,inputREFERENCEpath]=uigetfile('*.nii.gz','Reference file');
    inputREFERENCE = fullfile(inputREFERENCEpath,inputREFERENCEname);
    
    if isequal(0,inputREFERENCE)
        error('top level results directory not set')
    end
end

manualSegsdir = fullfile(PHOTO_RECON_HOME,'ManualSegmentations');

problems_flag = true(size(dlist_cases));

%% iterate through cases

for il = 1:length(dlist_cases)

    [~,caseID,~] = fileparts(dlist_cases(il).folder);
    inputPhotoDir = fullfile(dlist_cases(il).folder,dlist_cases(il).name);
    
    dlist_labels = dir(fullfile(manualSegsdir,[caseID,'*Labels*']));
    dlist_masks = dir(fullfile(manualSegsdir,[caseID,'*Mask*']));
    
    if isempty(dlist_labels)
        continue
    end
    
    niftiMask = fullfile(dlist_masks.folder,dlist_masks.name);
    niftiLabel = fullfile(dlist_labels.folder,dlist_labels.name);
    
    outputDir = fullfile(topOutDir,caseID,'soft');
    
    if ~exist(outputDir,'dir')
        mkdir(outputDir)
    end
    
    outputLabel = fullfile(outputDir,[caseID,'_manualLabel.mgz']);
    outputVolMask = fullfile(outputDir,[caseID,'_soft_mask.mgz']);
    outputWarpedRef = fullfile(outputDir,[caseID,'_soft_regatlas.mgz']);
    paramMat = fullfile(outputDir,[caseID,'_soft_history.mat']);
     
    try
        if forceFlag || ~exist(outputLabel,'file')
            UWphoto_function_segFromManualLabel(inputPhotoDir,inputREFERENCE,outputLabel,...
                paramMat,PHOTO_RES,SLICE_THICKNESS,...
                TARGET_RES,niftiMask,niftiLabel)
            
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