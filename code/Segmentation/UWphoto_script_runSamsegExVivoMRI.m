UWphoto_function_checkBashEnvironment();

flair_dir = '/home/henry/Documents/Brain/UWphoto/FLAIR_Scan_Data';
dlist_mris = dir(fullfile(flair_dir,'*.rotated.mgz'));


codedir = '/home/henry/Documents/Brain/UWphoto/code';

templatedir = '/home/henry/Documents/Brain/UWphoto/code/data';
registerdir = fullfile(templatedir,'MRItemplates');

segmnteddir = fullfile(flair_dir,'SAMSEG/');

if ~exist(segmnteddir,'dir')
    mkdir(segmnteddir)
end


atlasDir = ['/home/henry/Documents/Brain/Freesurfer/freesurfer/average/',...
    'samseg/20Subjects_smoothing2_down2_smoothingForAffine2'];

BSCBlistfile='/home/henry/Documents/Brain/UWphoto/code/data/BS_CB_cases.mat';

load('/home/henry/Documents/Brain/UWphoto/code/data/BS_CB_cases.mat',...
    'BS_CB_cases');

exvivoGMM = fullfile(templatedir,'sharedGMMParameters.exvivoMRI.txt');
noBSCBGMM = fullfile(templatedir,'sharedGMMParameters.NoBSorCB.exvivoMRI.txt');


%% run registration
problemsfile = '/home/henry/Documents/Brain/UWphoto/code/data/remainingsegmentations.mat';

if exist(problemsfile,'file')
    load(problemsfile,'problems_flag')
else
    problems_flag = true(size(dlist_mris));
end

dlist_mris=dlist_mris(problems_flag);


for il=1:length(dlist_mris)
   
    
    %% settup filenames for case
    [~,basename,~] = fileparts(dlist_mris(il).name);
    [~,caseID,~]   = fileparts(basename);
    
    MRI_mgz_version = fullfile(flair_dir,dlist_mris(il).name);
    movedhder = fullfile(registerdir,[caseID,'_template.mgz']);
    
    ptntdir = fullfile(segmnteddir,caseID);
    
    tempdir = fullfile(ptntdir,'temp');
    
    if ~exist(ptntdir,'dir')
        mkdir(ptntdir)
    end
    
    if ~endsWith(caseID,BS_CB_cases)
        rmdir(tempdir,'s')
    end
    
    if ~exist(tempdir,'dir')
        mkdir(tempdir)
        
        copyfile([atlasDir '/atlas_level1.txt.gz' ],tempdir);
        copyfile([atlasDir '/atlas_level2.txt.gz' ],tempdir);
        copyfile([atlasDir '/compressionLookupTable.txt' ],tempdir);
        copyfile([atlasDir '/modifiedFreeSurferColorLUT.txt' ],tempdir);
        
        if endsWith(caseID,BS_CB_cases)
            copyfile(exvivoGMM,fullfile(tempdir, 'sharedGMMParameters.txt'));
        else
            copyfile(noBSCBGMM,fullfile(tempdir, 'sharedGMMParameters.txt'));
        end 
           
    end
    
    
    SAMSEG_DATA_DIR=tempdir;
    setenv('SAMSEG_DATA_DIR',SAMSEG_DATA_DIR)
    
    
    
    %% 
    
    cmd = ['sh ',codedir,'/Segmentation/run_photoseg ',...
        '-o ',ptntdir,' -i ',MRI_mgz_version,' ',...
        '--threads 32 --template ',movedhder];
    
    dlist_segmentations = dir(fullfile(ptntdir,'*crispSegmentation*'));
    
    if isempty(dlist_segmentations) || ~endsWith(caseID,BS_CB_cases)
        try
            fprintf(['Segmenting case ',caseID,'\n'])
            [segstatus,segresult]=system(cmd);
            if segstatus
                warning(['problem with segmentation of ',caseID])
            else
                problems_flag(il)=false;
                fprintf(segresult(end-35:end))
            end
        catch
            warning(['problem with segmentation of ',caseID])
        end
    else
        problems_flag(il)=false;
        fprintf(['Skipping case ',caseID,'\n'])
    end
end