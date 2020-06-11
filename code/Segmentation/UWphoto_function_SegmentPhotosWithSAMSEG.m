% Function that segments a photo reconstruction with a modified version of SAMSEG
%
% PARAMETERS
%
% inputVol: a volume of photographs reconstructed with ReconPhotoVolume_joint*.m (either soft or hard)
% inputVolMask: the corresponding reconstructed mask
% inputWarpedRef: the corresponding deformed atlas (hard or soft, it doesn't matter!)
% outputDir: the output directory
% nthread: number of threads for SAMSEG
%
% VERY IMPORTANT NOTES:
% 1. Please set ROB_REG below to $FREESURFER_HOME/bin/mri_robust_register 
% 2. Please find a better alternative to the incredibly ugly code I have
% below (we can discuss over Skype). Essentially, right now I have a copy
% of freesurfer under EXVIVO_SAMSEG_FS_DIR (defined below) where I switched 
% off bias the field correction in samsegment.py (find estimateBiasField in
% the module getOptimizationOptions and set it to False) and where I bypass
% the affine registration in python/scripts/run_samseg, by commenting out 
% call to registerAtlas (and the immediate assignment to
% transformedTemplateFileName), and instead assigning transformedTemplateFileName
% directly to:
% transformedTemplateFileName = os.path.join(atlasDir, 'template.mgz')
% Of course, this requires that the template.mgz in the atlas directory is 
% already registered to the target volume to segment. But this is
% essentially what we do in this function: preparing an atlas directory
% with a registered template.mgz, and our customized
% sharedGaussianParameters.txt file, before calling our modified samseg.
%
function[segstatus,segresult]= UWphoto_function_SegmentPhotosWithSAMSEG(...
    caseID,inputVol,inputVolMask,inputWarpedRef,outputDir,nthread,varargin)

% clear
% clc
% 
% inputVol='/autofs/cluster/vive/UW_photo_recon/recons/outputsHardAtlasBin/18-1132.recon.mgz';
% inputVolMask='/autofs/cluster/vive/UW_photo_recon/recons/outputsHardAtlasBin/18-1132.mask.mgz';
% inputWarpedRef='/autofs/cluster/vive/UW_photo_recon/recons/outputsHardAtlasBin/18-1132.warped_ref.mgz';
% outputDir='/autofs/cluster/vive/UW_photo_recon/recons/outputsHardAtlasBin/18-1132-samseg/';

% ROB_REG='/usr/local/freesurfer/stable6_0_0/bin/mri_robust_register';
% EXVIVO_SAMSEG_FS_DIR='/autofs/space/panamint_005/users/iglesias/fsdevsamsegexvivo/';

%% get environment variables

FREESURFER_HOME = getenv('FREESURFER_HOME');
if isempty(FREESURFER_HOME)
    UWphoto_function_checkBashEnvironment();
    FREESURFER_HOME = getenv('FREESURFER_HOME');
end

PHOTO_RECON_HOME = getenv('PHOTO_RECON_HOME');
if isempty(PHOTO_RECON_HOME)
    error(['please set environment variable PHOTO_RECON_HOME.\n',...
        'run UWphoto_startup'])
end

codedir = fullfile(PHOTO_RECON_HOME,'code');

dlist_robreg = dir(fullfile(FREESURFER_HOME,'bin/mri_robust_register'));

if isempty(dlist_robreg)
    error('need to find robust register, try setting FREESURFER_HOME')
end

ROB_REG = fullfile(dlist_robreg.folder,dlist_robreg.name);

%% parse inputs


p=inputParser;

addRequired(p,'caseID',@(x) ischar(x))
addRequired(p,'inputVol',@(x) ischar(x) && exist(x,'file') &&...
    endsWith(x,{'.mgz','.nii','.nii.gz'}))
addRequired(p,'inputVolMask',@(x) ischar(x) && exist(x,'file') &&...
    endsWith(x,{'.mgz','.nii','.nii.gz'}))
addRequired(p,'inputWarpedRef',@(x) ischar(x) && exist(x,'file') &&...
    endsWith(x,{'.mgz','.nii','.nii.gz'}))
addRequired(p,'outputDir',@(x)  ischar(x))
addOptional(p,'nthread',maxNumCompThreads,@(x) isnumeric(x) && isreal(x) && ...
    x==floor(x) && x>0 && x<=2*maxNumCompThreads)
addParameter(p,'Hemisphere','Whole', @(x) any(strcmpi(x,...
    {'Whole','Left','Right'})))


parse(p,caseID,inputVol,inputVolMask,inputWarpedRef,outputDir,nthread,...
    varargin{:})

Hemi_switch = find(strcmpi(p.Results.Hemisphere,{'Whole','Left','Right'}));

%%


% '/home/henry/Documents/Brain/Freesurfer/freesurfer/bin/mri_robust_register';

% STUFF YOU PROBABLY SHOULDN'T TOUCH
BG_NOISE_MAX=3; % maximum level of noise that we introduce in background to make samseg work
% SAMSEG_CHANNELS=1;  % 1 for grayscale, 2 for PCA with 2 components, 3 for full RGB. Best setting is 1.

%%%%%%%%%%%

if nargin < 5
    nthread = 4;
end

%%%%%%%%%%%

% addpath functions


% Output and temp directories
if exist(outputDir,'dir')==0
    mkdir(outputDir);
end

tempdir=[outputDir '/temp/'];
if exist(tempdir,'dir')==0
    mkdir(tempdir);
end

disp('Reading in input volume and outputing greyscale file ...');
mri=MRIread(inputVol);
mriMask=MRIread(inputVolMask);
mri.vol(repmat(mriMask.vol<.5,[1 1 1 3]))=0;
mriCorr=mri;
mid=round(size(mri.vol,3)/2);
for c=1:3
    REF=uint8(mri.vol(:,:,mid,c));
    REF(mriMask.vol(:,:,mid)<0.5)=0;
    for z=[1:mid-1 mid+1:size(mri.vol,3)]
        I=uint8(mri.vol(:,:,z,c));
        I(mriMask.vol(:,:,z)<.5)=0;
        Imatched = matchHistoMasked(I, REF);
        mriCorr.vol(:,:,z,c)=Imatched;
    end
end
MRIwrite(mriCorr,fullfile(outputDir,[caseID, '_histogramcorrected.nii.gz']));

%% output grayscale scans for segmentation
disp('Preparing input scans...');
inputs=[];

se = strel('sphere',5);
dltdmask = imdilate(mriMask.vol>0.5,se);
% clsdmask = imclose(mriMask.vol>0.5,se);

mri_uncorrected = mri;
mri_uncorrected.vol = mean(mri.vol,4);
aux=mri_uncorrected.vol(mri_uncorrected.vol==0 & dltdmask);

% logmax=max(log(mri_uncorrected.vol(mri_uncorrected.vol(:)>0)));
% aux=3+logmax+(logmax/10)*randn(size(aux));
% aux=exp(aux);

aux=aux+BG_NOISE_MAX*rand(size(aux)); 

mri_uncorrected.vol(mri_uncorrected.vol==0 & dltdmask)=aux;
MRIwrite(mri_uncorrected,fullfile(outputDir,...
    [caseID,'_uncorrected_gray.nii.gz']));

inputs{1}=fullfile(outputDir,[caseID,'_uncorrected_gray.nii.gz']);


mri.vol=mean(mriCorr.vol,4);
aux=mri.vol(mri.vol==0 & dltdmask); 

aux=aux+BG_NOISE_MAX*rand(size(aux)); 

mri.vol(mri.vol==0 & dltdmask)=aux;

MRIwrite(mri,fullfile(outputDir,[caseID,'_histogramcorrected_gray.nii.gz']));
    

%% affine registration of template

disp('Initializing affine registration');

switch Hemi_switch
    case 2
        unregtemplate = fullfile(codedir,'data/template_left.nii');
    case 3
        unregtemplate = fullfile(codedir,'data/template_right.nii');
    otherwise
        unregtemplate = fullfile(codedir,'data/template.nii');
end

if ~exist(unregtemplate,'file')
    error(['unregistered photo volume template required in \n',unregtemplate])
end

% run rigid alignment
cmd=[ROB_REG, ' --leastsquares --dst ', inputWarpedRef,...
    ' --mov ',unregtemplate,' --nosym --mapmovhdr ', tempdir,...
    '/rigid.mgz --lta ' tempdir '/rigid.lta'];

a=system([cmd ' >/dev/null']);

% check for errors in rigid registration
if a, error('Error in rigid registration of template'); end

% run affine registration of mask
cmd=[ROB_REG, ' --leastsquares --affine --dst ', inputWarpedRef,...
    ' --mov ', tempdir, '/rigid.mgz  --nosym --mapmovhdr ', tempdir,...
    '/template.mgz --lta ' tempdir '/affine.lta'];

a=system([cmd ' >/dev/null']);

% check for errors in affine registration
if a, error('Error in affine registration of template'); end

%% copy rest of atlas files
atlasDir=[FREESURFER_HOME '/average/samseg/20Subjects_smoothing2_down2_smoothingForAffine2/'];
copyfile([atlasDir '/atlas_level1.txt.gz' ],tempdir);
copyfile([atlasDir '/atlas_level2.txt.gz' ],tempdir);
copyfile([atlasDir '/compressionLookupTable.txt' ],tempdir);
copyfile([atlasDir '/modifiedFreeSurferColorLUT.txt' ],tempdir);

switch Hemi_switch
    case 2
        copyfile(fullfile(codedir,'data/sharedGMMParameters.LeftNoBSorCB.photo.txt'),...
            [tempdir '/sharedGMMParameters.txt']);
    case 3
        copyfile(fullfile(codedir,'data/sharedGMMParameters.RightNoBSorCB.photo.txt'),...
            [tempdir '/sharedGMMParameters.txt']);
    otherwise
        copyfile(fullfile(codedir,'data/sharedGMMParameters.NoBSorCB.photo.txt'),...
            [tempdir '/sharedGMMParameters.txt']);
end
%% samseg call

disp('Running samseg')

SAMSEG_DATA_DIR=tempdir;
setenv('SAMSEG_DATA_DIR',SAMSEG_DATA_DIR)

movedhder = fullfile(tempdir,'template.mgz');

cmd = ['sh ',codedir,'/Segmentation/run_photoseg ',...
    '-o ',outputDir,' -i ',inputs{1},' ',...
    '--threads ',num2str(nthread),' --photo-seg ',movedhder];


fprintf(['Segmenting case ',caseID,'\n'])
[segstatus,segresult]=system(cmd);


%% correct for brightness
dlist_scale = dir(fullfile(outputDir,'*scaling-factor.txt'));
dlist_BF = dir(fullfile(outputDir,'*_biasField*'));

scal_factr_path = fullfile(dlist_scale.folder,dlist_scale.name);
bias_field_path = fullfile(dlist_BF.folder,dlist_BF.name);

scale_factor = dlmread(scal_factr_path);

vol_raw=MRIread(inputVol);

bias_field=MRIread(bias_field_path);

BC_out_path = fullfile(outputDir,[caseID,'_brightnessCorrected.nii.gz']);

%% correct photo volume
vol_corrected=vol_raw;

for il=1:3
vol_corrected.vol(:,:,:,il)=squeeze(vol_raw.vol(:,:,:,il))./bias_field.vol...
    ./scale_factor;
end

MRIwrite( vol_corrected,BC_out_path);


end



%% local functions





