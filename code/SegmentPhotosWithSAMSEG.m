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
function SegmentPhotosWithSAMSEG(inputVol,inputVolMask,inputWarpedRef,outputDir,nthread)

% clear
% clc
% 
% inputVol='/autofs/cluster/vive/UW_photo_recon/recons/outputsHardAtlasBin/18-1132.recon.mgz';
% inputVolMask='/autofs/cluster/vive/UW_photo_recon/recons/outputsHardAtlasBin/18-1132.mask.mgz';
% inputWarpedRef='/autofs/cluster/vive/UW_photo_recon/recons/outputsHardAtlasBin/18-1132.warped_ref.mgz';
% outputDir='/autofs/cluster/vive/UW_photo_recon/recons/outputsHardAtlasBin/18-1132-samseg/';

ROB_REG='/usr/local/freesurfer/stable6_0_0/bin/mri_robust_register';
EXVIVO_SAMSEG_FS_DIR='/autofs/space/panamint_005/users/iglesias/fsdevsamsegexvivo/';

% STUFF YOU PROBABLY SHOULDN'T TOUCH
BG_NOISE_MAX=3; % maximum level of noise that we introduce in background to make samseg work
SAMSEG_CHANNELS=1;  % 1 for grayscale, 2 for PCA with 2 components, 3 for full RGB. Best setting is 1.

%%%%%%%%%%%

if nargin < 5
    nthread = 4;
end

%%%%%%%%%%%

addpath functions

% Output and temp directories
if exist(outputDir,'dir')==0
    mkdir(outputDir);
end

tempdir=[outputDir '/temp/'];
if exist(tempdir,'dir')==0
    mkdir(tempdir);
end

disp('Reading in input volume and matching histogram of photos ...');
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
MRIwrite(mriCorr,[outputDir '/corrected.nii.gz']);


disp('Preparing input scans...');
inputs=[];
if SAMSEG_CHANNELS==1 % simply rgb->gray
    mri.vol=mean(mriCorr.vol,4);
    aux=mri.vol(mri.vol==0); aux=aux+BG_NOISE_MAX*rand(size(aux)); mri.vol(mri.vol==0)=aux;
    inputs{1}=[outputDir '/corrected_gray.nii.gz'];
    MRIwrite(mri,inputs{1});
    
elseif SAMSEG_CHANNELS==2 % reduce dimensionality with PCA
    X=zeros([numel(mriCorr.vol)/3 3]);
    for c=1:3
        aux=mriCorr.vol(:,:,:,c);
        X(:,c)=aux(:);
    end
    C=cov(X);
    [V,D]=eig(C);
    [~,idx]=sort(diag(D),'descend');
    idx=idx(1:2);
    P=V(:,idx);
    mu=mean(X)';
    b=P'*(X'-repmat(mu,[1 size(X,1)]));
    b=b+0.1*randn(size(b));
    
    aux=mriCorr;
    aux.vol=aux.vol(:,:,:,1);
    aux.vol(:)=b(1,:);
    aux.vol=aux.vol-min(aux.vol(:))+1e-3;
    inputs{1}=[outputDir '/corrected_c1.nii.gz'];
    MRIwrite(aux,inputs{1});
    
    aux=mriCorr;
    aux.vol=aux.vol(:,:,:,1);
    aux.vol(:)=b(2,:);
    aux.vol=aux.vol-min(aux.vol(:))+1e-3;
    inputs{2}=[outputDir '/corrected_c2.nii.gz'];
    MRIwrite(aux,inputs{2});
    
else  % run on RGB
    mri.vol=mriCorr.vol(:,:,:,1);
    aux=mri.vol(mri.vol==0); aux=aux+BG_NOISE_MAX*rand(size(aux)); mri.vol(mri.vol==0)=aux;
    inputs{1}=[outputDir '/corrected_red.nii.gz'];
    MRIwrite(mri,inputs{1});
    mri.vol=mriCorr.vol(:,:,:,2);
    aux=mri.vol(mri.vol==0); aux=aux+BG_NOISE_MAX*rand(size(aux)); mri.vol(mri.vol==0)=aux;
    inputs{2}=[outputDir '/corrected_green.nii.gz'];
    MRIwrite(mri,inputs{2});
    mri.vol=mriCorr.vol(:,:,:,3);
    aux=mri.vol(mri.vol==0); aux=aux+BG_NOISE_MAX*rand(size(aux)); mri.vol(mri.vol==0)=aux;
    inputs{3}=[outputDir '/corrected_blue.nii.gz'];
    MRIwrite(mri,inputs{3});
end


disp('Initializing affine registration');
cmd=[ROB_REG ' --leastsquares --dst ' inputWarpedRef ' --mov ./data/template.nii  --nosym --mapmovhdr ' tempdir '/rigid.mgz --lta ' tempdir '/rigid.lta'];
a=system([cmd ' >/dev/null']);
if a, error('Error in rigid registration of template'); end
cmd=[ROB_REG ' --leastsquares --affine --dst ' inputWarpedRef ' --mov ' tempdir '/rigid.mgz  --nosym --mapmovhdr ' tempdir '/template.mgz --lta ' tempdir '/affine.lta'];
a=system([cmd ' >/dev/null']);
if a, error('Error in affine registration of template'); end

% copy rest of atlas files
atlasDir=[EXVIVO_SAMSEG_FS_DIR '/average/samseg/20Subjects_smoothing2_down2_smoothingForAffine2/'];
copyfile([atlasDir '/atlas_level1.txt.gz' ],tempdir);
copyfile([atlasDir '/atlas_level2.txt.gz' ],tempdir);
copyfile([atlasDir '/compressionLookupTable.txt' ],tempdir);
copyfile([atlasDir '/modifiedFreeSurferColorLUT.txt' ],tempdir);
copyfile('./data/sharedGMMParameters.NoBSorCB.photo.txt',[tempdir '/sharedGMMParameters.txt']);

% samseg call
disp('Running samseg')
cmd=['setenv FREESURFER_HOME ' EXVIVO_SAMSEG_FS_DIR ' && source $FREESURFER_HOME/SetUpFreeSurfer.csh &&  samseg '];
for j=1:length(inputs)
    cmd=[cmd ' --i ' inputs{j}];
end
cmd=[cmd ' --o ' outputDir];
cmd=[cmd ' --threads ' num2str(nthread)];
cmd=[cmd ' --ssdd ' tempdir];
cmd=[cmd ' --save-posteriors'];
cmd=[cmd ' --no-save-warp'];
disp(cmd);
system(cmd);

disp('All done! Cleaning up...');
system(['rm -rf ' tempdir]);







