function UWphoto_function_segFromManualLabel(inputPhotoDir,inputREFERENCE,outputlabel,...
    paramMat,PHOTO_RES,SLICE_THICKNESS,TARGET_RES,niftiMask,niftiLabel,recontype)
%UWPHOTO_FUNCTION_SEGFROMMANUALLABEL Function to reconstruct a photo volume using a probabilistic atlas as reference
%
% PARAMETERS
%
% inputPhotoDir: a directory with .tif / .mat pairs with the photos and segmentations
%
% inputREFERENCE: the probabilistic atlas. It should normally be ../prob_atlases/onlyCerebrum.smoothed.nii.gz
%
% outputVol: the output volume, in mgz format (please don't use .nii.gz because it has trouble with shearing in header
%
% outputVolMask: the corresponding output mask, which may be useful for subsequent processing
%
% outputWarpedRef: the probabilistic atlas after registration to the photos, which is useful to compute the samseg affine
%                     registration. Again, please use mgz extension.
%
% outputMat: a mat file that will store the history of the parameter optimization and of the cost function. Useful for
%                     subsequent visuzalization / making movies.
%
% PHOTO_RES: resoluton of the photos, in mm. In this project, it's 0.1; don't use any other value
%
% SLICE_THICKNESS: approximate thickness of the slices, in mm. In this project, it's 4; don't use any other value
%
% TARGET_RES: resolution of the photos for processing. It is a vector, such that each element is a resolution in the
%                     multiresolution pyramid. I use [4 2 1 0.5], but you can use something coarser to play with the code,
%                     e.g., [4 2 1]
%

% % TARGET_RES=[4 2 1];
% % Schedule
% % Mode 1: rigid for images (3*Nim), similarity for atlas (7)
% % Mode 2: rigid for images (3*Nim), affine for atlas (12)
% % Mode 3: penalized affine for images (6*Nim), affine for atlas (12)

% FS_MATLAB_PATH='/usr/local/freesurfer/dev/matlab';

%%%%%

REL_DICE_INTER_WEIGHT = 10; % 50;  % mask of reference to mask of photo
REL_DICE_INTRA_WEIGHT = 2/50; % mask of photo: slice N to N+1
REL_NCC_INTRA_WEIGHT = 1/50;  % ncc of photo: slice N to N+1
REL_DETERMINANT_COST = 0.1/50;  % determinant of affine transform of photos


%%%%%%%%%%%%%

% Number of pre/post slices to add at the photo stack
Nphotos_pre = 2;
Nphotos_post = 2;

%%%%%%%%%%%%%

tic

%%%%%%%%%%%%%

FREESURFER_HOME = getenv('FREESURFER_HOME');

if isempty(FREESURFER_HOME)
    error('please initialise code with UWphoto_startup.m')
end

fs_matlab_path = fullfile(FREESURFER_HOME,'matlab');

pathCell = regexp(path, pathsep, 'split');
if ispc  % Windows is not case-sensitive
  onPath = any(strcmpi(fs_matlab_path, pathCell));
else
  onPath = any(strcmp(fs_matlab_path, pathCell));
end

if ~onPath
   
    addpath(genpath(fs_matlab_path));
    
end


%%%%%%%%%%%%%%

%% Load slice of interest mask and labels

mask_nii = niftiread(niftiMask);
mask_nii = flip(flip(mask_nii',1),2);

label_nii = niftiread(niftiLabel);
label_nii = flip(flip(label_nii',1),2);

% this bit and the tests later rely on a format for the label names that
% could be changed later. Format currently
% <brainID>_Labels_ImageFile<n in dlist_masks>slice<slice in Masks_original>

[~,slicename,~] = fileparts(niftiLabel);

split_array = strsplit(slicename,{'ImageFile','slice','.nii'});

target_image = str2double(split_array{2});
target_slice = str2double(split_array{3});

%% read in photos and masks

disp('Extracting slices from photographs')
dlist_masks=dir([inputPhotoDir '/*.mat']);
Nphotos=length(dlist_masks);
Slices_original=[];
Masks_original=[];
grouping=[]; % I don't use it right now, but maybe in the future...
for n=1:Nphotos_pre
    Slices_original{end+1}=zeros(3,1); %#ok<AGROW>
    Masks_original{end+1}=1; %#ok<AGROW>
end
for n=1:Nphotos
    %% read in group image
    group_imags=imread([inputPhotoDir '/' dlist_masks(n).name(1:end-4) '.tif']);
    
    % load in matlab masks for slices in image
    load([inputPhotoDir '/' dlist_masks(n).name(1:end)],'LABELS');
    group_masks=LABELS;
    clear LABELS
    
    % make note of which group each slice was in
    grouping=[grouping n*ones(1,max(group_masks(:)))]; %#ok<AGROW>
    
    %% split up group into slices
    for l=1:max(group_masks(:))
        [mask,cropping_mask]=cropLabelVol(group_masks==l,5/PHOTO_RES);
        mask=imfill(mask,'holes');
        cropping=cropping_mask;
        cropping(3)=1; cropping(6)=3;
        image=applyCropping(group_imags,cropping);
        image(repmat(mask,[1 1 3])==0)=0;
        Slices_original{end+1}=image; %#ok<AGROW>
        Masks_original{end+1}=mask; %#ok<AGROW>
        
        if n==target_image && length(Masks_original)==target_slice
            
            crpd_msk_nii = applyCropping(mask_nii,cropping_mask);
            crpd_msk_nii = imfill(crpd_msk_nii,'holes');
            
            if isequal(crpd_msk_nii~=0,mask~=0)
                Slices_original{end} = uint16(Slices_original{end});
                Slices_original{end}(:,:,1)=applyCropping(label_nii,cropping_mask);
            else
                error('UWphoto:manualLabelVol:slicemismatch',...
                    'The slice mask provided does not match the specified slice');
            end
        end
    end
end
for n=1:Nphotos_post
    Slices_original{end+1}=zeros(3,1); %#ok<AGROW>
    Masks_original{end+1}=1; %#ok<AGROW>
end

%% re-order/resample slices

Nscales = length(TARGET_RES);
Nslices=length(Slices_original);
if exist([inputPhotoDir filesep '..' filesep 'slice_order.mat'],'file')
    load([inputPhotoDir filesep '..' filesep 'slice_order.mat'], 'slice_order');
    slice_order = [1:Nphotos_pre slice_order+Nphotos_pre slice_order(end)+Nphotos_pre+1:slice_order(end)+Nphotos_pre+Nphotos_post];
else
    slice_order = 1:Nslices;
end

I=[];
M=[];
disp(['Resampling to highest target resolutioscheduleITsn: ' num2str(TARGET_RES(Nscales)) ' mm']);
for n=1:Nslices
    
    n_ordered = slice_order(n);
    if n_ordered==target_slice
        target_slice = n;
    end
    
    I{n}=imresize(Slices_original{n_ordered},PHOTO_RES/TARGET_RES(Nscales),'nearest'); %#ok<AGROW>
    M{n}=imresize(double(Masks_original{n_ordered}),PHOTO_RES/TARGET_RES(Nscales))>0.5; %#ok<AGROW>
    % In the registration code this is where the mask is applied. Applying
    % this here in the label import code can cause labels outside of the
    % mask area to be incorrectly removed so we've commented it out and
    % left it for reference.
    %     I{n}(M{n}==0)=0; %#ok<AGROW> 
    if length(size(I{n})) < 3
        I{n} = zeros(3,1); %#ok<AGROW>
    end
end

%% find COGs of the masks, center, and pad

disp('Coarse alignment and padding');
Imri=[];
Mmri=[];
cogs=zeros(Nslices,2);
for n=1:Nslices
    [r,c]=find(M{n});
    if isempty(r)
        cogs(n,1)=1;
        cogs(n,2)=1;
    else
        cogs(n,1)=round(mean(r));
        cogs(n,2)=round(mean(c));
    end
end

semiLen = round(1.4 * max(cogs));
siz=1+2*semiLen;
Imri{Nscales}=[];
Imri{Nscales}.volres=[TARGET_RES(Nscales) TARGET_RES(Nscales) SLICE_THICKNESS];
Imri{Nscales}.vox2ras0=[-TARGET_RES(Nscales) 0 0 0; 0 0 -SLICE_THICKNESS 0; 0 -TARGET_RES(Nscales) 0 0; 0 0 0 1];
Imri{Nscales}.vol=zeros([siz Nslices 3]);
Mmri{Nscales}=Imri{Nscales};
Mmri{Nscales}.vol=zeros([siz Nslices]);
for n=1:Nslices
    idx1=semiLen-cogs(n,:);
    idx2=idx1+size(M{n})-1;
    Imri{Nscales}.vol(idx1(1):idx2(1),idx1(2):idx2(2),n,:)=reshape(I{n},[size(M{n}) 1 3]);
    Mmri{Nscales}.vol(idx1(1):idx2(1),idx1(2):idx2(2),n)=M{n};
end


%% left over resolution pyramid


disp('Building resolution pyramid');
for s=1:Nscales-1
    for n=1:Nslices
        % these downsamples could possibly be changed to nearest neighbour
        % as we are working with labels. I'm leaving them as the lower
        % resolutions don't appear to be used much.
        mri=Imri{Nscales}; mri.vol=mri.vol(:,:,n,1); mri=downsampleMRI2d(mri,TARGET_RES(s)/TARGET_RES(Nscales));
        if n==1
            Imri{s}.vol=zeros([size(mri.vol) Nslices 3]);
            Imri{s}.vox2ras0=mri.vox2ras0;
            Imri{s}.volres=mri.volres;
            Mmri{s}=Imri{s};
            Mmri{s}.vol=zeros([size(mri.vol) Nslices]);
        end
        Imri{s}.vol(:,:,n,1)=mri.vol;
        mri=Imri{Nscales}; mri.vol=mri.vol(:,:,n,2); mri=downsampleMRI2d(mri,TARGET_RES(s)/TARGET_RES(Nscales));
        Imri{s}.vol(:,:,n,2)=mri.vol;
        mri=Imri{Nscales}; mri.vol=mri.vol(:,:,n,3); mri=downsampleMRI2d(mri,TARGET_RES(s)/TARGET_RES(Nscales));
        Imri{s}.vol(:,:,n,3)=mri.vol;
        mri=Mmri{Nscales}; mri.vol=mri.vol(:,:,n); mri=downsampleMRI2d(mri,TARGET_RES(s)/TARGET_RES(Nscales));
        Mmri{s}.vol(:,:,n)=mri.vol>0.5;
    end
end


%% read in reference mri

REFmri=MRIread(inputREFERENCE);
REFmri.vol=REFmri.vol/max(REFmri.vol(:));
% Clean up fields other than volres, vox2ras0 and vol to avoid trouble...
aux=REFmri;
REFmri=[];
REFmri.vol=aux.vol;
REFmri.volres=aux.volres;
REFmri.vox2ras0=aux.vox2ras0;


%% Matching centres of gravity

% Let's start by matching the COGs
disp('Initializing with centers of gravity');
[IIref,JJref,KKref]=ndgrid(1:size(REFmri.vol,1),1:size(REFmri.vol,2),1:size(REFmri.vol,3));
rasRef=vox2ras([IIref(:) JJref(:) KKref(:)],REFmri.vox2ras0);
cogREF=sum((rasRef.*repmat(REFmri.vol(:)',[3 1])),2)/sum(REFmri.vol(:));

IIph=[]; JJph=[]; KKph=[];
for s=1:Nscales
    [a,b,c]=ndgrid(1:size(Imri{s}.vol,1),1:size(Imri{s}.vol,2),1:size(Imri{s}.vol,3));
    IIph{s}=a; JJph{s}=b; KKph{s}=c; %#ok<AGROW>
end
rasPH=vox2ras([IIph{Nscales}(:) JJph{Nscales}(:) KKph{Nscales}(:)],Imri{Nscales}.vox2ras0);
cogPH=sum((rasPH.*repmat(Mmri{Nscales}.vol(:)',[3 1])),2)/sum(Mmri{Nscales}.vol(:));

REFmri.vox2ras0(1:3,4)=REFmri.vox2ras0(1:3,4)+cogPH-cogREF;


%% Load optimised parameters from reconstruction run

load(paramMat,'paramsOptim')
mode = find(~cellfun('isempty',paramsOptim),1,'last');
x = paramsOptim{mode};

disp('Optimization done!');
if strcmpi(recontype,'soft')
    [~,~, warpedPhotos, warpedMasks, ~] = ...
        costFun_labels(x,cogREF,REFmri,Imri{Nscales},Mmri{Nscales},IIph{Nscales},JJph{Nscales},KKph{Nscales},...
        REL_NCC_INTRA_WEIGHT,REL_DICE_INTRA_WEIGHT,REL_DICE_INTER_WEIGHT,...
        REL_DETERMINANT_COST, mode);
else
    [~,~, warpedPhotos, warpedMasks, ~] = ...
        costFunHardRef_labels(x,cogREF,REFmri,Imri{Nscales},Mmri{Nscales},IIph{Nscales},JJph{Nscales},KKph{Nscales},...
        REL_NCC_INTRA_WEIGHT,REL_DICE_INTRA_WEIGHT,REL_DICE_INTER_WEIGHT,...
        REL_DETERMINANT_COST, mode);
end

warpedLabels = warpedMasks;
warpedLabels(:,:,target_slice)=round(squeeze(warpedPhotos(:,:,target_slice,1)));

disp('Writing results to disk...');
mri=Imri{Nscales};
mri.vol=warpedLabels;
MRIwrite(mri,outputlabel);

disp('All done!');
toc




