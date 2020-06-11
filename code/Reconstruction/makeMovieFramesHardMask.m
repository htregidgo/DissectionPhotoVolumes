% Script that writes to an output directory a bunch of files (volumes and surfaces) 
% that you can sequentially open with Freeview, to visualize the evolution
% of the fit bye ReconPhotoVolume_joint_hard_multires.m (and even make a 
% movie with it!).
% In full honesty, I haven't played to much with it since I switched to 
% multi-resolution, so it may need a bit of work / debugging here and there
clear
clc

% a directory with .tif / .mat pairs with the photos and segmentations
% inputPhotoDir='/autofs/cluster/vive/UW_photo_recon/Photo_data/18-1132/18-1132 MATLAB/';
inputPhotoDir='/home/henry/Documents/Brain/UWphoto/Photo_data_updated/17-0333/17-0333 MATLAB/';
% a reference binary mask volume, in correct anatomical orientation. You can use ../FLAIR_Scan_Data/*.rotated.mask.mgz.
% inputREFERENCE='/autofs/cluster/vive/UW_photo_recon/FLAIR_Scan_Data/NP18_1132.rotated.mask.mgz';
inputREFERENCE='/home/henry/Documents/Brain/UWphoto/FLAIR_Scan_Data/MaskComparison/NP17_0333/NP17_0333.rotated.ventCrctd.binary.mgz';
% output directory with movie files
outputMovieDir='/home/henry/Documents/Brain/UWphoto/figures/Movies/17-0333/';
% The mat file written by ReconPhotoVolume_joint_hard_multires.m
matFile='/home/henry/Documents/Brain/UWphoto/Results_hard/17-0333/17-0333.hard.mat';
% final resolution used in ReconPhotoVolume_joint_hard_multires.m
TARGET_RES=0.5;
% the path to the matlab directory of your freesurfer distrbution, i.e., $FREESURFER_HOME/matlab
FS_MATLAB_PATH='/usr/local/freesurfer/dev/matlab';
% Constants that you shouldn't need to touch
PHOTO_RES=0.1;
SLICE_THICKNESS=4;


%%%%%%%%%%%%%

% Number of pre/post slices to add at the photo stack
Nphotos_pre = 2;
Nphotos_post = 2;

%%%%%%%%%%%%%
% 
% addpath([pwd() '/functions']);
% addpath(FS_MATLAB_PATH);
% 
% %%%%%%%%%%%%%%

if exist(outputMovieDir,'dir')==0
    mkdir(outputMovieDir);
end

%%%%%%%%%%%%%%%%

disp('Extracting slices from photographs')
d=dir([inputPhotoDir '/*.mat']);
Nphotos=length(d);
Iorig=[];
Morig=[];
grouping=[]; % I don't use it right now, but maybe in the future...
for n=1:Nphotos_pre
    Iorig{end+1}=zeros(3,1);
    Morig{end+1}=1;
end
for n=1:Nphotos
    X=imread([inputPhotoDir '/' d(n).name(1:end-4) '.tif']);
    load([inputPhotoDir '/' d(n).name(1:end)],'LABELS'); Y=LABELS; clear LABELS
    grouping=[grouping n*ones(1,max(Y(:)))];
    for l=1:max(Y(:))
        [mask,cropping]=cropLabelVol(Y==l,5/PHOTO_RES);
        mask=imfill(mask,'holes');
        cropping(3)=1; cropping(6)=3;
        image=applyCropping(X,cropping);
        image(repmat(mask,[1 1 3])==0)=0;
        Iorig{end+1}=image;
        Morig{end+1}=mask;
    end
end
for n=1:Nphotos_post
    Iorig{end+1}=zeros(3,1);
    Morig{end+1}=1;
end

%%%%%%%%%%%%%%%

Nscales = length(TARGET_RES);
Nslices=length(Iorig);
if exist([inputPhotoDir filesep '..' filesep 'slice_order.mat'])
    load([inputPhotoDir filesep '..' filesep 'slice_order.mat'], 'slice_order');
    slice_order = [1:Nphotos_pre slice_order+Nphotos_pre slice_order(end)+Nphotos_pre+1:slice_order(end)+Nphotos_pre+Nphotos_post];
else
    slice_order = 1:Nslices;
end

I=[];
M=[];
disp(['Resampling to target resolution: ' num2str(TARGET_RES) ' mm']);
Nslices=length(Iorig);
for n=1:Nslices
    n_ordered = slice_order(n);
    I{n}=imresize(Iorig{n_ordered},PHOTO_RES/TARGET_RES);
%     M{n}=imdilate(imresize(double(M{n}),PHOTO_RES/TARGET_RES)>0.5,strel('disk',2/TARGET_RES));
    M{n}=imresize(double(Morig{n_ordered}),PHOTO_RES/TARGET_RES)>0.5;
    I{n}(M{n}==0)=0;
    if length(size(I{n})) < 3
        I{n} = zeros(3,1);
    end
end

%%%%%%%%%%%%%%%%%

disp('Coarse alignment and padding');
% find COGs of the masks
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

unalign=zeros(Nslices,2);
for n=1:Nslices
    [r,c]=find(true(size(M{n})));
    if isempty(r)
        unalign(n,1)=1;
        unalign(n,2)=1;
    else
        unalign(n,1)=round(mean(r));
        unalign(n,2)=round(mean(c));
    end
end

semiLen = round(1.4 * max(cogs));
siz=1+2*semiLen;
Imri=[];
Imri.volres=[TARGET_RES TARGET_RES SLICE_THICKNESS];
Imri.vox2ras0=[-TARGET_RES 0 0 0; 0 0 -SLICE_THICKNESS 0; 0 -TARGET_RES 0 0; 0 0 0 1];
Imri.vol=zeros([siz Nslices 3]);
Mmri=Imri;
Mmri.vol=zeros([siz Nslices]);
Umri=Imri;
for n=1:Nslices
    idx1=semiLen-cogs(n,:);
    idx2=idx1+size(M{n})-1;
    
    if (n>3 && unalign(n,1)~=1)
        idx4=idx3+size(M{n})-1;
    else
        idx3=semiLen-unalign(n,:);
        idx4=idx3+size(M{n})-1;
    end
    
    Imri.vol(idx1(1):idx2(1),idx1(2):idx2(2),n,:)=reshape(I{n},[size(M{n}) 1 3]);
    Mmri.vol(idx1(1):idx2(1),idx1(2):idx2(2),n)=M{n};
    Umri.vol(idx3(1):idx4(1),idx3(2):idx4(2),n,:)=reshape(I{n},[size(M{n}) 1 3]);
end


%%%%%%%%%%%%%%%%%%%%%
REFmri=MRIread(inputREFERENCE);
REFmri.vol=REFmri.vol/max(REFmri.vol(:));
% Clean up fields other than volres, vox2ras0 and vol to avoid trouble...
aux=REFmri;
REFmri=[];
REFmri.vol=aux.vol;
REFmri.volres=aux.volres;
REFmri.vox2ras0=aux.vox2ras0;

%%%%%%%%%%%%%%%%%%%
% THE ACTUAL WORK

% Let's start by matching the COGs
[IIref,JJref,KKref]=ndgrid(1:size(REFmri.vol,1),1:size(REFmri.vol,2),1:size(REFmri.vol,3));
rasRef=vox2ras([IIref(:) JJref(:) KKref(:)],REFmri.vox2ras0);
cogREF=sum((rasRef.*repmat(REFmri.vol(:)',[3 1])),2)/sum(REFmri.vol(:));

[IIph,JJph,KKph]=ndgrid(1:size(Imri.vol,1),1:size(Imri.vol,2),1:size(Imri.vol,3));
rasPH=vox2ras([IIph(:) JJph(:) KKph(:)],Imri.vox2ras0);
cogPH=sum((rasPH.*repmat(Mmri.vol(:)',[3 1])),2)/sum(Mmri.vol(:));

REFmri.vox2ras0(1:3,4)=REFmri.vox2ras0(1:3,4)+cogPH-cogREF;

% Here is where we diverge from the Recon script.
fn=0; % frame number\
factor=round(Mmri.volres(3)/2); 

% Initial
fn=fn+1;
MRIwrite(REFmri,[outputMovieDir '/REF_frame_' num2str(fn,'%.4d') '.mgz']);
MRIwrite(REFmri,'/tmp/kk.mgz'); system('mri_convert /tmp/kk.mgz /tmp/kk2.mgz --conform -odt float >/dev/null');
system('mri_mc /tmp/kk2.mgz 1 /tmp/kk >/dev/null');
system(['mris_smooth -nw /tmp/kk ' outputMovieDir '/REF_frame_' num2str(fn,'%.4d') '.surf >/dev/null']);

MRIwrite(Imri,[outputMovieDir '/SL_frame_' num2str(fn,'%.4d') '.mgz']);
MRIwrite(Umri,[outputMovieDir '/SL_frame_' num2str(0,'%.4d') '.mgz']);

MRIwrite(Mmri,'/tmp/kk.mgz'); system(['mri_convert /tmp/kk.mgz /tmp/kk2.mgz --voxsize 1.5 1.5 ' num2str(factor) ' -odt float -rt nearest >/dev/null']);
mri=MRIread('/tmp/kk2.mgz'); for z=1:size(mri.vol,3), mri.vol(:,:,z)=ceil((z+factor/2)/factor)*mri.vol(:,:,z); end;
MRIwrite(mri,[outputMovieDir '/SL_frame_' num2str(fn,'%.4d') '.labels.mgz']);



% Go around modes
Nims=size(Mmri.vol,3);
load(matFile,'historyX');
for mode=1:2
    
    for j=1:size(historyX{mode},1)
        
        [mode j size(historyX{mode},1)]
        
        [~,~, warpedPhotos, warpedMasks, REFvox2ras0New] = ...
            costFunHardRef(historyX{mode}(j,:),cogREF,REFmri,Imri,Mmri,IIph,JJph,KKph,...
            0,0,0,0, mode);
        
        fn=fn+1;
        
        mri=REFmri; mri.vox2ras0=REFvox2ras0New; mri.volres=sqrt(sum(REFvox2ras0New(1:3,1:3).^2));
        MRIwrite(mri,[outputMovieDir '/REF_frame_' num2str(fn,'%.4d') '.mgz']);
        MRIwrite(mri,'/tmp/kk.mgz'); system('mri_convert /tmp/kk.mgz /tmp/kk2.mgz --conform -odt float >/dev/null');
        system('mri_mc /tmp/kk2.mgz 1 /tmp/kk >/dev/null');
        system(['mris_smooth -nw /tmp/kk ' outputMovieDir '/REF_frame_' num2str(fn,'%.4d') '.surf  >/dev/null']);
        
        mri=Imri; mri.vol=warpedPhotos;
        MRIwrite(mri,[outputMovieDir '/SL_frame_' num2str(fn,'%.4d') '.mgz'])
        
        mri=Mmri; mri.vol=warpedMasks>.5;
        MRIwrite(mri,'/tmp/kk.mgz'); system(['mri_convert /tmp/kk.mgz /tmp/kk2.mgz --voxsize 1.5 1.5 ' num2str(factor) ' -odt float -rt nearest >/dev/null']);
        mri=MRIread('/tmp/kk2.mgz'); for z=1:size(mri.vol,3), mri.vol(:,:,z)=ceil((z+factor/2)/factor)*mri.vol(:,:,z); end;
        MRIwrite(mri,[outputMovieDir '/SL_frame_' num2str(fn,'%.4d') '.labels.mgz']);
    end
end


% -v SL_frame_0001.mgz -v SL_frame_0001.labels.mgz:colormap=lut  -f REF_frame_0001.surf:color=255,0,0

disp('All done!');





