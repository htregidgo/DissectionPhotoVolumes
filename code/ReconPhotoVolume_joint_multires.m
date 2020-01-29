% Function to reconstruct a photo volume using a probabilistic atlas as reference
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
% scheduleITs: this is a matrix specifying the schedule. Must have as many rows as resolutions (i.e., the number of elements
%                     in TARGET_RES), and 3 columns. Each element is a number of iterations. The column indices 1, 2 and 3
%                     correspond to the 3 different modes of complexity of the registration:
%                         Mode 1: rigid for images (3*Nim parameters), similarity for atlas (7)
%                         Mode 2: rigid for images (3*Nim), affine for atlas (12)
%                         Mode 3: penalized affine for images (6*Nim), affine for atlas (12)
%                     I normally use [50 40 30; 25 20 15; 12 10 8; 6 5 5]. If you are playing with fewer resolutions
%                     TARGET_RES = [4 2 1], then you need 3 rows instead, for instance [45 30 20; 15 10 10; 5 5 5]
%
% FS_MATLAB_PATH: the path to the matlab directory of your freesurfer distrbution, i.e., $FREESURFER_HOME/matlab (e.g., something like /usr/local/freesurfer/matlab)
%

function ReconPhotoVolume_joint_multires(inputPhotoDir,inputREFERENCE,outputVol,...
    outputVolMask,outputWarpedRef,outputMat,PHOTO_RES,SLICE_THICKNESS,...
    TARGET_RES,scheduleITs,FS_MATLAB_PATH)

% clear
% clc
%
% inputPhotoDir='/autofs/cluster/vive/UW_photo_recon/Photo_data/18-1132/18-1132 MATLAB/';
% inputREFERENCE='/autofs/cluster/vive/UW_photo_recon/prob_atlases/onlyCerebrum.smoothed.nii.gz';
% outputVol='/autofs/cluster/vive/UW_photo_recon/recons/outputsProbAtlasBin/18-1132.recon.mgz';
% outputVolMask='/autofs/cluster/vive/UW_photo_recon/recons/outputsProbAtlasBin/18-1132.mask.mgz';
% outputWarpedRef='/autofs/cluster/vive/UW_photo_recon/recons/outputsProbAtlasBin/18-1132.warped_ref.mgz';
% outputMat='/autofs/cluster/vive/UW_photo_recon/recons/outputsProbAtlasBin/18-1132.mat';
% PHOTO_RES=0.1;
% SLICE_THICKNESS=4;
% TARGET_RES=[4 2 1 0.5];
% % TARGET_RES=[4 2 1];
% % Schedule
% % Mode 1: rigid for images (3*Nim), similarity for atlas (7)
% % Mode 2: rigid for images (3*Nim), affine for atlas (12)
% % Mode 3: penalized affine for images (6*Nim), affine for atlas (12)
% scheduleITs = [50 40 30; 25 20 15; 12 10 8; 6 5 5];
% % scheduleITs = [45 30 20; 15 10 10; 5 5 5];
% FS_MATLAB_PATH='/usr/local/freesurfer/dev/matlab';
%

%%%%%

REL_DICE_INTER_WEIGHT = 10; % 50;  % mask of reference to mask of photo
REL_DICE_INTRA_WEIGHT = 2/50; % mask of photo: slice N to N+1
REL_NCC_INTRA_WEIGHT = 1/50;  % ncc of photo: slice N to N+1
REL_DETERMINANT_COST = 0.1/50;  % determinant of affine transform of photos

%%%%%%%%%%%%%

% DON'T TOUCH THIS OR YOU'LL MESS UP the OPTIMIZATION. OR IF YOU DO, MAKE
% SURE YOU ALSO CHANGE IT IN THE COST FUNCTION costFun.m
FACTOR_AFFINE_MAT=20;
FACTOR_SCALING = 20;

%%%%%%%%%%%%%

% Number of pre/post slices to add at the photo stack
Nphotos_pre = 2;
Nphotos_post = 2;

%%%%%%%%%%%%%

tic

%%%%%%%%%%%%%

addpath(FS_MATLAB_PATH);
addpath([fileparts(mfilename('fullpath')) filesep 'functions/lbfgsb3.0_mex1.2/']);
addpath([fileparts(mfilename('fullpath')) filesep 'functions']);

%%%%%%%%%%%%%%
if strcmp(outputWarpedRef(end-3:end),'.mgz')==0
    error('Output warped reference volume must be a mgz file to support shear in the vox2ras matrix');
end

%%%%%%%%%%%%%%


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
disp(['Resampling to highest target resolutioscheduleITsn: ' num2str(TARGET_RES(Nscales)) ' mm']);
for n=1:Nslices
    
    n_ordered = slice_order(n);
    I{n}=imresize(Iorig{n_ordered},PHOTO_RES/TARGET_RES(Nscales));
    % M{n}=imdilate(imresize(double(Morig{n}),PHOTO_RES/TARGET_RES(Nscales))>0.5,strel('disk',ceil(2/TARGET_RES(Nscales))));
    M{n}=imresize(double(Morig{n_ordered}),PHOTO_RES/TARGET_RES(Nscales))>0.5;
    I{n}(M{n}==0)=0;
        if length(size(I{n})) < 3
        I{n} = zeros(3,1);
    end
end

%%%%%%%%%%%%%%%%%

% find COGs of the masks, center, and pad
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


%%%%%%%%%%%%%%%%%%%%%


disp('Building resolution pyramid');
for s=1:Nscales-1
    for n=1:Nslices
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
disp('Initializing with centers of gravity');
[IIref,JJref,KKref]=ndgrid(1:size(REFmri.vol,1),1:size(REFmri.vol,2),1:size(REFmri.vol,3));
rasRef=vox2ras([IIref(:) JJref(:) KKref(:)],REFmri.vox2ras0);
cogREF=sum((rasRef.*repmat(REFmri.vol(:)',[3 1])),2)/sum(REFmri.vol(:));

IIph=[]; JJph=[]; KKph=[];
for s=1:Nscales
    [a,b,c]=ndgrid(1:size(Imri{s}.vol,1),1:size(Imri{s}.vol,2),1:size(Imri{s}.vol,3));
    IIph{s}=a; JJph{s}=b; KKph{s}=c;
end
rasPH=vox2ras([IIph{Nscales}(:) JJph{Nscales}(:) KKph{Nscales}(:)],Imri{Nscales}.vox2ras0);
cogPH=sum((rasPH.*repmat(Mmri{Nscales}.vol(:)',[3 1])),2)/sum(Mmri{Nscales}.vol(:));

REFmri.vox2ras0(1:3,4)=REFmri.vox2ras0(1:3,4)+cogPH-cogREF;

% Let's go
Nims=size(Imri{1}.vol,3);
Ngroups=Nphotos;
paramsOptim=cell(1,32);
historyX=cell(1,3);
historyCost=cell(1,3);

for mode=1:3
    
    disp('*****************');
    disp(['*    MODE  ' num2str(mode) '    *']);
    disp('*****************');
    
    historyX{mode}=[];
    historyCost{mode}=[];
    
    for s=1:Nscales
        opts=[];
        opts.maxIts = scheduleITs(s,mode); % scheduled
        opts.maxTotalIts = 30 * opts.maxIts;
        opts.printEvery     = 1;
        % opts.verbose = -1; % default is -1, i.e., no outpuversion with hard maskt from
        opts.m  = 5; % should be between 3 and 20; default is 5
        
        if s==1 % first scale
            
            if mode==1 % in first mode, simply take all zero (easy!)
                opts.x0=zeros([3*Nims+7,1]);
                
            elseif mode==2 % in mode 2, we need to move from ref-similarity to ref-affine
                
                params  = zeros([3*Nims+12,1]);
                
                % photos are the same, but we need to scale translations!
                params(1:3*Nims)=x(1:3*Nims);
                params(2:3*Nims)=x(2:3*Nims)/TARGET_RES(s)*TARGET_RES(Nscales);
                params(3:3*Nims)=x(3:3*Nims)/TARGET_RES(s)*TARGET_RES(Nscales);
                
                % translation of atlas is the same, no need to scale, as it
                % is in RAS already
                params(end-2:end) = x(end-2:end);
                
                % Turn rotation into M
                scaling=exp(x(end-6)/FACTOR_SCALING);
                rotx=x(end-5)/180*pi;
                roty=x(end-4)/180*pi;
                rotz=x(end-3)/180*pi;
                T1=[1 0 0; 0 cos(rotx) -sin(rotx); 0 sin(rotx) cos(rotx)];
                T2=[cos(roty) 0 sin(roty); 0 1 0; -sin(roty) 0 cos(roty)];
                T3=[cos(rotz) -sin(rotz) 0; sin(rotz) cos(rotz) 0; 0 0 1];
                M = scaling * T3 * T2 * T1;
                params(end-11:end-3) = reshape(M,[9 1])*FACTOR_AFFINE_MAT;
                
                opts.x0=params;
                
            else  % in mode 3, we need to move from rigid + affine to  affine + affine
                
                params  = zeros([6*Nims+12,1]);
                
                % Reference is the same
                params(end-11:end) = x(end-11:end);
                
                % now for the photos
                % first rotation -> affine
                theta=x(1:3:end-12)/180*pi;
                for i = 1:Nims
                    M=[cos(theta(i)) -sin(theta(i)); sin(theta(i)) cos(theta(i))];
                    params(4*i-3:4*i)=reshape(M,[4 1])*FACTOR_AFFINE_MAT;
                end
                % translation is the same, but we need to scale
                tr=x(2:3:end-12)/TARGET_RES(s)*TARGET_RES(Nscales);
                tc=x(3:3:end-12)/TARGET_RES(s)*TARGET_RES(Nscales);
                params(4*Nims+1:2:6*Nims)=tr;
                params(4*Nims+2:2:6*Nims)=tc;
                
                opts.x0=params;
            end
            
        else % rest of scales: we simply scale the translation parameters as needed!
            
            if mode==1
                idx=sort([2:3:length(x)-7 3:3:length(x)-7]);
            elseif mode==2
                idx=sort([2:3:length(x)-12 3:3:length(x)-12]);
            else
                idx=sort([4*Nims+1:2:length(x)-12 4*Nims+2:2:length(x)-12]);
            end
            opts.x0=x;
            opts.x0(idx)=opts.x0(idx)/TARGET_RES(s)*TARGET_RES(Nscales);
        end
        
        
        disp(['Running ' num2str(scheduleITs(s,mode)) ' iterations of BFGS at scale ' num2str(s) ' of ' num2str(Nscales)]);
        n=length(opts.x0);
        u = Inf*ones(n,1);
        l = -u;
        [x,~,info] = lbfgsb(  @(p)costFun(p,cogREF,REFmri,Imri{s},Mmri{s},...
            IIph{s},JJph{s},KKph{s},REL_NCC_INTRA_WEIGHT,REL_DICE_INTRA_WEIGHT,...
            REL_DICE_INTER_WEIGHT,REL_DETERMINANT_COST,mode),l, u, opts );
        
        
        % Scale translations if needed, as they are in pixels
        if s<Nscales
            if mode==1
                idx=sort([2:3:length(x)-7 3:3:length(x)-7]);
            elseif mode==2
                idx=sort([2:3:length(x)-12 3:3:length(x)-12]);
            else
                idx=sort([4*Nims+1:2:length(x)-12 4*Nims+2:2:length(x)-12]);
            end
            x(idx)=x(idx)*TARGET_RES(s)/TARGET_RES(Nscales);
%             info.xs(:,idx)=info.xs(:,idx)*TARGET_RES(s)/TARGET_RES(Nscales);
        else
            paramsOptim{mode}=x;
        end
        
%         historyX{mode}=[historyX{mode}; info.xs];
        historyCost{mode}=[historyCost{mode}; info.err(:,1) ];
    end
    disp(' ');
end

disp('Optimization done!');
[~,~, warpedPhotos, warpedMasks, REFvox2ras0New] = ...
    costFun(x,cogREF,REFmri,Imri{Nscales},Mmri{Nscales},IIph{Nscales},JJph{Nscales},KKph{Nscales},...
    REL_NCC_INTRA_WEIGHT,REL_DICE_INTRA_WEIGHT,REL_DICE_INTER_WEIGHT,...
    REL_DETERMINANT_COST, mode);

disp('Writing results to disk...');
mri=Imri{Nscales};
mri.vol=warpedPhotos;
MRIwrite(mri,outputVol);
mri.vol=warpedMasks;
MRIwrite(mri,outputVolMask);
REFmri.vox2ras0=REFvox2ras0New;
REFmri.volres=sqrt(sum(Imri{Nscales}.vox2ras0(1:3,1:3).^2));
MRIwrite(REFmri,outputWarpedRef);
save(outputMat,'paramsOptim','historyX','historyCost');

disp('All done!');
toc




