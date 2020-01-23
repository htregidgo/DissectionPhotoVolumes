%% setup freesurfer
FSdir=getenv('FREESURFER_HOME');

if isempty(FSdir)
    error('FreeSurfer needs to be installed and sourced');
end

addpath([FSdir '/matlab/']);

%% set directories

data_dir = '/home/henry/Documents/Brain/UWphoto/Results_hard';
output_dir = fullfile(data_dir,'Niftis4seg');

if~exist(output_dir,'dir')
    mkdir(output_dir);
end

dlist_masks  = dir(fullfile(data_dir,'*','Mask*png'));
dlist_photos = dir(fullfile(data_dir,'*','Segmentation*png'));

resolution = 0.1;

if isnan(resolution)
    error('Resolution must be numeric');
end

%% output niftis

for il=1:length(dlist_photos)
    
    [~,brainID,~] =fileparts(dlist_photos(il).folder);
    
    %% photo
    % Read file
    imageFile=fullfile(dlist_photos(il).folder, dlist_photos(il).name);
    
    %   try
    disp(imageFile)
    I = imread(imageFile);
    disp(['Working on file ' num2str(il) ' of ' num2str(length(dlist_photos))])
    
    % Rotate image
    I = imrotate(I,180);
    
    I = reshape(I,[size(I,1) size(I,2) 1 size(I,3)]);
    mri=[];
    mri.vox2ras0=eye(4);
    mri.vox2ras0(1,1)=resolution;
    mri.vox2ras0(2,2)=resolution;
    mri.volres=[resolution resolution 1];
    mri.vol=I;
    
    MRIwrite(mri,fullfile(output_dir,...
        [brainID,'_',dlist_photos(il).name(1:end-4),'.nii.gz']));
    
    
    %% mask
    % Read file
    imageFile=fullfile(dlist_masks(il).folder, dlist_masks(il).name);
    
    %   try
    disp(imageFile)
    I = imread(imageFile);
    disp(['Working on file ' num2str(il) ' of ' num2str(length(dlist_masks))])
    
    % Rotate image
    I = imrotate(I,180);
    
    I = reshape(I,[size(I,1) size(I,2) 1 size(I,3)]);
    mri=[];
    mri.vox2ras0=eye(4);
    mri.vox2ras0(1,1)=resolution;
    mri.vox2ras0(2,2)=resolution;
    mri.volres=[resolution resolution 1];
    mri.vol=I;
    
    MRIwrite(mri,fullfile(output_dir,...
        [brainID,'_',dlist_masks(il).name(1:end-4),'.nii.gz']));
    
end

disp('All done!');