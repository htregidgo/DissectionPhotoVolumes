% Script to overlay the samseg segmentations of the ex-vivo MRI volumes on
% to the Hard registered photo volumes. 


%% setup directories 

MRIseg_directory = '/home/henry/Documents/Brain/UWphoto/FLAIR_Scan_Data/SAMSEG';

Photo_directory = '/home/henry/Documents/Brain/UWphoto/Results_hard';


%% get file list 

% get list of warped references used for headers
dlist_photos = dir(fullfile(Photo_directory,'*/*.warped_ref.mgz'));


%% work through files 

for il = 1:length(dlist_photos)
    
    %% get matching segmentation
    [~,brainID,~] = fileparts(dlist_photos(il).folder);
    
    brainID_conv = strrep(brainID,'-','_');
    
    dlist_mrSamSeg = dir(fullfile(MRIseg_directory,...
        ['*',brainID_conv,'/*_crispSegmentation.nii']));
    
    outfile = fullfile(dlist_photos(il).folder,...
        [brainID,'.hard.warped_seg.mgz']);
    
    %% read in volumes
    photo_ref = MRIread(fullfile(dlist_photos(il).folder,...
        dlist_photos(il).name));
    
    mri_seg = MRIread(fullfile(dlist_mrSamSeg.folder,dlist_mrSamSeg.name));
    
    % transfer header information from photo volumes to output segmentation
    photo_seg = photo_ref;
    
    %% output warped segmentation
    
    % transfer mri segmentation voxel data to output segmentation
    photo_seg.vol = mri_seg.vol;
    
    MRIwrite(photo_seg,outfile);
    
    
end

