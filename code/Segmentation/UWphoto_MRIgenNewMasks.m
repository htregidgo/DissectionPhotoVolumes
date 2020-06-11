% script to create an MRI mask from the samseg segmentations.

%% get files

MRI_dir = '/home/henry/Documents/Brain/UWphoto/FLAIR_Scan_Data';
Seg_dir = '/home/henry/Documents/Brain/UWphoto/FLAIR_Scan_Data/SAMSEG';

Pht_dir = '/home/henry/Documents/Brain/UWphoto/Results_hard';

out_dir_top = fullfile(MRI_dir,'MaskComparison');

if ~exist(out_dir_top,'dir')
    mkdir(out_dir_top)
end

dlist_mris = dir(fullfile(MRI_dir,'*rotated.mgz'));

%%
[stt,~]=system('mri_convert -h');

if stt
    fprintf('freesurfer bin not on path. Adding now!\n')
    
    FREESURFER_HOME = getenv('FREESURFER_HOME');
    
    dlist_bin=dir(fullfile(FREESURFER_HOME,'bin','mri_convert'));
    
    if ~exist(fullfile(dlist_bin.folder,'mri_convert'),'file') 
        error('can''t find mri_convert')
    end

    
    setenv('PATH',[getenv('PATH'),':',dlist_bin.folder])
    
    
end


%% create masks

for il = 1:length(dlist_mris)
    
    %% get matching files
    brainID = dlist_mris(il).name(1:9);
    
    dlist_binary = dir(fullfile(MRI_dir,[brainID,'*binary*']));
    
    dlist_seg = dir(fullfile(Seg_dir,brainID,'*crispSegmentation*'));
    
    out_dir = fullfile(out_dir_top,brainID);
    
    if ~exist(out_dir,'dir')
        mkdir(out_dir)
    end
    
    if ~exist(fullfile(out_dir,[brainID,'.rotated.newbinary.mgz']),'file')
        %% copy existing files
        
        copyfile(fullfile(dlist_mris(il).folder,dlist_mris(il).name),...
            fullfile(out_dir,dlist_mris(il).name))
        
        copyfile(fullfile(dlist_binary.folder,dlist_binary.name),...
            fullfile(out_dir,dlist_binary.name))
        
        copyfile(fullfile(dlist_seg.folder,dlist_seg.name),...
            fullfile(out_dir,dlist_seg.name))
        
        %% open volumes
        
        in_binary = MRIread(fullfile(dlist_binary.folder,dlist_binary.name));
        
        in_seg = MRIread(fullfile(dlist_seg.folder,dlist_seg.name));
        
        in_mri = MRIread(fullfile(dlist_mris(il).folder,dlist_mris(il).name));
        
        out_mask = in_binary;
        
        %% generate new mask
        
        new_mask = (in_seg.vol~=0 & in_seg.vol<98);
        new_mask = new_mask & in_seg.vol~=7 & in_seg.vol~=8; % left cerebelum
        new_mask = new_mask & in_seg.vol~=46 & in_seg.vol~=47; % right cerebelum
        new_mask = new_mask & in_seg.vol~=15 & in_seg.vol~=16; % BS & 4th ventricle
        new_mask = new_mask & in_seg.vol~=28 & in_seg.vol~=60; % ventral DC
        new_bckg = in_seg.vol==0;
        
        new_CSF = in_seg.vol==24;
        
        %% fill holes
        
        mask_filled = imfill(new_mask&~new_CSF,'holes');
        
        %% output new mask
        
        fill_val=median(in_binary.vol(in_binary.vol(:)~=0));
        
        out_mask.vol(:)=0;
        out_mask.vol(mask_filled(:))=fill_val;
        
        MRIwrite(out_mask,fullfile(out_dir,[brainID,'.rotated.newbinary.mgz']));
        
    else
        
        %% set up for intersection of ventral-DC with photo-mask
        
        altID = strrep(brainID(3:end),'_','-');
        
        dlist_phtMask = dir(fullfile(Pht_dir,altID,[altID,'.hard.mask.mgz']));
             
        if isempty(dlist_phtMask)  
        
            warning('missing hard photo volume mask, continuing.')
            continue
        end
        
        photomask = fullfile(dlist_phtMask.folder,dlist_phtMask.name);
        
        mrimask = fullfile(dlist_phtMask.folder,[altID,'.hard.mask.resampled.mgz']);
        
        dlist_warped = dir(fullfile(dlist_phtMask.folder,'*.warped_ref.mgz'));
        
        warpedmask = fullfile(dlist_warped.folder,dlist_warped.name);
        
        cmnd = ['mri_convert ',photomask,' ',mrimask,' -rl ',warpedmask,...
            ' -rt nearest'];
        
        system(cmnd)
        
        %%
        in_binary = MRIread(fullfile(dlist_binary.folder,dlist_binary.name));
        
        in_seg = MRIread(fullfile(dlist_seg.folder,dlist_seg.name));
        
        in_phtmsk = MRIread(mrimask);
        
        out_mask = in_binary;
        
        %% generate new mask
        
        new_mask = (in_seg.vol~=0 & in_seg.vol<98);
        new_mask = new_mask & in_seg.vol~=7 & in_seg.vol~=8; % left cerebelum
        new_mask = new_mask & in_seg.vol~=46 & in_seg.vol~=47; % right cerebelum
        new_mask = new_mask & in_seg.vol~=15 & in_seg.vol~=16; % BS & 4th ventricle
        
        vent_dc_excld = (in_seg.vol==28 | in_seg.vol==60) & ~in_phtmsk.vol;
        
        new_mask = new_mask & ~vent_dc_excld; % ventral DC
        new_bckg = in_seg.vol==0;
        
        new_CSF = in_seg.vol==24;
        
        %% fill holes
        
        mask_filled = imfill(new_mask&~new_CSF,'holes');
        
        se = strel('sphere',1);
        
        mask_closed = imclose(mask_filled,se);
        
        %% output new mask
        
        fill_val=median(in_binary.vol(in_binary.vol(:)~=0));
        
        out_mask.vol(:)=0;
        out_mask.vol(mask_closed(:))=fill_val;
        
        MRIwrite(out_mask,fullfile(out_dir,[brainID,'.rotated.ventCrctd.binary.mgz']));
        
        
    end
    
end






