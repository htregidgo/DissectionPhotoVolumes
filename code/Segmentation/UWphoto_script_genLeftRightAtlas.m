%% set filepaths

fpth_template = ['/home/henry/Documents/Brain/UWphoto/prob_atlases/',...
    'onlyCerebrum.nii.gz'];
fpth_tmpltSeg = ['/home/henry/Documents/Brain/UWphoto/code/temp/',...
    'onlyCerebrum.nii_crispSegmentation.nii'];


%% read in volumes

template_in = MRIread(fpth_template);

template_full = template_in.vol;

seg_in = MRIread(fpth_tmpltSeg);

template_seg  = seg_in.vol;

%% get mask

template_bckg = template_full==0;
template_mask = ~template_bckg;

left_mask_original = template_mask & template_seg<40 & template_seg~=14 & ...
    template_seg~=24;

right_mask_original = template_mask & template_seg>39 & template_seg<72 ;

template_boarder = false(size(template_mask));


%% iteratively dilate contour
se = strel('sphere',5);
flag_continue = true;

left_mask_old = imopen(left_mask_original,se);
right_mask_old = imopen(right_mask_original,se);

while flag_continue
   
    right_mask_d = imdilate(right_mask_old,se);
    left_mask_d  = imdilate(left_mask_old,se);
    
    template_boarder = left_mask_d & right_mask_d;
    
    right_mask = right_mask_old | (right_mask_d & ~template_boarder);
    left_mask = left_mask_old | (left_mask_d & ~template_boarder);
    
    flag_continue = any(template_mask(:) & ~(left_mask(:) | right_mask(:)...
        | template_boarder(:)));
    
    right_mask_old = right_mask;
    left_mask_old = left_mask;
    
end

se = strel('sphere',1);
flag_continue = true;

while flag_continue
   
    right_mask_d = imdilate(right_mask_old,se);
    left_mask_d  = imdilate(left_mask_old,se);
    
    template_boarder = left_mask_d & right_mask_d;
    
    right_mask = right_mask_old | (right_mask_d & ~template_boarder);
    left_mask = left_mask_old | (left_mask_d & ~template_boarder);
    
    flag_continue = any(template_mask(:) & ~(left_mask(:) | right_mask(:)...
        | template_boarder(:)));
    
    right_mask_old = right_mask;
    left_mask_old = left_mask;
    
end

%% get output

right_atlas = template_full;
right_atlas(~right_mask)=0;

left_atlas = template_full;
left_atlas(~left_mask)=0;




%% output cropped template

[in_path,in_name,~]=fileparts(fpth_template);
[~,in_name,~]=fileparts(in_name);

left_fpth_out = fullfile(in_path,['left_',in_name]);
right_fpth_out = fullfile(in_path,['right_',in_name]);


left_out = template_in;
right_out = template_in;

left_out.vol=left_atlas;
right_out.vol=right_atlas;

MRIwrite(left_out,[left_fpth_out,'.nii.gz']);
MRIwrite(right_out,[right_fpth_out,'.nii.gz']);


left_out.vol=imgaussfilt3(left_atlas,1.5);
right_out.vol=imgaussfilt3(right_atlas,1.5);

MRIwrite(left_out,[left_fpth_out,'.smoothed.nii.gz']);
MRIwrite(right_out,[right_fpth_out,'.smoothed.nii.gz']);




