%% set filepaths

fpth_template = ['/home/henry/Documents/Brain/Freesurfer/freesurfer/',...
    'average/samseg/20Subjects_smoothing2_down2_smoothingForAffine2/',...
    'template.nii'];
fpth_tmpltSeg = ['/home/henry/Documents/Brain/UWphoto/test_template/',...
    'template_crispSegmentation.nii'];


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

%% get parts to output

template_mask = (template_seg~=0 & template_seg<98);
template_mask = template_mask & template_seg~=7 & template_seg~=8; % left cerebelum
template_mask = template_mask & template_seg~=46 & template_seg~=47; % right cerebelum
template_mask = template_mask & template_seg~=15 & template_seg~=16; % BS & 4th ventricle
template_mask = template_mask & template_seg~=28 & template_seg~=60; % ventral DC
template_bckg = template_seg==0;

template_CSF = template_seg==24;


se = strel('sphere',5);
template_opened = imopen(template_mask&~template_CSF,se);
template_closed = imclose(template_mask&~template_CSF,se);

%% separate left and right

right_atlas = template_closed;
right_atlas(~right_mask)=0;

left_atlas = template_closed;
left_atlas(~left_mask)=0;




%% output cropped template

PHOTO_RECON_HOME = getenv('PHOTO_RECON_HOME');

left_fpth_out = fullfile(PHOTO_RECON_HOME,'code','data','template_left.nii');
right_fpth_out = fullfile(PHOTO_RECON_HOME,'code','data','template_right.nii');


left_out = template_in;
right_out = template_in;

left_out.vol=left_atlas;
right_out.vol=right_atlas;

MRIwrite(left_out,left_fpth_out);
MRIwrite(right_out,right_fpth_out);





