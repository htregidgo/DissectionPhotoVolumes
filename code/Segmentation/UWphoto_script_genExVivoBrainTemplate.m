%% set filepaths

fpth_template = ['/home/henry/Documents/Brain/Freesurfer/freesurfer/',...
    'average/samseg/20Subjects_smoothing2_down2_smoothingForAffine2/',...
    'template.nii'];
fpth_tmpltSeg = ['/home/henry/Documents/Brain/UWphoto/test_template/',...
    'template_crispSegmentation.nii'];


%% read in volumes

ninfo_in = niftiinfo(fpth_template);

template_full = niftiread(fpth_template);
template_seg  = niftiread(fpth_tmpltSeg);

%% get mask

template_mask = (template_seg~=0 & template_seg<98);
template_bckg = template_seg==0;

%% close holes

se = strel('sphere',5);

template_clsd = imclose(template_mask,se);

%% get cropped template

template_cropped = template_full;

fill_val = median(template_full(template_bckg(:)));

template_cropped(~template_clsd)=fill_val;

%% output cropped template

fpth_out = fullfile(pwd,'/data/template_exvivo');
ninfo_out = rmfield(ninfo_in,{'Filemoddate','Filename','Filesize'});

niftiwrite(template_cropped,fpth_out,ninfo_out)

%% output mask

fpth_out = fullfile(pwd,'/data/template_exvivo_mask');

ninfo_out = rmfield(ninfo_in,{'Filemoddate','Filename','Filesize'});
ninfo_out.Datatype='uint8';

outvol = zeros(size(template_clsd),'uint8');
outvol(template_clsd)=1;

niftiwrite(outvol,fpth_out,ninfo_out)



%% same for no brainstem

template_mask = (template_seg~=0 & template_seg<98);
template_mask = template_mask & template_seg~=7 & template_seg~=8; % left cerebelum
template_mask = template_mask & template_seg~=46 & template_seg~=47; % right cerebelum
template_mask = template_mask & template_seg~=15 & template_seg~=16; % BS & 4th ventricle
template_mask = template_mask & template_seg~=28 & template_seg~=60; % ventral DC
template_bckg = template_seg==0;

template_CSF = template_seg==24;


%% close holes

se = strel('sphere',5);

template_clsd = imclose(template_mask&~template_CSF,se);


se = strel('sphere',3);

template_dltd = imdilate(template_clsd,se);

%% get cropped template

template_cropped = template_full;

fill_val = median(template_full(template_bckg(:)));

template_cropped(~template_dltd)=fill_val;

%% output cropped template

fpth_out = fullfile(pwd,'/data/template_cerebrumOnly');
ninfo_out = rmfield(ninfo_in,{'Filemoddate','Filename','Filesize'});

niftiwrite(template_cropped,fpth_out,ninfo_out)

%% output mask

fpth_out = fullfile(pwd,'/data/template_cerebrumOnly_mask');

ninfo_out = rmfield(ninfo_in,{'Filemoddate','Filename','Filesize'});
ninfo_out.Datatype='uint8';

outvol = zeros(size(template_dltd),'uint8');
outvol(template_dltd)=1;

niftiwrite(outvol,fpth_out,ninfo_out)

