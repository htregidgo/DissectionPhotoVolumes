% Script to pre-register the ex-vivo template file to the individual flair
% volumes. This will require niftyreg to have been added to the path using
% setenv. It will also require the use of mri_convert from freesurfer.

%% find files

flair_dir = '/home/henry/Documents/Brain/UWphoto/FLAIR_Scan_Data';

CB_BS_file = '/home/henry/Documents/Brain/UWphoto/code/data/BS_CB_cases.mat';
load(CB_BS_file,'BS_CB_cases')

dlist_mris = dir(fullfile(flair_dir,'*.rotated.mgz'));

MRI_CONVERT = '/home/henry/Documents/Brain/Freesurfer/freesurfer/bin/mri_convert';

templatedir = '/home/henry/Documents/Brain/UWphoto/code/data';
exvivotempt = fullfile(templatedir,'template_exvivo.nii');
cerbrmtempt = fullfile(templatedir,'template_cerebrumOnly.nii');
exvivomask = fullfile(templatedir,'template_exvivo_mask.nii');
cerbrmmask = fullfile(templatedir,'template_cerebrumOnly_mask.nii');
floatingtempt = fullfile(templatedir,'template_floating.nii');
registerdir = fullfile(templatedir,'MRItemplates');

if ~exist(registerdir,'dir')
    mkdir(registerdir)
end


if ~exist(exvivotempt,'file')
    error('construct ex-vivo template first')    
end
if ~exist(cerbrmtempt,'file')
    error('construct ex-vivo template without brain stem etc. first')    
end
    
%% convert and register files

for il=1:length(dlist_mris)
   
    %% settup filenames for case
    [~,basename,~] = fileparts(dlist_mris(il).name);
    [~,caseID,~]   = fileparts(basename);
    
    MRI_mgz_version = fullfile(flair_dir,dlist_mris(il).name);
    MRI_nii_version = fullfile(flair_dir,[basename,'.nii.gz']);
    Msk_mgz_version = fullfile(flair_dir,[basename,'.binary.mgz']);
    MSK_nii_version = fullfile(flair_dir,[basename,'_mask.nii.gz']);
    
    transform = fullfile(registerdir,[caseID,'_afftrans.txt']);
    resampled = fullfile(registerdir,[caseID,'_template_resample.nii.gz']);
    movedhder = fullfile(registerdir,[caseID,'_template.mgz']);
    
    %% get nifti version of flair for reg_aladin
    if ~exist(MRI_nii_version,'file')
        
        cmd = [MRI_CONVERT,' ',MRI_mgz_version,' ',MRI_nii_version];
        
        [cnvrtStatus,cnvrtResult]=system(cmd);
        
    end
    
    %% get nifti version of flair for reg_aladin
    if ~exist(MSK_nii_version,'file')
        
        cmd = [MRI_CONVERT,' ',Msk_mgz_version,' ',MSK_nii_version];
        
        [cnvrtStatus,cnvrtResult]=system(cmd);
        
    end
    
    %% do affine registration with reg_aladin
    % requires niftyreg directory added to path using setenv
    if ~exist(resampled,'file') || endsWith(caseID,BS_CB_cases)
        
        if endsWith(caseID,BS_CB_cases)
%             mri_a = MRIread(exvivotempt);
%             mri_b = [];
%             mri_b.vol = mri_a.vol;
%             mri_b.vox2ras0 = mri_a.vox2ras0;
%             
%             mri_b.vox2ras0(1:3,1:3)=0.9*mri_b.vox2ras0(1:3,1:3);
%             
%             mri_b.volres=sqrt(sum(mri_b.vox2ras0(1:3,1:3).^2));
%             
%             MRIwrite(mri_b,floatingtempt);
            
            cmd = ['reg_aladin -ref ',MRI_nii_version,' -flo ',exvivotempt,...
                ' -fmask ',cerbrmmask, ' -rmask ',MSK_nii_version,...
                ' -aff ',transform,' -res ',resampled];
            
        else
%             mri_a = MRIread(cerbrmtempt);
%             mri_b = [];
%             mri_b.vol = mri_a.vol;
%             mri_b.vox2ras0 = mri_a.vox2ras0;
%             
%             mri_b.vox2ras0(1:3,1:3)=0.9*mri_b.vox2ras0(1:3,1:3);
%             
%             mri_b.volres=sqrt(sum(mri_b.vox2ras0(1:3,1:3).^2));
%             
%             MRIwrite(mri_b,floatingtempt)
            
            cmd = ['reg_aladin -ref ',MRI_nii_version,' -flo ',cerbrmtempt,...
                ' -fmask ',cerbrmmask, ' -rmask ',MSK_nii_version,...
                ' -aff ',transform,' -res ',resampled];
        end
        
        fprintf(['Registering case ',caseID,'\n'])
        [regStatus,regResult]=system(cmd);
        
        if regStatus
            warning(['problem with registration of ',caseID])
        end
    else
        fprintf(['Skipping case ',caseID,'\n'])
    end
   
    
    %% apply affine transform to template header
    % read in template and output moved version as mgz
    if ~exist(movedhder,'file') || endsWith(caseID,BS_CB_cases)
        fprintf('outputing registered template\n')
        T = dlmread(transform);
        
        if endsWith(caseID,BS_CB_cases)
            A = MRIread(exvivotempt);
        else
            A = MRIread(cerbrmtempt);
        end

%         A = MRIread(floatingtempt);

        B = [];
        
        B.vol = A.vol;
        B.vox2ras0 = T\A.vox2ras0;
        
        B.volres=sqrt(sum(B.vox2ras0(1:3,1:3).^2));
        
        MRIwrite(B,movedhder);
    end
end