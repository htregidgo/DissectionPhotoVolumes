UWphoto_function_checkBashEnvironment();

nthread = 32;
% 
% photovols_dir =  '/home/henry/Documents/Brain/UWphoto/Results';
% dlist_Photovols = dir(fullfile(photovols_dir,'*/soft','*_soft.mgz'));

photovols_dir =  '/home/henry/Documents/Brain/UWphoto/Results_hard';
dlist_Photovols = dir(fullfile(photovols_dir,'*','*.hard.recon.mgz'));

codedir = '/home/henry/Documents/Brain/UWphoto/code';

templatedir = '/home/henry/Documents/Brain/UWphoto/code/data';
registerdir = fullfile(templatedir,'MRItemplates');

segmnteddir = fullfile(photovols_dir,'SAMSEGMRBIAS/');
oldsegmnteddir = fullfile(photovols_dir,'SAMSEG');

if ~exist(segmnteddir,'dir')
    mkdir(segmnteddir)
end

%% run registration
problems_flag = true(size(dlist_Photovols));
segstatus_cell = cell(size(dlist_Photovols));
segresult_cell = cell(size(dlist_Photovols));

for il=1:length(dlist_Photovols)
   
    
    %% settup filenames for case
    [~,basename,~] = fileparts(dlist_Photovols(il).name);
    [~,caseID,~]   = fileparts(basename);
    
    ptntdir = fullfile(segmnteddir,caseID);
    prevdir = fullfile(oldsegmnteddir,caseID);
    
    
    if ~exist(ptntdir,'dir')
        mkdir(ptntdir)
    end
        
    inputVol = fullfile(dlist_Photovols(il).folder,dlist_Photovols(il).name);
    
    dlist_mask = dir(fullfile(dlist_Photovols(il).folder,'*.hard.mask.mgz'));
    
    inputVolMask = fullfile(dlist_mask.folder,dlist_mask.name);
    
    
    dlist_ref = dir(fullfile(dlist_Photovols(il).folder,'*.hard.warped_ref.mgz'));
    
    inputWarpedRef = fullfile(dlist_ref.folder,dlist_ref.name);
    
    dlist_output = dir(fullfile(ptntdir,'*_brightnessCorrected.nii.gz'));
    if isempty(dlist_output) || (exist('forceflag','var')&&forceflag==true)
        %%
        try
            [segstatus,segresult]=UWphoto_function_SegmentPhotosWithSAMSEG_3Dbias(...
                caseID,inputVol,inputVolMask,inputWarpedRef,ptntdir,nthread,...
                'PrevRunDir',prevdir);
            if segstatus
                warning(['problem with segmentation of ',caseID])
                
                segstatus_cell{il} = segstatus;
                segresult_cell{il} = segresult;
            else
                problems_flag(il)=false;
                fprintf(segresult(end-35:end))
            end
        catch
            warning(['problem with segmentation of ',caseID])
        end
    else 
        fprintf(['skiping ',caseID,'\n'])
        problems_flag(il)=false;
    end
end


%% soft segmentation

photovols_dir =  '/home/henry/Documents/Brain/UWphoto/Results';
dlist_Photovols_soft = dir(fullfile(photovols_dir,'*/soft','*_soft.mgz'));

codedir = '/home/henry/Documents/Brain/UWphoto/code';

templatedir = '/home/henry/Documents/Brain/UWphoto/code/data';
registerdir = fullfile(templatedir,'MRItemplates');

segmnteddir = fullfile(photovols_dir,'SAMSEGMRBIAS/');
oldsegmnteddir = fullfile(photovols_dir,'SAMSEG');

if ~exist(segmnteddir,'dir')
    mkdir(segmnteddir)
end

%% run registration
problems_flag_soft = true(size(dlist_Photovols_soft));
segstatus_cell_soft = cell(size(dlist_Photovols_soft));
segresult_cell_soft = cell(size(dlist_Photovols_soft));

for il=1:length(dlist_Photovols_soft)
   
    
    %% settup filenames for case
    [~,basename,~] = fileparts(dlist_Photovols_soft(il).name);
    [~,caseID,~]   = fileparts(basename);
    
    ptntdir = fullfile(segmnteddir,caseID);
    prevdir = fullfile(oldsegmnteddir,caseID);
    
    
    if ~exist(ptntdir,'dir')
        mkdir(ptntdir)
    end
        
    inputVol = fullfile(dlist_Photovols_soft(il).folder,dlist_Photovols_soft(il).name);
    
    dlist_mask = dir(fullfile(dlist_Photovols_soft(il).folder,'*soft_mask.mgz'));
    
    inputVolMask = fullfile(dlist_mask.folder,dlist_mask.name);
    
    
    dlist_ref = dir(fullfile(dlist_Photovols_soft(il).folder,'*soft_regatlas.mgz'));
    
    inputWarpedRef = fullfile(dlist_ref.folder,dlist_ref.name);
    
    dlist_output = dir(fullfile(ptntdir,'*_brightnessCorrected.nii.gz'));
    if isempty(dlist_output) || (exist('forceflag','var')&&forceflag==true)
        %%
        try
            [segstatus,segresult]=UWphoto_function_SegmentPhotosWithSAMSEG_3Dbias(...
                caseID,inputVol,inputVolMask,inputWarpedRef,ptntdir,nthread,...
                'PrevRunDir',prevdir);
            if segstatus
                warning(['problem with segmentation of ',caseID])
                
                segstatus_cell_soft{il} = segstatus;
                segresult_cell_soft{il} = segresult;
            else
                problems_flag_soft(il)=false;
                fprintf(segresult(end-35:end))
            end
        catch
            warning(['problem with segmentation of ',caseID])
        end
    else 
        fprintf(['skiping ',caseID,'\n'])
        problems_flag_soft(il)=false;
    end
end