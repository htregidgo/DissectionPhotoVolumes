UWphoto_function_checkBashEnvironment();

nthread = 32;

%% soft segmentation

photovols_dir =  uigetdir(pwd,'Top level results directory');
dlist_Photovols_soft = dir(fullfile(photovols_dir,'*/soft','*_soft.mgz'));

PHOTO_RECON_HOME = getenv('PHOTO_RECON_HOME');

codedir = fullfile(PHOTO_RECON_HOME,'code');

templatedir = fullfile(codedir,'data');

segmnteddir = fullfile(photovols_dir,'SAMSEG/');

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
    
    
    if ~exist(ptntdir,'dir')
        mkdir(ptntdir)
    end
        
    inputVol = fullfile(dlist_Photovols_soft(il).folder,dlist_Photovols_soft(il).name);
    
    dlist_mask = dir(fullfile(dlist_Photovols_soft(il).folder,'*soft_mask.mgz'));
    
    inputVolMask = fullfile(dlist_mask.folder,dlist_mask.name);
        
    dlist_ref = dir(fullfile(dlist_Photovols_soft(il).folder,'*soft_regatlas.mgz'));
    
    inputWarpedRef = fullfile(dlist_ref.folder,dlist_ref.name);
    
    %% 
    try
        [segstatus,segresult]=UWphoto_function_SegmentPhotosWithSAMSEG(...
            caseID,inputVol,inputVolMask,inputWarpedRef,ptntdir,nthread,...
            'Hemisphere','Whole');
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
end