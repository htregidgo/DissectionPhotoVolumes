% Script to take the Manual label volumes and samseg labels and produce
% label specific dice scores. Changes to the file structure may require
% modification of the 'setup for file structure' section. Script requires
% the code base be initialised with UWphoto_startup.m.
%
% dice scores are produced in a vector with scores placed according to the
% label value. ie dice_scores(<label value WM>).

if ~exist('forceFlag','var')
    forceFlag=false;
end

%% setup for file structure

PHOTO_RECON_HOME=getenv('PHOTO_RECON_HOME');

top_soft_dir = fullfile(PHOTO_RECON_HOME,'Results');
top_hard_dir = fullfile(PHOTO_RECON_HOME,'Results_hard/');

general_soft_path = fullfile(top_soft_dir,'*','*','*soft_manualLabel.mgz');
general_hard_path = fullfile(top_hard_dir,'*','*hard_manualLabel.mgz');

top_scores_dir = fullfile(PHOTO_RECON_HOME,'figures','diceScores');

if ~exist(top_scores_dir,'dir')
    mkdir(top_scores_dir)
end



%% find segmentation volumes
% - manual label volumes
% - matching samseg volumes

dlist_manLabelsSft = dir(general_soft_path);
dlist_manLabelsHrd = dir(general_hard_path);

dlist_manLabelsAll = vertcat(dlist_manLabelsSft,dlist_manLabelsHrd);


flag_isSftRecon = false(size(dlist_manLabelsAll));
flag_isSftRecon(1:length(dlist_manLabelsSft)) = true;

%% go through cases

for il=1:length(dlist_manLabelsAll)
    
    %% get matching files
    
    [~,infile,~] = fileparts(dlist_manLabelsAll(il).name);
    brainID = infile(1:7);
    
    pth_label_manual = fullfile(dlist_manLabelsAll(il).folder,...
        dlist_manLabelsAll(il).name);
    
    dlist_masks  = dir(fullfile(dlist_manLabelsAll(il).folder,'*mask.mgz'));
    
    pth_masks = fullfile(dlist_masks.folder,dlist_masks.name);
    
    if flag_isSftRecon(il)
        dlist_samseg = dir(fullfile(top_soft_dir,'SAMSEG',[brainID,'*'],...
            '*crispSegmentation*'));
    else
        dlist_samseg = dir(fullfile(top_hard_dir,'SAMSEG',[brainID,'*'],...
            '*crispSegmentation*'));
    end
    
    pth_label_samseg = fullfile(dlist_samseg.folder,dlist_samseg.name);
    
    
    vol_label_manual = MRIread(pth_label_manual);
    vol_label_samseg = MRIread(pth_label_samseg);
    vol_masks        = MRIread(pth_masks);
    
    
    pth_merged_manual = fullfile(dlist_manLabelsAll(il).folder,...
        [infile,'_merged.mgz']);
    
    pth_merged_samseg = fullfile(dlist_manLabelsAll(il).folder,...
        [infile(1:12),'_samseg_merged.mgz']);
    
    if exist(pth_merged_manual,'file') && exist(pth_merged_samseg,'file')...
            && ~forceFlag
        vol_merged_manual = MRIread(pth_merged_manual);
        vol_merged_samseg = MRIread(pth_merged_samseg);
    else
        
        %% merge labels - pre-defined
        % - ventral DC & WM
        
        values_merged_manual = vol_label_manual.vol;
        values_merged_samseg = vol_label_samseg.vol;
        
        % left-ventral-dc   28
        % left-cerebral-WM  2
        
        values_merged_manual(values_merged_manual==28)=2;
        values_merged_samseg(values_merged_samseg==28)=2;
        
        % right-ventral-dc  60
        % right-cerebral-WM 41
        
        values_merged_manual(values_merged_manual==60)=41;
        values_merged_samseg(values_merged_samseg==60)=41;
        
        %% merge labels - find nearest mach
        % - WM hypo-intensities
        % - choroid plexus and ventricles (or background)
        
        % WM-hyp-intensities 77
        % left-wm            2
        % right-wm           41
        
        se = strel('sphere',1);
        
        WM_hypo_full = values_merged_samseg==77;
        
        if any(WM_hypo_full(:))
            
            ConnComp = bwconncomp(WM_hypo_full);
            
            for jl = 1: length(ConnComp.PixelIdxList)
                
                comp_test_vol = false(size(WM_hypo_full));
                comp_test_vol(ConnComp.PixelIdxList{jl})=true;
                
                comp_test_dltd = imdilate(comp_test_vol,se);
                
                sum_left = sum(values_merged_samseg(comp_test_dltd(:))==2);
                sum_right = sum(values_merged_samseg(comp_test_dltd(:))==41);
                
                if sum_left>sum_right
                    values_merged_samseg(comp_test_vol)=2;
                else
                    values_merged_samseg(comp_test_vol)=41;
                end
            end
            
        end
        
        % left-choroid-plexus    31
        % left-lat-ventricle     4
        
        LCP_full = values_merged_samseg==31;
        
        if any(LCP_full(:))
            
            ConnComp = bwconncomp(LCP_full);
            
            for jl = 1: length(ConnComp.PixelIdxList)
                
                comp_test_vol = false(size(LCP_full));
                comp_test_vol(ConnComp.PixelIdxList{jl})=true;
                
                comp_test_dltd = imdilate(comp_test_vol,se);
                
                if any(values_merged_samseg(comp_test_dltd(:))==4)
                    values_merged_samseg(comp_test_vol)=4;
                else
                    values_merged_samseg(comp_test_vol)=5;
                end
            end
            
        end
        
        % right-choroid-plexus   63
        % right-lat-ventricle    43
        
        RCP_full = values_merged_samseg==63;
        
        if any(RCP_full(:))
            
            ConnComp = bwconncomp(RCP_full);
            
            for jl = 1: length(ConnComp.PixelIdxList)
                
                comp_test_vol = false(size(RCP_full));
                comp_test_vol(ConnComp.PixelIdxList{jl})=true;
                
                comp_test_dltd = imdilate(comp_test_vol,se);
                
                if any(values_merged_samseg(comp_test_dltd(:))==43)
                    values_merged_samseg(comp_test_vol)=43;
                else
                    values_merged_samseg(comp_test_vol)=44;
                end
            end
        end
        
        %% write out merged volumes
        
        vol_merged_manual = vol_label_manual;
        vol_merged_samseg = vol_label_samseg;
        
        vol_merged_manual.vol = values_merged_manual;
        vol_merged_samseg.vol = values_merged_samseg;
        
        MRIwrite(vol_merged_manual,pth_merged_manual);
        MRIwrite(vol_merged_samseg,pth_merged_samseg);
        
    end
    
    %% generate Dice by species
    
    slice_location = find(vol_merged_manual.vol==2,1);
    [~,~,slice_number] = ind2sub(size(vol_merged_manual.vol),slice_location);
    
    slice_dice_manual = squeeze(vol_merged_manual.vol(:,:,slice_number));
    slice_dice_samseg = squeeze(vol_merged_samseg.vol(:,:,slice_number));
    
    
    dice_scores = dice(slice_dice_manual,slice_dice_samseg);
    
    
    %% output stats to files
    % - in figures
    
    pth_diceScores_mat = fullfile(top_scores_dir,[infile(1:12),'_dice.mat']);
    save(pth_diceScores_mat,'dice_scores','slice_dice_*')
    
    
end