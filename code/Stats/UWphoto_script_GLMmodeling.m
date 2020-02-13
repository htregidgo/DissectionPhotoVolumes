% Script to generate GLMs with the matched case info and volume statistics
% from both MRI and photo segmentations

%% load data

PHOTO_RECON_HOME=getenv('PHOTO_RECON_HOME');

figuresDir = fullfile(PHOTO_RECON_HOME,'figures');

load(fullfile(figuresDir,'AdjustedCaseStats.mat'),'segVolumeInfo','matchedInfo')


%% sort out cases to remove

removalLevel = 2;

removalFlag = (0<[segVolumeInfo.removaltype]) &...
    ([segVolumeInfo.removaltype] <=removalLevel);

segVolumeInfo_kept=segVolumeInfo(~removalFlag);


%% get target demographics

flag_mrs = strcmp({segVolumeInfo_kept.segtype},'MRI');
flag_hrd = strcmp({segVolumeInfo_kept.segtype},'Hard');
flag_sft = strcmp({segVolumeInfo_kept.segtype},'Soft');


flag_female_disease = strcmp(matchedInfo.sex,'Female') & ...
    strcmp(matchedInfo.CognitiveStatus,'Dementia');

fd_demographics = matchedInfo(flag_female_disease,{'caseID','age'});

flag_male_disease = strcmp(matchedInfo.sex,'Male') & ...
    strcmp(matchedInfo.CognitiveStatus,'Dementia');

md_demographics = matchedInfo(flag_male_disease,{'caseID','age'});

flag_male_control = strcmp(matchedInfo.sex,'Male') & ...
    strcmp(matchedInfo.CognitiveStatus,'No dementia');

mc_demographics = matchedInfo(flag_male_control,{'caseID','age'});

%% make info structures

mri_strct_fd = struct();
mri_strct_fd(size(fd_demographics,1)) = struct();

mri_strct_md = struct();
mri_strct_md(size(md_demographics,1)) = struct();

mri_strct_mc = struct();
mri_strct_mc(size(mc_demographics,1)) = struct();


hrd_strct_fd = struct();
hrd_strct_fd(size(fd_demographics,1)) = struct();

hrd_strct_md = struct();
hrd_strct_md(size(md_demographics,1)) = struct();

hrd_strct_mc = struct();
hrd_strct_mc(size(mc_demographics,1)) = struct();


crct_strct_fd = struct();
crct_strct_fd(size(fd_demographics,1)) = struct();

crct_strct_md = struct();
crct_strct_md(size(md_demographics,1)) = struct();

crct_strct_mc = struct();
crct_strct_mc(size(mc_demographics,1)) = struct();


sft_strct_fd = struct();
sft_strct_fd(size(fd_demographics,1)) = struct();

sft_strct_md = struct();
sft_strct_md(size(md_demographics,1)) = struct();

sft_strct_mc = struct();
sft_strct_mc(size(mc_demographics,1)) = struct();


strfields = fieldnames(segVolumeInfo_kept);

%% get info fd

keep_flag = true(size(fd_demographics,1),1);

for il = 1:size(fd_demographics,1)
    
    selection_flag = strcmp({segVolumeInfo_kept.caseID},fd_demographics.caseID(il));
    
    if any(selection_flag & flag_mrs)
        mri_strct_fd(il).age=fd_demographics.age(il);
        
        hrd_strct_fd(il).age=fd_demographics.age(il);
        
        crct_strct_fd(il).age=fd_demographics.age(il);
        
        sft_strct_fd(il).age=fd_demographics.age(il);
        
        for jl=1:length(strfields)
            mri_strct_fd(il).(strfields{jl})=segVolumeInfo_kept(...
                selection_flag & flag_mrs).(strfields{jl});
            
            hrd_strct_fd(il).(strfields{jl})=segVolumeInfo_kept(...
                selection_flag & flag_hrd).(strfields{jl});
            
            if jl<5
                crct_strct_fd(il).(strfields{jl})=segVolumeInfo_kept(...
                    selection_flag & flag_hrd).(strfields{jl});
            else
                crct_strct_fd(il).(strfields{jl})=segVolumeInfo_kept(...
                    selection_flag & flag_hrd).(strfields{jl})...
                    .*segVolumeInfo_kept(selection_flag & flag_hrd);
            end
            
            
            sft_strct_fd(il).(strfields{jl})=segVolumeInfo_kept(strcmp({segVolumeInfo_kept.caseID},...
                fd_demographics.caseID(il)) & flag_sft).(strfields{jl});
        end
    else
        keep_flag(il)=false;
    end
    
    
    
end

mri_strct_fd=mri_strct_fd(keep_flag);

hrd_strct_fd=hrd_strct_fd(keep_flag);

sft_strct_fd=sft_strct_fd(keep_flag);

%% get info md

keep_flag = true(size(md_demographics,1),1);

for il = 1:size(md_demographics,1)
          
    if any(strcmp({segVolumeInfo_kept.caseID},md_demographics.caseID(il)) &...
            flag_mrs)
        mri_strct_md(il).age=md_demographics.age(il);
        
        hrd_strct_md(il).age=md_demographics.age(il);
        
        sft_strct_md(il).age=md_demographics.age(il);
        
        for jl=1:length(strfields)
            mri_strct_md(il).(strfields{jl})=segVolumeInfo_kept(strcmp({segVolumeInfo_kept.caseID},...
                md_demographics.caseID(il)) & flag_mrs).(strfields{jl});
            
            hrd_strct_md(il).(strfields{jl})=segVolumeInfo_kept(strcmp({segVolumeInfo_kept.caseID},...
                md_demographics.caseID(il)) & flag_hrd).(strfields{jl});
            
            sft_strct_md(il).(strfields{jl})=segVolumeInfo_kept(strcmp({segVolumeInfo_kept.caseID},...
                md_demographics.caseID(il)) & flag_sft).(strfields{jl});
        end
    else
        keep_flag(il)=false;
    end
    
    
    
end

mri_strct_md=mri_strct_md(keep_flag);

hrd_strct_md=hrd_strct_md(keep_flag);

sft_strct_md=sft_strct_md(keep_flag);

%% get info mc

keep_flag = true(size(mc_demographics,1),1);

for il = 1:size(mc_demographics,1)
          
    if any(strcmp({segVolumeInfo_kept.caseID},mc_demographics.caseID(il)) &...
            flag_mrs)
        mri_strct_mc(il).age=mc_demographics.age(il);
        
        hrd_strct_mc(il).age=mc_demographics.age(il);
        
        sft_strct_mc(il).age=mc_demographics.age(il);
        
        for jl=1:length(strfields)
            mri_strct_mc(il).(strfields{jl})=segVolumeInfo_kept(strcmp({segVolumeInfo_kept.caseID},...
                mc_demographics.caseID(il)) & flag_mrs).(strfields{jl});
            
            hrd_strct_mc(il).(strfields{jl})=segVolumeInfo_kept(strcmp({segVolumeInfo_kept.caseID},...
                mc_demographics.caseID(il)) & flag_hrd).(strfields{jl});
            
            sft_strct_mc(il).(strfields{jl})=segVolumeInfo_kept(strcmp({segVolumeInfo_kept.caseID},...
                mc_demographics.caseID(il)) & flag_sft).(strfields{jl});
        end
    else
        keep_flag(il)=false;
    end
    
    
    
end

mri_strct_mc=mri_strct_mc(keep_flag);

hrd_strct_mc=hrd_strct_mc(keep_flag);

sft_strct_mc=sft_strct_mc(keep_flag);


%% build GLM design matrix
% columns correspond to age, gender and disease respectively. For gender 1
% corresponds to female 0 to male. For disease 1 is dementia 0 is control.
% There is also a column of ones to account for the constant.

% numbers in each cohort 
% N_fd - female diseased
N_fd = length(mri_strct_fd);
% N_md - male diseased  
N_md = length(mri_strct_md);
% N_mc - male control
N_mc = length(mri_strct_mc);

N_total = N_fd+N_md+N_mc;

% predictors will be age, gender and disease, with a constant term
Design_mat = zeros(N_total,4);
Design_mat(:,4)=1;

% females with disease
for il= 1:N_fd
   
    Design_mat(il,1) = mri_strct_fd(il).age;
    Design_mat(il,2) = 1;
    Design_mat(il,3) = 1;
    
end

% males with disease
for il= 1:N_md
    
    Ientry = il+N_fd;
   
    Design_mat(Ientry,1) = mri_strct_md(il).age;
    Design_mat(Ientry,2) = 0;
    Design_mat(Ientry,3) = 1;
    
end

% males without disease
for il= 1:N_mc
    
    Ientry = il+N_fd+N_md;
   
    Design_mat(Ientry,1) = mri_strct_mc(il).age;
    Design_mat(Ientry,2) = 0;
    Design_mat(Ientry,3) = 0;
    
end

% remove column for gender due to lack of female controls
Design_mat_maleOnly=Design_mat(N_fd+1:end,[1,3,4]);

%% run model 
% not going to use volume corrections at the moment because I'm not sure
% they are working well.

volume_fields = fieldnames(mri_strct_fd);

volume_flag = ~ismember(volume_fields,{'age','caseID','segtype',...
    'volumeAdjustmentFactor'}');

volume_fields = volume_fields(volume_flag);


MRI_GLM_resultsStruct = struct();
Hrd_GLM_resultsStruct = struct();
Sft_GLM_resultsStruct = struct();

for il=1:length(volume_fields)
    
    Observations = [mri_strct_fd(:).(volume_fields{il}),...
        mri_strct_md(:).(volume_fields{il}),...
        mri_strct_mc(:).(volume_fields{il})]';
        
%     [b,dev,stats] = glmfit(Design_mat,Observations);
    Design_mat_maleOnly=[Design_mat(N_fd+1:end,[1,3]) ones(size(Design_mat_maleOnly,1),1)];
    [b,dev,stats] = glmfit(Design_mat_maleOnly,Observations(N_fd+1:end),'normal',...
        'constant','off');
    
    MRI_GLM_resultsStruct.(volume_fields{il}).b = b;
    MRI_GLM_resultsStruct.(volume_fields{il}).dev = dev;
    MRI_GLM_resultsStruct.(volume_fields{il}).stats = stats;
    
    
    Observations = [hrd_strct_fd(:).(volume_fields{il}),...
        hrd_strct_md(:).(volume_fields{il}),...
        hrd_strct_mc(:).(volume_fields{il})]';
        
%     [b,dev,stats] = glmfit(Design_mat,Observations);
    Design_mat_maleOnly=Design_mat(N_fd+1:end,[1,3]);
    Design_mat_maleOnly=[Design_mat_maleOnly ones(size(Design_mat_maleOnly,1),1)];
    [b,dev,stats] = glmfit(Design_mat_maleOnly,Observations(N_fd+1:end),'normal',...
        'constant','off');
    
    Hrd_GLM_resultsStruct.(volume_fields{il}).b = b;
    Hrd_GLM_resultsStruct.(volume_fields{il}).dev = dev;
    Hrd_GLM_resultsStruct.(volume_fields{il}).stats = stats;
    
    temp = Observations(N_fd+1:end);
    
    stats.p = ranksum(temp(Design_mat_maleOnly(:,2)==0),temp(Design_mat_maleOnly(:,2)==1));
    
    Observations = [sft_strct_fd(:).(volume_fields{il}),...
        sft_strct_md(:).(volume_fields{il}),...
        sft_strct_mc(:).(volume_fields{il})]';
        
%     [b,dev,stats] = glmfit(Design_mat,Observations);
    Design_mat_maleOnly=Design_mat(N_fd+1:end,[1,3]);
    Design_mat_maleOnly=[Design_mat_maleOnly ones(size(Design_mat_maleOnly,1),1)];
    [b,dev,stats] = glmfit(Design_mat_maleOnly,Observations(N_fd+1:end),'normal',...
        'constant','off');
    
    temp = Observations(N_fd+1:end);
    
    stats.p = ranksum(temp(Design_mat_maleOnly(:,2)==0),temp(Design_mat_maleOnly(:,2)==1));
    
    Sft_GLM_resultsStruct.(volume_fields{il}).b = b;
    Sft_GLM_resultsStruct.(volume_fields{il}).dev = dev;
    Sft_GLM_resultsStruct.(volume_fields{il}).stats = stats;
    
    
    
    
    
end