% Script to generate GLMs with the matched case info and volume statistics
% from both MRI and photo segmentations

clear all

%% load data
removalLevel = 1;
GLM_RESULTS = ['/home/acasamitjana/Data/UWphoto/PowerAnalysis/GLM/GLMresults_removalLevel_' num2str(removalLevel) '_new.mat'];
load(GLM_RESULTS, 'MRI_GLM_resultsStruct', 'Hrd_GLM_resultsStruct', 'Sft_GLM_resultsStruct', 'Design_mat'); X=Design_mat; clear Design_mat;


ALPHA = 0.05;
POWER = 0.8;
N_COVARIATES = 4;
GAMMAX = inv(X'*X);

%% MRI power analysis
disp('MRI analysis')
PA_MRI.ALPHA = ALPHA;
PA_MRI.POWER = POWER;

N_MRI = MRI_GLM_resultsStruct.N_total;
MRI_results = MRI_GLM_resultsStruct.full;
volume_fields = fieldnames(MRI_results);

MRI_col = zeros(length(volume_fields),1);
for it_volume=1:length(volume_fields)
    stats = MRI_results.(volume_fields{it_volume}).stats;
    
    beta = stats.beta(3);
    gammaX=GAMMAX(3,3);
    sigma2 = sum(stats.resid.*stats.resid);
     
    [nout, tval, pval, power] = computeSampleSize(0, beta, sigma2, gammaX, N_COVARIATES, POWER, ALPHA);
    
    PA_MRI.full.(volume_fields{it_volume}).power = power;
    PA_MRI.full.(volume_fields{it_volume}).pval = pval;
    PA_MRI.full.(volume_fields{it_volume}).tval = tval;
    PA_MRI.full.(volume_fields{it_volume}).sample_size = nout;
    
%     if pval>ALPHA
%         disp(['MRI_ALPHA_' volume_fields{it_volume} ])
%     end
%     if power<POWER
%         disp(['MRI_POWER_' volume_fields{it_volume} ])
%     end
    
    MRI_col(it_volume) = nout; 
    
end

%% Hard power analysis
disp('Hard analysis')
PA_Hard.ALPHA = ALPHA;
PA_Hard.POWER = POWER;

N_Hard = Hrd_GLM_resultsStruct.N_total;
Hard_results = Hrd_GLM_resultsStruct.full;
volume_fields = fieldnames(Hard_results);

Hard_col = zeros(length(volume_fields),1);
for it_volume=1:length(volume_fields)
    stats = Hard_results.(volume_fields{it_volume}).stats;
    
    beta = stats.beta(3);
    gammaX=GAMMAX(3,3);
    sigma2 = sum(stats.resid.*stats.resid);
     
    [nout, tval, pval, power] = computeSampleSize(0, beta, sigma2, gammaX, N_COVARIATES, POWER, ALPHA);
    
    PA_Hard.full.(volume_fields{it_volume}).power = power;
    PA_Hard.full.(volume_fields{it_volume}).pval = pval;
    PA_Hard.full.(volume_fields{it_volume}).tval = tval;
    PA_Hard.full.(volume_fields{it_volume}).sample_size = nout;

%     if pval>ALPHA
%         disp(['Hard_ALPHA_' volume_fields{it_volume} ])
%     end
%     if power<POWER
%         disp(['Hard_POWER_' volume_fields{it_volume} ])
%     end
    
    
    Hard_col(it_volume) = nout; 
end



%% Soft power analysis
disp('Soft analysis')
PA_Soft.ALPHA = ALPHA;
PA_Soft.POWER = POWER;

N_Soft = Sft_GLM_resultsStruct.N_total;
Soft_results = Sft_GLM_resultsStruct.full;
volume_fields = fieldnames(Soft_results);

Soft_col = zeros(length(volume_fields),1);
for it_volume=1:length(volume_fields)
    stats = Soft_results.(volume_fields{it_volume}).stats;
    
%     if strcmp(volume_fields{it_volume}, 'Average_Amygdala')
%         pause(0.1)
%     end
    beta = stats.beta(3);
    gammaX=GAMMAX(3,3);
    sigma2 = sum(stats.resid.*stats.resid);
     
    [nout, tval, pval, power] = computeSampleSize(0, beta, sigma2, gammaX, N_COVARIATES, POWER, ALPHA);
    
    PA_Soft.full.(volume_fields{it_volume}).power = power;
    PA_Soft.full.(volume_fields{it_volume}).pval = pval;
    PA_Soft.full.(volume_fields{it_volume}).tval = tval;
    PA_Soft.full.(volume_fields{it_volume}).sample_size = nout;
    
%     if pval>ALPHA
%         disp(['Soft_ALPHA_' volume_fields{it_volume} ])
%     end
%     if power<POWER
%         disp(['Soft_POWER_' volume_fields{it_volume} num2str(power) ])
%     end
    
    
    Soft_col(it_volume) = nout; 
end




%% Plots

INTERESTING_VOLUMES = {'Cerebral_White_Matter', 'Cerebral_Cortex', 'Lateral_Ventricle', 'Thalamus', 'Caudate', 'Putamen', 'Pallidum', 'Hippocampus', 'Amygdala'};

INTERESTING_VOLUMES_PLOT = {'Cerebral-White-Matter', 'Cerebral-Cortex', 'Lateral-Ventricle', 'Thalamus', 'Caudate', 'Putamen', 'Pallidum', 'Hippocampus', 'Amygdala'};


% Interesting table
MRI_interesting_col = zeros(length(INTERESTING_VOLUMES), 1);
Soft_interesting_col = zeros(length(INTERESTING_VOLUMES), 1);
Hard_interesting_col = zeros(length(INTERESTING_VOLUMES), 1);
for it_interesting=1:length(INTERESTING_VOLUMES)
    MRI_interesting_col(it_interesting) = PA_MRI.full.(['Average_' INTERESTING_VOLUMES{it_interesting}]).sample_size;
    Soft_interesting_col(it_interesting) = PA_Soft.full.(['Average_' INTERESTING_VOLUMES{it_interesting}]).sample_size;
    Hard_interesting_col(it_interesting) = PA_Hard.full.(['Average_' INTERESTING_VOLUMES{it_interesting}]).sample_size;
end
Tint = table(MRI_interesting_col,Hard_interesting_col,Soft_interesting_col, 'RowNames', INTERESTING_VOLUMES );

figure
plot(MRI_interesting_col,'b'), hold on
plot(Hard_interesting_col,'r'), hold on
plot(Soft_interesting_col,'k'), hold on
set(gca, 'YScale', 'log')
ylabel('Sample size')
xticklabels(INTERESTING_VOLUMES_PLOT)
xtickangle(90)
legend('MRI', 'Hard', 'Soft')
title(['Power analysis: removal level ' num2str(removalLevel)])

% % Full table
% T = table(MRI_col,Hard_col,Soft_col, 'RowNames', volume_fields );
% 
% figure
% plot(MRI_col,'b'), hold on
% plot(Hard_col,'r'), hold on
% plot(Soft_col,'k'), hold on
% set(gca, 'YScale', 'log')
% ylabel('Sample size')
% legend('MRI', 'Hard', 'Soft')

