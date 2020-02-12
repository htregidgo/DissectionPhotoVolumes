% Script to generate correlations between volumes from segmentations on
% both photo reconstructions and the MRI used as ground truth. The script
% relies upon saved .mat files from the script UWphoto_script_extractStats.

clearvars

%% load data 

PHOTO_RECON_HOME=getenv('PHOTO_RECON_HOME');

figuresDir = fullfile(PHOTO_RECON_HOME,'figures');
corltnsDir = fullfile(figuresDir,'correlations');

if ~exist(corltnsDir,'dir')
    mkdir(corltnsDir)
end

load(fullfile(figuresDir,'AdjustedCaseStats.mat'),'segVolumeInfo')

strfields = fieldnames(segVolumeInfo);

softdir = fullfile(corltnsDir,'soft');
if ~exist(softdir,'dir')
    mkdir(softdir)
end

harddir = fullfile(corltnsDir,'hard');
if ~exist(harddir,'dir')
    mkdir(harddir)
end

crctdir = fullfile(corltnsDir,'hardCorrected');
if ~exist(crctdir,'dir')
    mkdir(crctdir)
end


%% setup variables

mriFlag = ismember({segVolumeInfo.segtype},'MRI');
N_mri   = sum(mriFlag);

hrdFlag = ismember({segVolumeInfo.segtype},'Hard');
N_hrd   = sum(hrdFlag);

sftFlag = ismember({segVolumeInfo.segtype},'Soft');
N_sft   = sum(sftFlag);


CorrelationStruct(length(strfields)-2) = struct();

%% iterate through labels

for il=3:length(strfields)
        
    strctI = il-2;
    
    MRIsplot = [segVolumeInfo(mriFlag).(strfields{il})];
    hardplot = [segVolumeInfo(hrdFlag).(strfields{il})];
    hardCrct = [segVolumeInfo(hrdFlag).(strfields{il})]...
        .*[segVolumeInfo(hrdFlag).volumeAdjustmentFactor];
    softplot = [segVolumeInfo(sftFlag).(strfields{il})];
    
    
    CorrelationStruct(strctI).Label = strfields{il};
    CorrelationStruct(strctI).MRI.raw = MRIsplot;
    
    %% soft statistics
    
    CorrelationStruct(strctI).soft.raw = softplot;
    
    [rho,p,rhoLower,rhoUpper]=corrcoef(MRIsplot,softplot);
    
    CorrelationStruct(strctI).soft.corrcoef.rho = rho(1,2);
    CorrelationStruct(strctI).soft.corrcoef.p = p(1,2);
    CorrelationStruct(strctI).soft.corrcoef.rhoLower = rhoLower(1,2);
    CorrelationStruct(strctI).soft.corrcoef.rhoUpper = rhoUpper(1,2);
    
    [b,fitstats] = robustfit(MRIsplot,softplot);
    
    CorrelationStruct(strctI).soft.robustfit.b = b;
    CorrelationStruct(strctI).soft.robustfit.stats = fitstats;
    
    
    %% plot figures
    
    fig_soft=figure;
    plot(MRIsplot,softplot,'*')
    hold on
    plot(MRIsplot,b(1)+b(2)*MRIsplot,'r','LineWidth',2)
    
    pltTitle = strrep(['Soft: ', strfields{il}],'_',' ');
    
    title(pltTitle,'Interpreter','none')
    
    str_r=[' r= ',num2str(rho(1,2)),', p= ',num2str(p(1,2))];
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str_r);
    set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    
    xlabel('MRI label volume')
    ylabel('Photo label volume')
    
    
    %% hard statistics
    
    CorrelationStruct(strctI).soft.raw = hardplot;
    
    [rho,p,rhoLower,rhoUpper]=corrcoef(MRIsplot,hardplot);
    
    CorrelationStruct(strctI).hard.corrcoef.rho = rho(1,2);
    CorrelationStruct(strctI).hard.corrcoef.p = p(1,2);
    CorrelationStruct(strctI).hard.corrcoef.rhoLower = rhoLower(1,2);
    CorrelationStruct(strctI).hard.corrcoef.rhoUpper = rhoUpper(1,2);
    
    [b,fitstats] = robustfit(MRIsplot,hardplot);
    
    CorrelationStruct(strctI).hard.robustfit.b = b;
    CorrelationStruct(strctI).hard.robustfit.stats = fitstats;
    
    
    %% plot figures
    
    fig_hard=figure;
    plot(MRIsplot,hardplot,'*')
    hold on
    plot(MRIsplot,b(1)+b(2)*MRIsplot,'r','LineWidth',2)
    
    pltTitle = strrep(['Hard: ', strfields{il}],'_',' ');
    
    title(pltTitle,'Interpreter','none')
    
    str_r=[' r= ',num2str(rho(1,2)),', p= ',num2str(p(1,2))];
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str_r);
    set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    
    xlabel('MRI label volume')
    ylabel('Photo label volume')
    
    
    %% volume corrected hard statistics
    
    CorrelationStruct(strctI).soft.raw = hardCrct;
    
    [rho,p,rhoLower,rhoUpper]=corrcoef(MRIsplot,hardCrct);
    
    CorrelationStruct(strctI).hardCrctd.corrcoef.rho = rho(1,2);
    CorrelationStruct(strctI).hardCrctd.corrcoef.p = p(1,2);
    CorrelationStruct(strctI).hardCrctd.corrcoef.rhoLower = rhoLower(1,2);
    CorrelationStruct(strctI).hardCrctd.corrcoef.rhoUpper = rhoUpper(1,2);
    
    [b,fitstats] = robustfit(MRIsplot,hardCrct);
    
    CorrelationStruct(strctI).hardCrctd.robustfit.b = b;
    CorrelationStruct(strctI).hardCrctd.robustfit.stats = fitstats;
    
    
    %% plot figures
    
    fig_crct=figure;
    plot(MRIsplot,hardCrct,'*')
    hold on
    plot(MRIsplot,b(1)+b(2)*MRIsplot,'r','LineWidth',2)
    
    pltTitle = strrep(['Hard: ', strfields{il}],'_',' ');
    
    title(pltTitle,'Interpreter','none')
    
    str_r=[' r= ',num2str(rho(1,2)),', p= ',num2str(p(1,2))];
    T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str_r);
    set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
    
    xlabel('MRI label volume')
    ylabel('Photo label volume')
    
    
    %% save figures
    
    saveas(fig_soft,fullfile(softdir,strfields{il}))
    saveas(fig_soft,fullfile(softdir,[strfields{il},'.tiff']))
    
    saveas(fig_hard,fullfile(harddir,strfields{il}))
    saveas(fig_hard,fullfile(harddir,[strfields{il},'.tiff']))
    
    saveas(fig_crct,fullfile(crctdir,strfields{il}))
    saveas(fig_crct,fullfile(crctdir,[strfields{il},'.tiff']))
    
    close all


end


%% save correlations statistics

save(fullfile(corltnsDir,'CorrelationStructure'),'CorrelationStruct')


