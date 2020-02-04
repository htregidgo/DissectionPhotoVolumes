% Script to take the dice score values produced by
% UWphoto_script_generateDiceScores and plot them in a scatter plot for
% ease of comparison. 

%% set up directories

PHOTO_RECON_HOME=getenv('PHOTO_RECON_HOME');

path_label_table = fullfile(PHOTO_RECON_HOME,'code','data',...
    'UWphoto_samsegusedlabels.csv');

top_scores_dir = fullfile(PHOTO_RECON_HOME,'figures','diceScores');


%% find the dice files

dlist_soft_scores = dir(fullfile(top_scores_dir,'*soft_dice.mat'));

dlist_hard_scores = dir(fullfile(top_scores_dir,'*hard_dice.mat'));

%% read in values

full_soft_dice = [];
full_hard_dice = [];

for il=1:length(dlist_soft_scores)
    
    load(fullfile(dlist_soft_scores(il).folder,dlist_soft_scores(il).name),...
        'dice_scores')
    
    full_soft_dice=[full_soft_dice,dice_scores]; %#ok<AGROW>
    
end

for il=1:length(dlist_hard_scores)
    
    load(fullfile(dlist_hard_scores(il).folder,dlist_hard_scores(il).name),...
        'dice_scores')
    
    full_hard_dice=[full_hard_dice,dice_scores]; %#ok<AGROW>
    
end

table_labels = readtable(path_label_table);

%% setup plot 

label_present_soft = any(~isnan(full_soft_dice),2)& any(full_soft_dice>0,2);

label_present_hard = any(~isnan(full_hard_dice),2)& any(full_hard_dice>0,2);


label_values_soft = find(label_present_soft);

label_values_hard = find(label_present_hard);

soft_xticks = table_labels.Labels_present2(...
    ismember(table_labels.Labels_present1,label_values_soft));

hard_xticks = table_labels.Labels_present2(...
    ismember(table_labels.Labels_present1,label_values_hard));

%% make plots

plotx =repmat((1:length(soft_xticks))',[1,size(full_soft_dice,2)]);
ploty =full_soft_dice(label_present_soft,:);

h_soft = figure;
scatter(plotx(:),ploty(:))
set(gca,'xtick',1:length(soft_xticks),'xticklabel',soft_xticks)
xtickangle(60)
set(h_soft,'Position',[760,155,1280,720])
xlabel('Label')
ylabel('Dice score')
title('Soft reconstruction Dice scores')

drawnow
pause(2)

plotx =repmat((1:length(hard_xticks))',[1,size(full_hard_dice,2)]);
ploty =full_hard_dice(label_present_hard,:);

h_hard = figure;
scatter(plotx(:),ploty(:))
set(gca,'xtick',1:length(hard_xticks),'xticklabel',hard_xticks)
xtickangle(60)
set(h_hard,'Position',[760,155,1280,720])
xlabel('Label')
ylabel('Dice score')
title('Hard reconstruction Dice scores')

drawnow
pause(2)

%% save plots

saveas(h_soft,fullfile(top_scores_dir,'soft_diceScores.fig'))
saveas(h_soft,fullfile(top_scores_dir,'soft_diceScores.png'))

saveas(h_hard,fullfile(top_scores_dir,'hard_diceScores.fig'))
saveas(h_hard,fullfile(top_scores_dir,'hard_diceScores.png'))
