% Script to get a list of the label numbers used in the segmentations of
% the photos and the MRIs. List needed for reference during manual
% segmentation. 

%% set directories

sft_dir = '/home/henry/Documents/Brain/UWphoto/Results_no_extra_slices/SAMSEG';

hrd_dir = '/home/henry/Documents/Brain/UWphoto/Results_hard_no_extra_slices/SAMSEG';

mri_dir = '/home/henry/Documents/Brain/UWphoto/FLAIR_Scan_Data/SAMSEG';

%% get files

dlist_sft = dir(fullfile(sft_dir,'*','*crisp*'));

dlist_hrd = dir(fullfile(hrd_dir,'*','*crisp*'));

dlist_mri = dir(fullfile(mri_dir,'N*','*crisp*'));

load('/home/henry/Documents/Brain/UWphoto/code/data/FreeviewLUTvaluesMatlab.mat',...
    'FreeviewLUTvalues')

%% go through files soft

labels_sft = cell(length(dlist_sft),1);

for il=1:length(dlist_sft)
    
    temp=MRIread(fullfile(dlist_sft(il).folder,dlist_sft(il).name));
    
    labels_sft{il} = unique(temp.vol(:));
    
end

labels_sft_total = unique(vertcat(labels_sft{:}));


%% go through files hard

labels_hrd = cell(length(dlist_hrd),1);

for il=1:length(dlist_hrd)
    
    temp=MRIread(fullfile(dlist_hrd(il).folder,dlist_hrd(il).name));
    
    labels_hrd{il} = unique(temp.vol(:));
    
end

labels_hrd_total = unique(vertcat(labels_hrd{:}));


%% go through files mri

labels_mri = cell(length(dlist_mri),1);

for il=1:length(dlist_mri)
    
    temp=MRIread(fullfile(dlist_mri(il).folder,dlist_mri(il).name));
    
    labels_mri{il} = unique(temp.vol(:));
    
end

labels_mri_total = unique(vertcat(labels_mri{:}));

%% get ready for output 

if ~isequal(labels_mri_total,labels_sft_total) || ...
         ~isequal(labels_mri_total,labels_hrd_total)
     warning('some reconstructions have different label lists')
end

FreeviewLUTnum2=zeros(size(FreeviewLUTvalues,1),1);

for il=1:size(FreeviewLUTvalues,1)
FreeviewLUTnum2(il)=str2double(FreeviewLUTvalues{il,1});
end

membs_ind = ismember(FreeviewLUTnum2,labels_hrd_total);

Labels_present = FreeviewLUTvalues(membs_ind,:);

%% write out labels 

writetable(cell2table(Labels_present),'UWphoto_samsegusedlabels.csv')
