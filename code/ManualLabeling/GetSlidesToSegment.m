
%%%%%%%%%% Parameters
filepath = [pwd() filesep 'data' filesep 'segmentation_slices.csv'];
MRI_PATH='/home/acasamitjana/Results/UWPhoto_mod/hard/';
if ~exist(MRI_PATH,'dir')
    
    MRI_PATH = '/home/henry/Documents/Brain/UWphoto/Results_hard';
    
end
    
PHOTO_PATH='/home/acasamitjana/Data/UWphoto/Photo_data';
if ~exist(PHOTO_PATH,'dir')
    
    PHOTO_PATH = '/home/henry/Documents/Brain/UWphoto/Photo_data_updated';
    
end


% Set up these parameters according to the reconstruction used in freeivew
% to get init/fi slices.
Nphotos_pre = 2;
Nphotos_post = 2;


%%%%%%%%%% Read file
fileID = fopen(filepath,'r');
g = textscan(fileID,'%s','delimiter','\n');
num_lines = length(g{1}) - 1;
fclose(fileID);

fileID = fopen(filepath,'r');
fgetl(fileID);
for it_lines=1:num_lines
    tline = fgetl(fileID);
    tline = strsplit(tline, ',');

    subject = tline{1};
    init_slice = str2num(tline{2});
    fi_slice = str2num(tline{3});
    
    if length(tline) > 3  %it means that we have a fixed slice number
        num_slice = str2num(tline{4});
    else
        num_slice = randi([init_slice fi_slice]);
    end
    
    
    %Load slice from MRI.
    MRI = MRIread([MRI_PATH filesep subject filesep subject '.hard.recon.mgz']);
    
    inputPhotoDir = [PHOTO_PATH filesep subject filesep subject ' MATLAB'];
    d=dir([inputPhotoDir '/*.mat']);
    Nphotos=length(d);
    Nslices = Nphotos_pre;
    
    if exist(fullfile(PHOTO_PATH,subject,'slice_order.mat'),'file')
        
        load(fullfile(PHOTO_PATH,subject,'slice_order.mat'),'slice_order')
        
        num_slice = Nphotos_pre+slice_order(num_slice-Nphotos_pre);
        
    end
    
    %Load LABELS until you get the file you want.
    for n=1:Nphotos
        
        load([inputPhotoDir '/' d(n).name(1:end)],'LABELS'); Y=LABELS; clear LABELS
        unique_y = unique(Y);
        
        if num_slice <= Nslices + length(unique_y) - 1  %1 is for background
            clear Y; load([inputPhotoDir '/' d(n).name(1:end)],'LABELS'); Y=LABELS; clear LABELS
            
            X=imread([inputPhotoDir '/' d(n).name(1:end-4) '.tif']);
            
            label_slice = num_slice - Nslices;
        
            Y_new = uint8(255*(Y == label_slice));
            imwrite(Y_new, [MRI_PATH filesep subject filesep 'Mask_ImageFile' d(n).name(end-4:end-4) 'slice' num2str(num_slice) '.png'])
            Y_new = Y_new/255;
            X_new = zeros(size(Y), 'uint8');
            for c=1:3
                X_new(:, :, c) = X(:, :, c).*Y_new;
            end
            
            imwrite(X_new, [MRI_PATH filesep subject filesep 'Segmentation_ImageFile' d(n).name(end-4:end-4) 'slice' num2str(num_slice) '.png'])
        
            break
        else
            Nslices = Nslices + length(unique_y) - 1;
        end
               
    end    

end
fclose(fileID);