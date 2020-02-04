% Startup script for photo reconstruction code
%
% Adds paths and makes sure that freesurfer has been properly initialised.
% Plesae modify the defualt paths to match your system.
%
% H. Tregidgo h.tregidgo@ucl.ac.uk Jan 2020

%% Set base path

basePath_default = '/home/henry/Documents/Brain/UWphoto';

% compiled version of freesurfer containing a python packages directory
default_FS = '/home/henry/Documents/Brain/Freesurfer/freesurfer';
% version of freesurfer containing samseg capable of segmenting photos
default_DV = '/home/henry/Documents/Brain/Freesurfer/GitVersion/freesurfer';

if exist(basePath_default,'dir')
    codepath = fullfile(basePath_default,'code');
    addpath(genpath(codepath));
    
    setenv('PHOTO_RECON_HOME',basePath_default);
    
else
    
    error('Please set basePath_defualt to existing directory')
    
end


%% initialise freesurfer paths 
% Includes adding a dev path for where uncompiled FS files are held. This
% allows the python scripts in this version of python to be called using
% the precompiled python packages from a compiled dev version of
% Freesurfer.

[FREESURFER_HOME,FS_PYTHON_DEV] = UWphoto_function_checkBashEnvironment(...
    default_FS,default_DV);

%% add freesurfer matlab path

fs_matlab_path = fullfile(FREESURFER_HOME,'matlab');

pathCell = regexp(path, pathsep, 'split');
if ispc  % Windows is not case-sensitive
  onPath = any(strcmpi(fs_matlab_path, pathCell));
else
  onPath = any(strcmp(fs_matlab_path, pathCell));
end

if ~onPath
   
    addpath(genpath(fs_matlab_path));
    
end

%% add freesurfer bin path to bash path

bash_path = getenv('PATH');

bin_path = fullfile(FREESURFER_HOME,'bin');

pathCell = regexp(bash_path, pathsep, 'split');
if ispc  % Windows is not case-sensitive
  onPath = any(strcmpi(fs_matlab_path, pathCell));
else
  onPath = any(strcmp(fs_matlab_path, pathCell));
end


if ~onPath
   
    setenv('PATH',[getenv('PATH'),pathsep,bin_path]);
    
end
