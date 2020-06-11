% Startup script for photo reconstruction code
%
% Adds paths and makes sure that freesurfer has been properly initialised.
% Plesae modify the defualt paths to match your system.
%
% Note this has only been tested on bash shells
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
    
    warning('Please set basePath_defualt to existing directory')
    warning('Now asking for base directory')
    
    basepath_guess = mfilename('fullpath');
    
    [basepath_guess,~,~] = fileparts(basepath_guess);
    [basepath_guess,~,~] = fileparts(basepath_guess);
    [basepath_guess,~,~] = fileparts(basepath_guess);
    
    
    basePath_default = uigetdir(basepath_guess,'Select UWPhoto directory');
    
    codepath = fullfile(basePath_default,'code');
    addpath(genpath(codepath));
    
    setenv('PHOTO_RECON_HOME',basePath_default);
    
end


%% initialise freesurfer paths
% Includes adding a dev path for where uncompiled FS files are held. This
% allows the python scripts in this version of python to be called using
% the precompiled python packages from a compiled dev version of
% Freesurfer.

if exist(default_FS,'dir') && exist(default_DV,'dir')
    [FREESURFER_HOME,FS_PYTHON_DEV] = UWphoto_function_checkBashEnvironment(...
        default_FS,default_DV);
else
    
    FREESURFER_HOME = getenv('FREESURFER_HOME');
    FS_PYTHON_DEV   = getenv('FS_PYTHON_DEV');
    
    warning('Please set default_FS to existing directory')
    warning('Please set default_DV to existing directory')
    
    if isempty(FREESURFER_HOME)
        warning('Now asking for freesurfer directories')
        
        default_FS = uigetdir(pwd,'Select Compiled Freesurfer directory');
        default_DV = uigetdir(pwd,'Select Compiled Freesurfer directory');
    else
        default_FS = FREESURFER_HOME;
        if isempty(FS_PYTHON_DEV)
            default_DV = default_FS;
        end
    end
    
    
    [FREESURFER_HOME,FS_PYTHON_DEV] = UWphoto_function_checkBashEnvironment(...
        default_FS,default_DV);
    
end

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
