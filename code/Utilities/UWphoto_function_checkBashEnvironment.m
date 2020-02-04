function [FREESURFER_HOME,FS_PYTHON_DEV] = UWphoto_function_checkBashEnvironment(varargin)
%UWPHOTO_FUNCTION_CHECKBASHENVIRONMENT Summary of this function goes here
%   Detailed explanation goes here


%% parse inputs

default_FS = '/home/henry/Documents/Brain/Freesurfer/freesurfer';
default_DV = '/home/henry/Documents/Brain/Freesurfer/GitVersion/freesurfer';

p=inputParser;

addOptional(p,'FS_HOME',default_FS,@(x) exist(x,'dir'))
addOptional(p,'FS_DEV',default_DV,@(x) exist(x,'dir'))
addOptional(p,'verbose',false,@(x) islogical(x) ||...
    (isnumeric(x) && (x==1 || x==0) ))

parse(p,varargin{:})

FS_HOME = p.Results.FS_HOME;
FS_DEV = p.Results.FS_DEV;
verbose = p.Results.verbose;

clearvars p


%% set pathss


FREESURFER_HOME = getenv('FREESURFER_HOME');
FS_PYTHON_DEV   = getenv('FS_PYTHON_DEV');

if isempty(FREESURFER_HOME) || (exist(FS_HOME,'dir') && ~strcmp(FREESURFER_HOME,FS_HOME))
    if exist(FS_HOME,'dir')
        setenv('FREESURFER_HOME',FS_HOME);
        FREESURFER_HOME = FS_HOME;
    else
        error('Please set FREESURFER_HOME environment variable using setenv.')
    end
elseif verbose 
    disp(FREESURFER_HOME)
end

if isempty(FS_PYTHON_DEV) || (exist(FS_DEV,'dir') && ~strcmp(FS_PYTHON_DEV,FS_DEV))
    if exist(FS_DEV,'dir')
        setenv('FS_PYTHON_DEV',FS_DEV);
        FS_PYTHON_DEV = FS_DEV;
    else
        error(['Please set FS_PYTHON_DEV environment variable using setenv.',...
            '\nThis should point to the python freesurfer codebase'])
    end
elseif verbose
    disp(FS_PYTHON_DEV)
end


end

