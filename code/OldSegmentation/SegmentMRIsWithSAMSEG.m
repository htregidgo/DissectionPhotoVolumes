% script to run SAMSEG on the MRIs
% For each volume:
% First, it tries to register the SAMSEG atlas to the ex vivo MRI, and
% shows you the result with freeview. If you like it, then it uses this
% registration to initilize SAMSEG. If you don't like it, then it asks you
% to register the template directly
% After that, you need to use a modified version of SAMSEG as in
% SegmentPhotosWithSAMSEG, but without switching off the bias field
% correction. I know, this is superugly. This is why you should totally
% come up with something better, e.g., by adding personalized flags to SAMSEG,
% to a) bypass bias field correction and b) skipping affine registration
% and using a provided affinely registered atlas directly.
% In addition, I must say that the automatic registration almost never
% works, and you need to manually align in most cases. This is bad, we need
% something better and more reproducible. I was thinking that a possibility
% would be to kill the extracerebral tissue from the SAMSEG tamplate
% (Nellie can help you with that; she would know what semiautomated tools
% to use, and you don't need anything superperfect)
clear

inputDir='/autofs/cluster/vive/UW_photo_recon/FLAIR_Scan_Data/';
outputDir='/autofs/cluster/vive/UW_photo_recon/FLAIR_Scan_Data/SAMSEG/';
EXVIVO_SAMSEG_FS_DIR='/autofs/space/panamint_005/users/iglesias/fsdevsamsegexvivo/';
ALADIN='/autofs/space/panamint_005/users/iglesias/software/niftyreg-kcl/build/reg-apps/reg_aladin -ln 3 -lp 2 -speeeeed ';
nthread=16;

cases={'17-0333', % no missing photos in these ones
    '18-1045',
    '18-1132',
    '18-1196',
    '18-1680',
    '18-1913',
    '18-1930',
    '18-2128',
    '18-2259'};

if exist(outputDir,'dir')==0
    mkdir(outputDir);
end

% pass 1: register atlas and prepare SAMSEG dirs
T=MRIread([EXVIVO_SAMSEG_FS_DIR '/average/samseg/20Subjects_smoothing2_down2_smoothingForAffine2/template.orig.nii']);
atlasDir=[EXVIVO_SAMSEG_FS_DIR '/average/samseg/20Subjects_smoothing2_down2_smoothingForAffine2/'];
template=[tempdir '/template.nii.gz'];
MRIwrite(T,template);
T=MRIread(template);
for c=1:length(cases)
    
    disp(['Preparing grounds for case ' num2str(c) ' of ' num2str(length(cases))]);
    
    aux=cases{c};
    aux(aux=='-')='_';
    inputScan=[inputDir '/NP' aux '.rotated.nii.gz'];
    
    odir=[outputDir '/' cases{c} '/'];
    if exist(odir,'dir')==0
        mkdir(odir);
    end
    adir=[odir '/atlas/'];
    if exist(adir,'dir')==0
        mkdir(adir);
    end
    
    if exist([adir '/template.mgz'],'file')==0
        
        cmd=[ALADIN ' -ref ' inputScan ' -flo ' template ' -res /tmp/res.nii.gz -aff /tmp/aff.txt' ];
        system([cmd ' >/dev/null']);
        
        system(['freeview ' inputScan ' /tmp/res.nii.gz']);
        
        a=input('Do you like the output? ');
        
        if a>0
            
            fid = fopen('/tmp/aff.txt');
            M = fscanf(fid,'%f',[4,4])';
            fclose(fid);
            TT=[];
            TT.vol=T.vol;
            TT.vox2ras0=inv(M)*T.vox2ras0;
            TT.volres=sum(sqrt(TT.vox2ras0(1:3,1:3).^2));
            MRIwrite(TT,'/tmp/template.mgz');
            
        else
            system(['mri_convert ' template ' /tmp/template.mgz >/dev/null']);
            disp(['Align in freeview and save using only header']);
            system(['freeview ' inputScan ' /tmp/template.mgz']);
        end
        
        movefile('/tmp/template.mgz',adir);
    end
    
    % copy files to atlas directory
    copyfile([atlasDir '/atlas_level1.txt.gz' ],adir);
    copyfile([atlasDir '/atlas_level2.txt.gz' ],adir);
    copyfile([atlasDir '/compressionLookupTable.txt' ],adir);
    copyfile([atlasDir '/modifiedFreeSurferColorLUT.txt' ],adir);
    if exist([adir '/sharedGMMParameters.txt'],'file')==0
        a=input('Does the case have cerebellum and brainstem? ');
        if a > 0
            copyfile('./data/sharedGMMParameters.exvivoMRI.txt',[adir '/sharedGMMParameters.txt']);
        else
            copyfile('./data/sharedGMMParameters.NoBSorCB.exvivoMRI.txt',[adir '/sharedGMMParameters.txt']);
        end
    end
    
    % samseg call
    dd=dir([odir '/*crispSegmentation.nii']);
    if length(dd)==0
        cmd=['setenv FREESURFER_HOME ' EXVIVO_SAMSEG_FS_DIR ' && source $FREESURFER_HOME/SetUpFreeSurfer.csh &&  samseg '];
        cmd=[cmd ' --i ' inputScan];
        cmd=[cmd ' --o ' odir];
        cmd=[cmd ' --threads ' num2str(nthread)];
        cmd=[cmd ' --ssdd ' adir];
        cmd=[cmd ' --no-save-warp'];
        disp(cmd);
    else
        cmd=[];
    end
        
    cmds{c}=cmd;
    
end

scriptName=[pwd() '/script.sh'];
fid = fopen(scriptName,'w');
for j=1:length(cmds)
    if ~isempty(cmds{j})
        fprintf(fid,'%s\n',cmds{j});
    end
end
fclose(fid);
system(['chmod u+x ' scriptName]);
edit(scriptName);
disp(scriptName);



