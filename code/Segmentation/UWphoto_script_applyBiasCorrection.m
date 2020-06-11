
%% Paths
Raw_volume_path = '17-0333.hard.recon.mgz';
bias_field_path = '17-0333.hard_uncorrected_gray.nii_biasField.nii';
grey_scale_path = '17-0333.hard_uncorrected_gray.nii.gz';
scal_factr_path = '17-0333.hard_uncorrected_gray.nii_scaling-factor.txt';
Sseg_bCrct_path = '17-0333.hard_uncorrected_gray.nii_biasCorrected.nii';

BC_out_path    = '17-0333.hard_brightnessCorrected.mgz';
BF_out_path    = '17-0333.hard_brightnessField.mgz';
Gry_out_path   = '17-0333.hard_greyscale.mgz';
BCgry_sseg_out_path = '17-0333.hard_greyscale_brightnessCorrected_sseg.mgz';
BCgry_out_path = '17-0333.hard_greyscale_brightnessCorrected.mgz';
BFcrct_out_path    = '17-0333.hard_brightnessField_adjusted.mgz';

%% read in data
scale_factor = dlmread(scal_factr_path);

vol_raw=MRIread(Raw_volume_path);

bias_field=MRIread(bias_field_path);


%% correct photo volume
vol_corrected=vol_raw;

for il=1:3
vol_corrected.vol(:,:,:,il)=squeeze(vol_raw.vol(:,:,:,il))./bias_field.vol...
    ./scale_factor;
end

MRIwrite( vol_corrected,BC_out_path);

%% correct greyscale ims

vol_sseg_BFcrct=MRIread(Sseg_bCrct_path);
MRIwrite( vol_sseg_BFcrct,BCgry_sseg_out_path);

vol_BFcrct=vol_sseg_BFcrct;
vol_BFcrct.vol=vol_sseg_BFcrct.vol./scale_factor;
MRIwrite( vol_BFcrct,BCgry_out_path);

vol_grey = MRIread(grey_scale_path);
MRIwrite( vol_grey,Gry_out_path);

%% correct biasfields


MRIwrite(bias_field,BF_out_path);

crct_bias = bias_field;
crct_bias.vol=bias_field.vol.*scale_factor;
MRIwrite(crct_bias,BFcrct_out_path);