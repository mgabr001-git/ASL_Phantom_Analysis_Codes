clear; 
clc;

codepath = pwd;

direction = [0 0 1];
% 
% V1 = niftiread('../Nifti/SUB01/3DSp_ubpCASL_TR6s_LD1500_PLD500_BS0_PUMP400_ASL_0018/meanPerf_SUB01.nii');
% V2 = niftiread('../Nifti/SUB01/3DSp_ubpCASL_TR6s_LD1500_PLD500_BS0_PUMP400_ASL_0018/STDERR_Perf_5pct_SUB01.nii');
% V3 = niftiread('../Nifti/SUB01/3DSp_ubpCASL_TR6s_LD1500_PLD500_BS0_PUMP400_ASL_0018/z_score_perf_5pct_SUB01.nii');
% V4 = niftiread('../Nifti/SUB01/3DSp_ubpCASL_TR6s_LD1500_PLD500_BS0_PUMP400_ASL_0018/RATIO_STDERRPerf_Perf_5pct_SUB01.nii');
% 
% M0 = niftiread('../Nifti/SUB01/3DSp_ubpCASL_TR6s_LD1500_PLD500_BS0_PUMP400_ASL_0018/meanM0_SUB01.nii');
% 
% 
% M0_rot = imrotate3(M0,45, direction);
% niftiwrite(M0_rot, './par_maps/results/NoSmooth_2M0_PUMP400_BS0/Rotated_Maps/Rotated_M0_PUMP400.nii')
% 
% 
%  V_rot = imrotate3(V1,45, direction);
%  niftiwrite(V_rot, './par_maps/results/NoSmooth_2M0_PUMP400_BS0/Rotated_Maps/meanPerf_SUB01.nii')
% 
% 
%   V_rot = imrotate3(V2,45, direction);
%  niftiwrite(V_rot, './par_maps/results/NoSmooth_2M0_PUMP400_BS0/Rotated_Maps/STDERR_Perf_5pct_SUB01.nii')
%  
%  
%   V_rot = imrotate3(V3,45, direction);
%  niftiwrite(V_rot, './par_maps/results/NoSmooth_2M0_PUMP400_BS0/Rotated_Maps/z_score_perf_5pct_SUB01.nii')
%  
%  
%   V_rot = imrotate3(V4,45, direction);
%  niftiwrite(V_rot, './par_maps/results/NoSmooth_2M0_PUMP400_BS0/Rotated_Maps/RATIO_STDERRPerf_Perf_5pct_SUB01.nii')
%  
%  
 
 
V1 = niftiread('../Nifti/SUB01/PCASL_Experiments_Set1/3DSp_ubpCASL_TR6s_LD1500_PLD500_ASYM_PUMP400_ASL_0080/meanPerf_SUB01.nii');

 V_rot = imrotate3(V1,45, direction);
 niftiwrite(V_rot, 'Rotated_ASYM_meanPerf_SUB01.nii')


 V2= niftiread('../Nifti/SUB01/3DSp_ubpCASL_TR6s_LD1500_PLD500_BS0_PUMP400_ASL_0103/meanPerf_SUB01.nii');
  V_rot = imrotate3(V2,45, direction);
 niftiwrite(V_rot, 'Rotated_SYMM_meanPerf_SUB01.nii')
 
 
 