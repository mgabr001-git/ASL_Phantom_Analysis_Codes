clear; 
clc;

codepath = pwd;

% direction = [0 0 1];
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
%   V_rot = imrotate3(V1,45, direction);
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
 
 
 

%%% Rotate ATT Maps

V = niftiread('./par_maps/results/NoSmooth_2M0_PUMP400_BS0/ATT_map_GeomMask_NoSmooth_PUMP400.nii');
whos V

% subplot(2,1,1)
%imshow(V(31,:,:),[0, 3],'InitialMagnification',1000, 'Colormap', jet); 

direction = [0 0 1];
V_rot = imrotate3(V,45, direction);
whos V_rot

% subplot(2,1,2)
%imshow(V_rot(31,:,:),[0, 3],'InitialMagnification',1000, 'Colormap', jet); 

%volshow(V_rot)

niftiwrite(V_rot, 'Rotated_Maps/Rotated_ATT_map_GeomMask_NoSmooth_PUMP400.nii')



V = niftiread('./par_maps/results/NoSmooth_2M0_PUMP300_BS0/ATT_map_GeomMask_NoSmooth_PUMP300.nii');
V_rot = imrotate3(V,45, direction);
niftiwrite(V_rot, 'Rotated_Maps/Rotated_ATT_map_GeomMask_NoSmooth_PUMP300.nii')

V = niftiread('./par_maps/results/NoSmooth_2M0_PUMP200_BS0/ATT_map_GeomMask_NoSmooth_PUMP200.nii');
V_rot = imrotate3(V,45, direction);
niftiwrite(V_rot, 'Rotated_Maps/Rotated_ATT_map_GeomMask_NoSmooth_PUMP200.nii')




%GOF rotate

gofmap = niftiread('./par_maps/results/NoSmooth_2M0_PUMP400_BS0/GoodnessOfFit_R2_GeomMask_NoSmooth_PUMP400.nii');
gofmap_rot = imrotate3(gofmap,45, direction);
%whos gofmap_rot
niftiwrite(gofmap_rot, 'Rotated_Maps/Rotated_GoodnessOfFit_R2_GeomMask_NoSmooth_PUMP400.nii')

gofmap = niftiread('./par_maps/results/NoSmooth_2M0_PUMP300_BS0/GoodnessOfFit_R2_GeomMask_NoSmooth_PUMP300.nii');
gofmap_rot = imrotate3(gofmap,45, direction);
niftiwrite(gofmap_rot, 'Rotated_Maps/Rotated_GoodnessOfFit_R2_GeomMask_NoSmooth_PUMP300.nii')

gofmap = niftiread('./par_maps/results/NoSmooth_2M0_PUMP200_BS0/GoodnessOfFit_R2_GeomMask_NoSmooth_PUMP200.nii');
gofmap_rot = imrotate3(gofmap,45, direction);
niftiwrite(gofmap_rot, 'Rotated_Maps/Rotated_GoodnessOfFit_R2_GeomMask_NoSmooth_PUMP200.nii')





M0 = niftiread('../Nifti/SUB01/3DSp_ubpCASL_TR6s_LD1500_PLD500_BS0_PUMP400_ASL_0018/meanM0_SUB01.nii');
M0_rot = imrotate3(M0,45, direction);
niftiwrite(M0_rot, 'Rotated_Maps/Rotated_M0_PUMP400.nii')
 


M0 = niftiread('../Nifti/SUB01/3DSp_ubpCASL_TR6s_LD1500_PLD500_BS0_PUMP300_ASL_0058/meanM0_SUB01.nii');
M0_rot = imrotate3(M0,45, direction);
niftiwrite(M0_rot, 'Rotated_Maps/Rotated_M0_PUMP300.nii')
 

M0 = niftiread('../Nifti/SUB01/3DSp_ubpCASL_TR6s_LD1500_PLD500_BS0_PUMP200_ASL_0098/meanM0_SUB01.nii');
M0_rot = imrotate3(M0,45, direction);
niftiwrite(M0_rot, 'Rotated_Maps/Rotated_M0_PUMP200.nii')
 



