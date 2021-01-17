clear; 
clc;

codepath = pwd;

pld_msec=[100,300,500,700,1000,1300,1700,2200,2700,3200,3500];   % in milisec
pld=[0.1,0.3,0.5,0.7,1.0,1.3,1.7,2.2,2.7,3.2,3.5];   % in sec

%%% Rotate CBF Maps

for i =1:length(pld_msec)
    fname = strcat('CALC_CBF_LD1500_PLD',num2str(pld_msec(i)),'.nii');

    V = niftiread(fname);
    whos V

    direction = [0 0 1];
    V_rot = imrotate3(V,45, direction);
    whos V_rot

    niftiwrite(V_rot, strcat('Rotated_CBF_map_GeomMask_NoSmooth_PLD',num2str(pld_msec(i)),'_PUMP200.nii'));

end

