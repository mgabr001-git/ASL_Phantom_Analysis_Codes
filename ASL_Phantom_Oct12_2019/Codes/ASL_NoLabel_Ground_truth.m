%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ASL analysis and fitting macro, based on % Perfusion Kinetic Model (Buxton et al 1998 )
%%%%% It requires DeltaM.m and ReturnTriplet.m
%%%%%
%%%%% Author: M. Gabrielyan

%%%%% Phantom Analysis Steps
%
% 1.  Convert_to_Nifti: Use NEW_Dicom2Nifti.m
% 2.  Create Geometric mask from an M0 image. Create ROI mask with respect to Geometric markers near the top of the phantom. 
% 3.  Run ASL_FullAnalysis_v3.m   to generate the perfusion maps, fit and  plot the data.
%
%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear; clc;
close all;

codepath = pwd;

cd ../Nifti
Niftipath = pwd;

cd(codepath);

A = dir(Niftipath); 
A = {A.name};
A(ismember(A,{'.','..','.DS_Store'})) = [];

% digits(10);

%%%% asl is a structure holding the constants used in the calculations.
%%%% passed as an argument to ReturnTriplet()

asl.pld_msec=[100,300,500,700,1000,1300,1700,2200,2700,3200,3500];   % in milisec
asl.pld=[0.1,0.3,0.5,0.7,1.0,1.3,1.7,2.2,2.7,3.2,3.5];   % in sec

asl.LD = 1500;        % Label Duration LD in msec
asl.TR =6;            % sec
asl.lambda = 1;
asl.T1blood =2.1;
asl.alpha=0.85;
asl.tau = asl.LD/1000;        %  tau = Label Duration LD in sec 

%%%% Choose ROI over which you want to average the signal
ROIflag = 2;  %%%  1 = Spher. Segment top , 2 = Full Phantom, 3 = Jet ROI 

unit_factor = 6000 ;   %%%   Flow unit conversion from mL/Vvox/sec into mL/100g/min 
smooth_pld=0.1:0.02:5.2;

sesno=length(asl.pld);

slc = 1:34;     % slice
col = 1:64;     % columns
row = 1:64;     % rows
   

perf_intensity=zeros(length(row),length(col),length(slc),sesno);
flow_intensity=zeros(length(row),length(col),length(slc),sesno);
M0_intensity=zeros(length(row),length(col),length(slc),sesno);
num_intensity=zeros(length(row),length(col),length(slc),sesno);
denom_intensity=zeros(length(row),length(col),length(slc),sesno);

mean_lbl=zeros(sesno,1);
mean_ctrl=zeros(sesno,1);
meanM0=zeros(sesno,1); 

%pump = [200, 300, 400] ;
pump = 0;

%  
% figure(5);
% set(gcf, 'Position',  [50, 900, 900, 550]);
% txt = {};  % Used for Legend generation
% mycol=['g','b','r','k','m'];
% Mean = zeros(1,sesno);
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%      Create ROI:                %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


R = 16      % Radius of the Spherical ROI
slice = 23
Oz_max = 33

%%% Reference Points

%%% B------C
%%% |      |
%%% |      |
%%% A------D


AA = [18 14 12];
BB = [18 44 12];
CC = [48 44 12];
DD = [48 14 12];


Ox = BB(1) + floor((CC(1)-BB(1))/2);
Oy = AA(2) + floor((BB(2)-AA(2))/2);
Oz = Oz_max - R ;

origin = [Ox Oy Oz]

roimsk=zeros(length(row),length(col),length(slc));


 for c=1:64
   for r=1:64
      for s=1:34
          
         vox=[r c s];
            roimsk(r, c, s) = (norm(vox - origin))^2 <= R^2  &    s > slice;   % hemisperical ROI
      end
   end
 end
   
 roipath = fullfile(Niftipath,A{1}) ;
 
%%%  To preserve the DICOM header info, with voxel size, coordinates  etc,
%%%  Read in an original nii volume. 

  dataFolder = fullfile(roipath,'3DSp_ubpCASL_TR6s_NoLabel_PUMP400_M0_0197/');
  firstfile = 's2019-02-17_12-27-142831-00001-00000-1.nii';

 P=spm_select('ExtFPlist',dataFolder,firstfile);
 v = spm_vol(P);   
 v0 = v(1); % Copy from volume v, in order to preserve the header info

 
 v0.fname = fullfile(strcat(roipath,'/SphereSegmentROI_R',num2str(R),'_Slice',num2str(slice),'_OzMax',num2str(Oz_max),'.nii'));  
 spm_write_vol( v0,roimsk);


 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%     Read in Jet ROI Mask generated from  M0 images
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    P=spm_select('ExtFPlist', roipath,'^JetROI_OnlyJets.nii');  % Pump rate ana 
%    JetMask = spm_read_vols(spm_vol(P));

   %  [-0.15,0.05]   //  initial pars for the fit
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%     Read in Geometric Mask generated from  M0 images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % P=spm_select('ExtFPlist', roipath,'^Phantom3mmROI_3Dsegment_LD1500_Pump200.nii');  % Pump rate ana 
 %  P=spm_select('ExtFPlist', roipath,'^PCASL_Exp_3Dsegment_LD1500_Pump400.nii');  % PCASL Experiments
   P=spm_select('ExtFPlist', roipath,'^Full_phantom_NoLabel_NoPump.nii');  % Pump rate ana 
   GeomMask1 = spm_read_vols(spm_vol(P));


   
   
   P=spm_select('ExtFPlist', roipath,'^Full_phantom_NoLabel_Pump400.nii');  % Pump rate ana 
   GeomMask2 = spm_read_vols(spm_vol(P));


 
%for pp = 1: length(pump) 
    
    
    
for sb = 1:length(A)    % Loop over the number of Subjects 
    subid = A{sb};
 

% Labeling is turned off
  aslsesdir_pumpOff = dir(fullfile(Niftipath,A{1},strcat('3DSp_ubpCASL_TR6s_NoLabel_PLD500_BS0_PUMP0_ASL*')));
  m0sesdir_pumpOff = dir(fullfile(Niftipath,A{1},strcat('3DSp_ubpCASL_TR6s_NoLabel_PLD500_BS0_PUMP0_M0*')));
  
  
  aslsesdir_pump300 = dir(fullfile(Niftipath,A{1},strcat('3DSp_ubpCASL_TR6s_NoLabel_PUMP300_ASL*')));
  m0sesdir_pump300 = dir(fullfile(Niftipath,A{1},strcat('3DSp_ubpCASL_TR6s_NoLabel_PUMP300_M0*')));
    

  aslsesdir_pump400 = dir(fullfile(Niftipath,A{1},strcat('3DSp_ubpCASL_TR6s_NoLabel_PUMP400_ASL*')));
  m0sesdir_pump400 = dir(fullfile(Niftipath,A{1},strcat('3DSp_ubpCASL_TR6s_NoLabel_PUMP400_M0*')));
  




%%%%%%  Run this after creating the Geometric Mask from M0 images as well
%%%%%%  as ROI mask with respect to Geometric Markers.

        [meanperf1, M01, weight1, Control1, Label1] = ReturnTriplet( aslsesdir_pumpOff, m0sesdir_pumpOff, subid,asl.pld_msec(1),asl,GeomMask1);
        
        
        
        [meanperf2, M02, weight2, Control2, Label2] = ReturnTriplet(aslsesdir_pump300, m0sesdir_pump300, subid,asl.pld_msec(1),asl,GeomMask2);
          
          
          
        [meanperf3, M03, weight3, Control3, Label3] = ReturnTriplet(aslsesdir_pump400, m0sesdir_pump400, subid,asl.pld_msec(1),asl,GeomMask2);
        

        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%  No smoothing, No Variance Weighteing 

         perf_pumpOff(:,:,:) = meanperf1(:,:,:)./M01(:,:,:) .* GeomMask1(:,:,:); 
         perf_pump300(:,:,:) = meanperf2(:,:,:)./M02(:,:,:) .* GeomMask2(:,:,:); 
         perf_pump400(:,:,:) = meanperf3(:,:,:)./M03(:,:,:) .* GeomMask2(:,:,:); 
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           

        perf_pumpOff(isnan(perf_pumpOff))=0;    % replace NaN-s with 0s
        perf_pump300(isnan(perf_pump300))=0;
        perf_pump400(isnan(perf_pump400))=0;
        
        
        
         msk = perf_pumpOff(:) ~= 0;
         mean(mean(mean(perf_pumpOff(msk))))
              
         msk = perf_pump300(:) ~= 0;
         mean(mean(mean(perf_pump300(msk))))
              
         msk = perf_pump400(:) ~= 0;
         mean(mean(mean(perf_pump400(msk))))
                         
        
           

lambda = asl.lambda;
T1b = asl.T1blood;
alpha = asl.alpha;
tau = asl.tau ;   
        
      for iip=1:sesno
              pld = asl.pld(iip); % convert to sec
        
             %%% Calculate CBF using White paper formula. 
  
              cbf_factor = 6000.0 .* lambda .* exp(pld ./ T1b) ./ ( 2 .* alpha  .* T1b .* (1 - exp(-tau ./ T1b))) ;
             
             
              Calc_CBF  = perf_pumpOff(:,:,:).* cbf_factor; 
              mask = Calc_CBF(:) ~= 0; 
              mean_CBF_pumpOff(iip) = mean(mean(mean(Calc_CBF(mask))));
              
            
             
              Calc_CBF  = perf_pump300(:,:,:).* cbf_factor;
              mask = Calc_CBF(:) ~= 0; 
              mean_CBF_pump300(iip) = mean(mean(mean(Calc_CBF(mask))));
              
            
              Calc_CBF3  = perf_pump400(:,:,:).* cbf_factor;
              mask = Calc_CBF3(:) ~= 0; 
              mean_CBF_pump400(iip) = mean(mean(mean(Calc_CBF3(mask))));
              
            
              
          
        
      end
        
      mean_CBF_pumpOff
      mean_CBF_pump300
      mean_CBF_pump400
      
     
end


