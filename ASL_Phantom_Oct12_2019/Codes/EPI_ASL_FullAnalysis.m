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

asl.pld_msec=[100,300,500,700,1000,1300,1700,2200,2700];   % in milisec
asl.pld=[0.1,0.3,0.5,0.7,1.0,1.3,1.7,2.2,2.7];

slice_times = [0.8525,      0, 0.9025, 0.0500, 0.9525, 0.1000, 1.0025, 0.1500, ...
               1.0525, 0.2000, 1.1025, 0.2525, 1.1525, 0.3025, 1.2025, 0.3525, ...
               1.2525, 0.4025, 1.3025, 0.4525, 1.3525, 0.5025, 1.4025, 0.5525, ...
               1.4525, 0.6025, 1.5025, 0.6525, 1.5525, 0.7025, 1.6025, 0.7525, ...
               1.6550, 0.8025];   % converted to sec

asl.LD = 1500;        % Label Duration LD in msec
asl.TR =6;            % sec
asl.lambda = 1;
asl.T1blood =2.1;
asl.alpha=0.85;
asl.tau = asl.LD/1000;        %  tau = Label Duration LD in sec 
asl.M0scale = 10; 

unit_factor = 6000;   %%%   (Flow unit conversion from mL/Vvol/sec into mL/100g/min  or fractional exchange per sec)
smooth_pld=0.1:0.02:5.2;

slidx = 28;


sesno=length(asl.pld);

slc = 1:34;     % slice
col = 1:64;     % columns
row = 1:64;     % rows
   

perf_intensity=zeros(length(row),length(col),length(slc),sesno);
flow_intensity=zeros(length(row),length(col),length(slc),sesno);
M0_intensity=zeros(length(row),length(col),length(slc),sesno);
num_intensity=zeros(length(row),length(col),length(slc),sesno);
denom_intensity=zeros(length(row),length(col),length(slc),sesno);

 
 pump = [200, 300, 400] ;

% pump = 400;

 
figure(5);
set(gcf, 'Position',  [50, 900, 900, 550]);
txt = {};  % Used for Legend generation
mycol=['g','b','r','k','m'];
Mean = zeros(1,sesno);
ROIflag = 1;  %%%  1 = Conical ROI, 0 = Thresholded
threshval = 0.03;   %%% Threshold Value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%      Create ROI:   Conical base + Spherical top                %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Optimized for Non-smooth case (SPIRAL)
R = 16      % Radius of the Spherical ROI
angle = 60  % Opening of the conus angle 
slice = 23
Oz_max = 33


% %%%%   Optimized for EPI 
% R = 13     % Radius of the Spherical ROI
% angle = 60 % Opening of the conus angle 
% slice = 23
% Oz_max = 32
% 




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
        % roimsk(r, c, s) = (norm(vox - origin))^2 <= R^2  &   (r - Ox)^2+(c - Oy)^2 <= (s-Oz)^2 * tan(deg2rad(angle)) & s > slice ;
        roimsk(r, c, s) = (norm(vox - origin))^2 <= R^2  & s > slice ; 
      end
   end
 end
   
 roipath = fullfile(Niftipath,A{1}) ;
 
%%%  To preserve the DICOM header info, with voxel size, coordinates  etc,
%%%  Read in an original nii volume. 

  dataFolder = fullfile(roipath,'ep2d_pcasl_TR6_LD1500_PLD100_PUMP200_0002/');
  firstfile = 'f2019-10-11_17-49-180106-00001-00001-1.nii';

 P=spm_select('ExtFPlist',dataFolder,firstfile);
 v = spm_vol(P);   
 v0 = v(1); % Copy from volume v, in order to preserve the header info

 
 v0.fname = fullfile(strcat(roipath,'/SphereSegmentROI_R',num2str(R),'Slice',num2str(slice),'OzMax',num2str(Oz_max),'.nii'));  
 spm_write_vol( v0,roimsk);

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   Read in Thresholded ROI:  need to create first if it does not exist          %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   After creating the thresholded mask at highest pump rate,
%%%%%%%%%%%%%%%   uncomment this section to read in the mask file and apply to all PUMP rates and PLDs.
%%%%%%%%%%%%%%%   Comment out the creation part in lines 295-308. Follow
%%%%%%%%%%%%%%%   the instructions therein.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    file = strcat('ThresholdedROI_pump400_PLDvol6_th',num2str(threshval),'.nii');
%    P=spm_select('ExtFPlist', roipath,file);  % Pump rate ana 
%    thresholdmsk = spm_read_vols(spm_vol(P));
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%     Read in Geometric Mask generated from  M0 images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % P=spm_select('ExtFPlist', roipath,'^Phantom3mmROI_3Dsegment_LD1500_Pump200.nii');  % Pump rate ana 
 %  P=spm_select('ExtFPlist', roipath,'^PCASL_Exp_3Dsegment_LD1500_Pump400.nii');  % PCASL Experiments
 %  P=spm_select('ExtFPlist', roipath,'^FullPhantom_pump200.nii');  % Pump rate ana 
   
   P=spm_select('ExtFPlist', roipath,'^EPI_FullPhantom_pump400.nii');  % Pump rate ana 
   
   GeomMask = spm_read_vols(spm_vol(P));


 
for pp = 1: length(pump) 
    
    
    
for sb = 1:length(A)    % Loop over the number of Subjects 
    subid = A{sb};
 
    for i =1:sesno
   
%%%% Labeled Set
% 
%         aslsesdir = dir(fullfile(Niftipath,A{sb},strcat('3DSp_ubpCASL_TR6s_LD',num2str(asl.LD),'_PLD',num2str(asl.pld_msec(i)),'_BS0_PUMP',num2str(pump(pp)),'_ASL*')));
%         m0sesdir = dir(fullfile(Niftipath,A{sb},strcat('3DSp_ubpCASL_TR6s_LD',num2str(asl.LD),'_PLD',num2str(asl.pld_msec(i)),'_BS0_PUMP',num2str(pump(pp)),'_M0*')));
%      
          
%%%%%%%%%   EPI   Danny's sequence

        aslsesdir = dir(fullfile(Niftipath,A{sb},strcat('ep2d_pcasl_TR6_LD',num2str(asl.LD),'_PLD',num2str(asl.pld_msec(i)),'_PUMP',num2str(pump(pp)),'_*')));
       % m0sesdir = dir(fullfile(Niftipath,A{sb},strcat('ep2d_pcasl_TR6_NoLabel_BS0_PUMP',num2str(pump(pp)),'_*')));   % using Unlabeled ASL images as M0 for each pump rate.

  



%%%%%%  Run this after creating the Geometric Mask from M0 images as well
%%%%%%  as ROI mask with respect to Geometric Markers.
        
        [meanperf,M0,meanCBF,weight] = EPI_ReturnTriplet(aslsesdir,subid,asl.pld_msec(i),asl); 
      
        
        numerator = meanperf .* weight; 
        denominator = weight; 
        
        % Clean up data
        numerator(isnan(numerator))=0;
        denominator(isnan(denominator))=0;
        numerator(isinf(numerator))=0;
        denominator(isinf(denominator))=0;
        
        perf_intensity(:,:,:,i) = meanperf(:,:,:);
        M0_intensity(:,:,:,i) = M0(:,:,:);
        num_intensity(:,:,:,i) = numerator(:,:,:);
        denom_intensity(:,:,:,i) = denominator(:,:,:);
        
    end
       
end

        

 sigma = 3.0 ;

% DM_map=zeros(length(r),length(c),length(s),sesno);
 
 M0_filterg3=zeros(length(row),length(col),length(slc),sesno);
 num_filterg3=zeros(length(row),length(col),length(slc),sesno);
 denom_filterg3=zeros(length(row),length(col),length(slc),sesno);
 weighted_perf=zeros(length(row),length(col),length(slc),sesno);


for iip=1:sesno
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%      %%%%%%  SMOOTHING.    Uncomment this section if want to smooth and varicne-weight the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           M0_filterg3(:,:,:,iip) = imgaussfilt3(squeeze(M0_intensity(:,:,:,iip)),sigma);         
%           num_filterg3(:,:,:,iip) = imgaussfilt3(squeeze(num_intensity(:,:,:,iip)),sigma);           
%           denom_filterg3(:,:,:,iip) = imgaussfilt3(squeeze(denom_intensity(:,:,:,iip)),sigma);
% 
%      %%%%%% Variance weighted and SMOOTHED perfusion .  (perf .* Weight .* GeomMask) * K / (Weight .* GeomMask) * K     (K is the Gaussian convolution kernel with 3D isotropic sigma) 
%         
%         weighted_perf(:,:,:,iip) = (num_filterg3(:,:,:,iip) ./ denom_filterg3(:,:,:,iip) ./ M0_filterg3(:,:,:,iip) ) .* GeomMask(:,:,:); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%  No smoothing, No Variance Weighteing 
         weighted_perf(:,:,:,iip) = perf_intensity(:,:,:,iip)./M0_intensity(:,:,:,iip) .* GeomMask(:,:,:); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
end
        weighted_perf(isnan(weighted_perf))=0;    % replace NaN-s with 0s

  

for  iip=1:sesno

       
       
%%%%      NO SMOOTHING  %%%%
%        v0.fname = fullfile(strcat('par_maps/results/EPI_NoSmooth_2M0_PUMP',num2str(pump(pp)),'_BS0/'),strcat('EPI_NoSmooth_Perf_LD',num2str(asl.LD),'_PLD',num2str(asl.pld_msec(iip)),'.nii'));         
%        spm_write_vol( v0,squeeze(weighted_perf(:,:,:,iip)) );
%        


  

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%   Thresholded ROI, Specific Value   %%%%%%%%%%%%%%%%
%%%%%% Create this mask at highest pump rate and at PLD volume where you expect 
%%%%%%  the Max signal (in this case Pump=400 mL/min, PLD vol 6 = 1300 ms )
%%%%%% Just run this code once then comment out this section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  thresholdmsk = weighted_perf(:,:,:,6) >= threshval;
%      
%  v0.fname = fullfile(strcat(roipath,'/ThresholdedROI_pump400_PLDvol6_th',num2str(threshval),'.nii'));         
%  spm_write_vol( v0,thresholdmsk);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Apply  ROI to calculate the Mean and STDev.
%%% ROIflag:  1 = Conical ROI, 0 = Thresholded
   
slpld = [];
slmean = [];



        im = squeeze(weighted_perf(:,:,:,iip));
                 
         for isl= 1:34
           tempsl = im(:,:,isl);
   
           slmean = [slmean, mean(tempsl(tempsl~=0))];
           Slice_mean(isl,iip) = mean(tempsl(tempsl~=0));
           
           pld_sl = asl.pld(iip) + slice_times(isl);
           slpld = [slpld, pld_sl] ;
           
         end

    
       
       sl_pld(iip) = asl.pld(iip) + slice_times(slidx) ;
       Mean(iip) = slmean(slidx);
            
       
       
         
end

%   figure(2);
% 
%   subplot(1,3,pp);
%  plot(sl_pld,sl_mean); 



  eval(sprintf('EPI_sig_pump%d = [Mean]', pump(pp))) ;
 % eval(sprintf('Oct12_err_pump%d = [STDev]', pump(pp))) ;
 % eval(sprintf('Oct12_wgt_pump%d = [1./(Oct12_err_pump%d .^2)]', pump(pp),pump(pp))) ;


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% dM/M0 vs PLD for each Pump Rate: OVERLAY
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5);

  
  sig = eval(sprintf('EPI_sig_pump%d', pump(pp)));
 
   ft = fittype('DeltaM(x,lambda,tau,T1blood,p1,p2)', 'problem', {'lambda','tau','T1blood'},'coefficients',{'p1','p2'});
    [ff, gof]= fit((asl.pld+slice_times(slidx)+asl.tau)',(sig)',ft,'StartPoint',[2.,0.04],'problem',{asl.lambda,asl.tau,asl.T1blood});
   
   %%%% Use Fit parameters to plot Buxton Model. 
   DM =DeltaM(smooth_pld,asl.lambda,asl.tau,asl.T1blood,ff.p1,ff.p2);
  h=plot(asl.pld + asl.tau  ,sig,'+','Color',mycol(pp));  
    
  % h=plot(asl.pld + asl.tau + slice_times(fit_sliceN),sig,'+','Color',mycol(pp));  
  
      
   tt = sprintf('Data: PR = %d mL/min',pump(pp));
   txt = [txt; tt];
   
   hold on 
  hf = plot(smooth_pld-slice_times(slidx), DM,'LineWidth',2,'LineStyle', '--', 'Color', mycol(pp)); % plot best fit curve, shifted x-axis values
 
 % hf = plot(smooth_pld - slice_times(fit_sliceN), DM,'LineWidth',2,'LineStyle', '--', 'Color', mycol(pp)); % plot best fit curve, shifted x-axis values

  
%    tt = sprintf('BM: PR = %d mL/min',pump(pp));
%    txt = [txt; tt];
   hf.Annotation.LegendInformation.IconDisplayStyle = 'off';
   
   
   dt = sprintf('%4.2f',ff.p1+asl.tau - slice_times(slidx) )
  flow = sprintf('%4.2f',ff.p2*unit_factor);
    str = {strcat('dt=', num2str(dt)),strcat('flow=', num2str(flow))};
  % text(0.25, 0.01 + pp*0.01, str,'Color',mycol(pp),'FontSize',14);
%   
  set(gca,'FontSize',24);   % axis font size
  set(h,'LineWidth',2,'MarkerSize',13);

      
  hold on
  

end

hold off
  % title('3 mm Phantom, BS = 0','FontSize',16);
   xlabel('Time (s)');
   ylabel('\Delta M /M0');
   xlim([0 5.5]);
   ylim([0 0.05]);
 
   y1=get(gca,'ylim');
   LL=line([1.5 1.5], y1);
   set(LL,'LineWidth',2);
 legend(txt);
 hold off
% 
%  if ROIflag == 0 
%    fname=fullfile(codepath,strcat('par_maps/results/FitResults/ConicalBase_Reproduce_EPI_Oct12_data/Thresholded/OVERLAY_EPI_Perf_vs_PLD_th',num2str(threshval)));
%     print(strcat(fname),'-djpeg');  
%  else  
%     fname=fullfile(codepath,strcat('par_maps/results/FitResults/ConicalBase_Reproduce_EPI_Oct12_data/OVERLAY_EPI_Perf_vs_PLD_R',num2str(R),'The',num2str(angle),'Slice',num2str(slice),'OzMax',num2str(Oz_max)));
%     print(strcat(fname),'-djpeg');  
%      
%  end 



    fname=fullfile(codepath,strcat('OVERLAY_EPI_Slice',num2str(slidx),'_Perf_vs_PLD'));
    print(strcat(fname),'-djpeg');  




    