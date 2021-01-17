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


slidx = 27;


unit_factor = 6000;   %%%  Flow unit conversion from mL/Vvol/sec into mL/100g/min
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

 
 pump = [200, 300, 400] ;

% pump = 400;

 
figure(5);
set(gcf, 'Position',  [50, 900, 900, 550]);
txt = {};  % Used for Legend generation
txt1 = {};
mycol=['g','b','r','k','m'];
Mean = zeros(1,sesno);
ROIflag = 1;  %%%  1 = Conical ROI, 0 = Thresholded


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
          roimsk(r, c, s) = (norm(vox - origin))^2 <= R^2  & s > slice ; 
      end
   end
 end
   
 roipath = fullfile(Niftipath,A{1}) ;
 
%%%  To preserve the DICOM header info, with voxel size, coordinates  etc,
%%%  Read in an original nii volume. 

  dataFolder = fullfile(roipath,'3DSp_ubpCASL_TR6s_LD1500_PLD100_BS0_PUMP200_ASL_0092/');
  firstfile = 's2019-10-12_13-25-175425-00001-00001-1.nii';

 P=spm_select('ExtFPlist',dataFolder,firstfile);
 v = spm_vol(P);   
 v0 = v(1); % Copy from volume v, in order to preserve the header info

 
 v0.fname = fullfile(strcat(roipath,'/SphericalSegmentROI_R',num2str(R),'Slice',num2str(slice),'OzMax',num2str(Oz_max),'.nii'));  
 spm_write_vol( v0,roimsk);

 

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%     Read in Geometric Mask generated from  M0 images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % P=spm_select('ExtFPlist', roipath,'^Phantom3mmROI_3Dsegment_LD1500_Pump200.nii');  % Pump rate ana 
 %  P=spm_select('ExtFPlist', roipath,'^PCASL_Exp_3Dsegment_LD1500_Pump400.nii');  % PCASL Experiments
   P=spm_select('ExtFPlist', roipath,'^FullPhantom_pump200.nii');  % Pump rate ana 
   GeomMask = spm_read_vols(spm_vol(P));


 
for pp = 1: length(pump) 
    
    
    
for sb = 1:length(A)    % Loop over the number of Subjects 
    subid = A{sb};
 
    for i =1:sesno
   
%%%% Labeled Set

        aslsesdir = dir(fullfile(Niftipath,A{sb},strcat('3DSp_ubpCASL_TR6s_LD',num2str(asl.LD),'_PLD',num2str(asl.pld_msec(i)),'_BS0_PUMP',num2str(pump(pp)),'_ASL*')));
        m0sesdir = dir(fullfile(Niftipath,A{sb},strcat('3DSp_ubpCASL_TR6s_LD',num2str(asl.LD),'_PLD',num2str(asl.pld_msec(i)),'_BS0_PUMP',num2str(pump(pp)),'_M0*')));
     
            
     

%%%%%%  Run this after creating the Geometric Mask from M0 images as well
%%%%%%  as ROI mask with respect to Geometric Markers.

        [meanperf,M0,meanCBF,weight] = ReturnTriplet(aslsesdir,m0sesdir,subid,asl.pld_msec(i),asl,GeomMask);
        
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
  Calc_CBF= zeros(length(row),length(col),length(slc),sesno);

for iip=1:sesno

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%  No smoothing, No Variance Weighteing 
         weighted_perf(:,:,:,iip) = perf_intensity(:,:,:,iip)./M0_intensity(:,:,:,iip) .* GeomMask(:,:,:); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
end
        weighted_perf(isnan(weighted_perf))=0;    % replace NaN-s with 0s

  
       % To preserve the DICOM header info, with voxel size, coordinates  etc, read in an original nii volume       
 
        P=spm_select('ExtFPlist',dataFolder,['meanPerf_' A{1} '.nii']);
        v = spm_vol(P); 
        
        v0 = v(1);   % Copy from volume v, in order to preserve the header info: meanperf datatype is converted to 64 bit float in ReturnTriplet.m
         
        
        

for  iip=1:sesno

    
pld = asl.pld(iip); % convert to sec
lambda = asl.lambda;
T1b = asl.T1blood;
alpha = asl.alpha;
tau = asl.tau ;    

 cbf_factor = 6000.0 .* lambda .* exp(pld ./ T1b) ./ ( 2 .* alpha  .* T1b .* (1 - exp(-tau ./ T1b))) ;


%%%%      NO SMOOTHING, No Weighting  %%%%
%        v0.fname = fullfile(strcat('par_maps/results/NoSmooth_2M0_PUMP',num2str(pump(pp)),'_BS0/'),strcat('NoSmooth_Perf_LD',num2str(asl.LD),'_PLD',num2str(asl.pld_msec(iip)),'.nii'));         
%        spm_write_vol( v0,squeeze(weighted_perf(:,:,:,iip)) );


slpld = [];
slmean = [];


        im = squeeze(weighted_perf(:,:,:,iip));

    
     
         for isl= 1:34
           tempsl = im(:,:,isl);
%            mean(tempsl(tempsl~=0)); % slice_mean for the fit
   
           slmean = [slmean, mean(tempsl(tempsl~=0))];
           
           Slice_mean(isl,iip) = mean(tempsl(tempsl~=0));
           
           
         end

    
       
       
       Mean(iip) = slmean(slidx);
            
       

%%%%      NO SMOOTHING, No Weighting  %%%%
%        v0.fname = fullfile(strcat('par_maps/results/NoSmooth_2M0_PUMP',num2str(pump(pp)),'_BS0/'),strcat('CALC_CBF_LD',num2str(asl.LD),'_PLD',num2str(asl.pld_msec(iip)),'.nii'));         
%        spm_write_vol( v0,squeeze(Calc_CBF(:,:,:,iip)) );
% 
%        
       
       
end

% size(Slice_mean)
%  Slice_mean 



  eval(sprintf('Oct12_sig_pump%d = [Mean]', pump(pp))) ;
 % eval(sprintf('Oct12_err_pump%d = [STDev]', pump(pp))) ;
 % eval(sprintf('Oct12_wgt_pump%d = [1./(Oct12_err_pump%d .^2)]', pump(pp),pump(pp))) ;


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% dM/M0 vs PLD for each Pump Rate: OVERLAY
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
%  figure(4);
%  set(gcf, 'Position',  [850, 1000, 1000, 850]);
%  
%  subplot(3,1,pp)
%  
%   for isl = 25:2:29
%   sig = Slice_mean(isl,:);
%  
%    ft = fittype('DeltaM(x,lambda,tau,T1blood,p1,p2)', 'problem', {'lambda','tau','T1blood'},'coefficients',{'p1','p2'}); 
%    [ff, gof]= fit((asl.pld+asl.tau)',(sig)',ft,'StartPoint',[1.0,0.05],'problem',{asl.lambda,asl.tau,asl.T1blood});
%    gof
%    
%    %%%% Use Fit parameters to plot Buxton Model. 
%    DM =DeltaM(smooth_pld,asl.lambda,asl.tau,asl.T1blood,ff.p1,ff.p2);
%    h=plot(asl.pld+asl.tau,sig,'+','Color',[1-0.15*(isl-25) 0 0.6+0.1*(isl-25)]);  
%    tt = sprintf('Slice: %d',isl);
%    txt1 = [txt1; tt];
%   
%    hold on 
%    hf=plot(smooth_pld, DM,'LineWidth',2,'LineStyle', '--', 'Color', [1-0.15*(isl-25) 0 0.6+0.1*(isl-25)]); % plot best fit curve, shifted x-axis values
% %     tt = sprintf('BM: PR = %d mL/min',pump(pp));
% %     txt = [txt; tt];
%   hf.Annotation.LegendInformation.IconDisplayStyle = 'off';
%    dt = sprintf('%4.2f',ff.p1+asl.tau)
%   flow = sprintf('%4.2f',ff.p2*unit_factor) ;
%     str = {strcat('dt=', num2str(dt)),strcat('flow=', num2str(flow))};
%    text(0.25, 0.02 + (pp-1)*0.02, strcat('Pump Rate = ' , num2str(pump(pp)),' mL/min'),'Color','k','FontSize',18);
% %   
%   set(gca,'FontSize',24);   % axis font size
%   set(h,'LineWidth',2,'MarkerSize',13);
% 
% 
%   
%   hold on
%   
%   end
%  
% % hold off
% %   title('3 mm Phantom, BS = 0','FontSize',16);
%    xlabel('PLD (s)');
%    ylabel('\Delta M /M0');
%    xlim([0 5.5]);
%   
%     ylim([0 0.034*pp]);
%   
%    y1=get(gca,'ylim');
%    LL=line([1.5 1.5], y1);
%    set(LL,'LineWidth',2);
%    legend(txt1);
%    txt1 = {}; 
%    
%   hold off
% 
%     fname=fullfile(codepath,strcat('OVERLAY_SPIRAL_Slice_ATTchange_Perf_vs_PLD'));
%     print(strcat(fname),'-djpeg');  
%  
figure(5);

  
  sig = eval(sprintf('Oct12_sig_pump%d', pump(pp)));
 
   ft = fittype('DeltaM(x,lambda,tau,T1blood,p1,p2)', 'problem', {'lambda','tau','T1blood'},'coefficients',{'p1','p2'});  
  [ff, gof]= fit((asl.pld+asl.tau)',(sig)',ft,'StartPoint',[1.,0.05],'problem',{asl.lambda,asl.tau,asl.T1blood});
  % gof
   
   %%%% Use Fit parameters to plot Buxton Model. 
   
   DM =DeltaM(smooth_pld,asl.lambda,asl.tau,asl.T1blood,ff.p1,ff.p2);
   h=plot(asl.pld+asl.tau,sig,'+','Color',mycol(pp));  
   
   tt = sprintf('Data: PR = %d mL/min',pump(pp));
   txt = [txt; tt];
   
   hold on 
   hf=plot(smooth_pld, DM,'LineWidth',2,'LineStyle', '--', 'Color', mycol(pp)); % plot best fit curve, shifted x-axis values
%    tt = sprintf('BM: PR = %d mL/min',pump(pp));
%    txt = [txt; tt];
  hf.Annotation.LegendInformation.IconDisplayStyle = 'off';
   dt = sprintf('%4.2f',ff.p1+asl.tau)
  flow = sprintf('%4.2f',ff.p2*unit_factor) ;
    str = {strcat('dt=', num2str(dt)),strcat('flow=', num2str(flow))};
%   text(0.25, 0.01+ pp*0.01, str,'Color',mycol(pp),'FontSize',14);
%   
  set(gca,'FontSize',24);   % axis font size
  set(h,'LineWidth',2,'MarkerSize',13);


  
  hold on
  

end

%hold off
%   title('3 mm Phantom, BS = 0','FontSize',16);
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
%    fname=fullfile(codepath,strcat('par_maps/results/FitResults/2M0_Analysis/Thresholded/OVERLAY_Perf_vs_PLD_th',num2str(threshval)));
%     print(strcat(fname),'-djpeg');  
%  else  
%     fname=fullfile(codepath,strcat('par_maps/results/FitResults/2M0_Analysis/NoSmooth/OVERLAY_Perf_vs_PLD_R',num2str(R),'The',num2str(angle),'Slice',num2str(slice),'OzMax',num2str(Oz_max)));
%     print(strcat(fname),'-djpeg');  
    
%             fname=fullfile(codepath,'par_maps/results/FitResults/2M0_Analysis/NoSmooth/OVERLAY_JET_Perf_vs_PLD');
%     print(strcat(fname),'-djpeg');  
%      
    
%     
%         fname=fullfile(codepath,'par_maps/results/FitResults/2M0_Analysis/NoSmooth/OVERLAY_FULPahntom_Perf_vs_PLD');
%     print(strcat(fname),'-djpeg');  
%      
%     
%  end 



    fname=fullfile(codepath,strcat('OVERLAY_SPIRAL_Slice',num2str(slidx),'_Perf_vs_PLD'));
    print(strcat(fname),'-djpeg');  
    