function [meanperf,M0,meanCBF,weight] = EPI_ReturnTriplet(aslsesdir,subid,pld_val,asl)
 
       aslpath = fullfile(aslsesdir.folder,aslsesdir.name);
       ROIpath = fullfile(aslsesdir.folder);
       
        %%% Shoul be the same names as in Convert_to_Nifti.m
        aslprefix = strcat('ASL_1shot_LD',num2str(asl.LD),'_PLD',num2str(pld_val),'msec');
         M0prefix = strcat('M0_1shot_LD',num2str(asl.LD),'_PLD',num2str(pld_val),'msec');  % 100 meas from unlabeled  L/C dataset merged to create 100 pseudo M0
         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Comment this section aafter the first run , no need to remearge the
%%%% data files, once it's already done
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%        origaslprefix = 'f';
%         origM0prefix = 'f';
% 
%       
%        files = spm_select('FPlist',aslpath,['^' origaslprefix '.*nii$']);
%       
%         if ~isempty(files)
%             spm_file_merge(files,[aslprefix '_' subid '.nii']);  % merges all repeatition (10 files in my case), and creates a new files with aslprefix names containing 4D volumes. 
%         end 
%        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        

%%%%%  Read ROI mask      
        P=spm_select('ExtFPlist', ROIpath,'^ConicalBaseROI_R16The60Slice20OzMax31.nii');  
        ROImask = spm_read_vols(spm_vol(P));
                
%%%%   Read mearged 4D volumes: ASL 
   
        P=spm_select('ExtFPlist',aslpath,['^' aslprefix '.*nii']);  
        v = spm_vol(P);
        Y = spm_read_vols(v);
%%%%%  Control - Label with step size = 2; (2:end) means go through entire volume   
        perf = Y(:,:,:,2:2:end) - Y(:,:,:,1:2:end); 
       % perf=vpa(perf,10);
        meanperf = mean(perf,4);   % mean along the 4th dimention volume index
        var_perf = var(perf,1,4);    % variance between 10 repetitions        
        std_perf = std(perf,1,4);    % Standard Deviation between 10 repetitions
          weight = 1./var_perf;

          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    These histograms are intended for Data quality and Calibration checks 
%%%%%%%    Histo: Mean Label, Control, M0 and dM of the time series avg-ed within the ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    figure(1)
%         edges = [0:15:1500];
         
%          subplot(4,1,1); 
%          hh = mean(Y(:,:,:,1:2:end), 4).*ROImask;
%          mskhh= hh~= 0;
%          
%          hLab = histogram(hh(mskhh), edges)
%          title ('Label');
%           
%           
%          subplot(4,1,2);    
%          hh = mean(Y(:,:,:,2:2:end), 4).*ROImask;
%          mskhh = hh~= 0;
%          
%     
%          hCon = histogram(hh(mskhh), edges) 
%          title ('Control');
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
 %%%% Read mearged 4D volumes: M0 
   
%         P = spm_select('ExtFPlist',m0path,['^' M0prefix '.*nii']);
%         M0 = mean(spm_read_vols(spm_vol(P)),4);
%            
%         var_M0 = var(spm_read_vols(spm_vol(P)),1,4);    % variance between 100 repetitions        
%         std_M0 = std(spm_read_vols(spm_vol(P)),1,4);    % Standard Deviation between 100 repetitions
% 
%          
         
     %%%  EPI  --- Control is taken as MO
        M0 = mean(Y(:,:,:,2:2:end),4);
 
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%   
%         figure(1)
%          subplot(4,1,3);            
%          hh = M0.*ROImask;
%          mskhh = hh~= 0;  
%          hM0 = histogram(hh(mskhh), 100) 
%          title ('M0');
%        
%                  
%          subplot(4,1,4);
%             
%          hh = meanperf.*ROImask;
%          mskhh = hh~= 0;  
%          hPerf = histogram(hh(mskhh), 100) 
%          title ('Mean Perf \Delta M');
%          xlabel ('Voxel Intensity inside Connical ROI');
%             
%          
%         fname=fullfile(aslpath,'hist_LabConPerf');
%         print(strcat(fname),'-djpeg');  
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%    These histograms are intended to Data quality checks 
%%%%%%%%%   Histo: dM/M0 avg-ed within the ROI for each measurement 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  figure(1)       
%   edges = [-0.05:0.001:0.05];
%   hh = perf .* ROImask;
%  % size(hh);       
% 
% for i = 1:50
% 
%   m = squeeze(hh(:,:,:,i))./M0;
%   
%   m(isnan(m))=0;
%   m(isinf(m))=0;
%   
%   msk = m~=0;
%   mm(i) = mean(m(msk));
%     
% end
% 
% size(mm);
% histogram(mm(:),edges);
% 
% 
%   title ('Mean Perf \Delta M / M0');
%   xlabel ('Mean \Delta M / M0 avg-ed in ROI for time series');
%                      
%   fname=fullfile(aslpath,'hist_Perf');
%   print(strcat(fname),'-djpeg');  
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        
%%%%%  Scaling Magnitude images by a factor of 10, since ASL signal is
%%%%%  tiny, it will be lost in the nose otherwise. Creating mask to avoid
%%%%%  dividing by 0.
        M0 = M0 / 10;
        msk = M0 > (max(M0(:))/10); 
       
          
         
        %%%%%%  CBF factor from Alsop et al. Recommended implementation for
        %%%%%%  PCASL  %MRM White Paper 2014
%%%%%%  CBF_factor=6000*lambda*exp(PLD/T1blood)/(2*alpha*T1blood*(1-exp(-tau/T1blood))); 
%%%%%% For Phantom: lambda = 1, T1=1.5sec (tap water instead of blood) 


pld = pld_val/1000.; % convert to sec
lambda = asl.lambda;
T1b = asl.T1blood;
alpha = asl.alpha;
tau = asl.tau ;    
M0scale = asl.M0scale; 

        CBF_factor = 6000.0 .* lambda .* exp(pld ./ T1b) ./ ( 2 .* alpha .* M0scale .* T1b .* (1 - exp(-tau ./ T1b))); 
    
        
%%%%%  Only divide to non-zero M0 masked region, to avoid dividing by 0. 
        meanCBF = zeros(size(meanperf));   
        meanCBF(msk) = CBF_factor .* (meanperf(msk)./M0(msk));   % mL/100g/min
       % meanCBF(msk) = (meanperf(msk)./M0(msk));  
      
        
        
        v0 = v(1);
        fields={'pinfo'};    % Removing this field to extract correct info from memory when writing info as 64-bit float.  
        v0=rmfield(v0,fields);
        v0.dt = [64 0];  %  DT_FLOAT   16   /* float (32 bits/voxel) */ from NIFTI header
         
        v0.fname = fullfile(aslpath,['meanPerf_' subid '.nii']);   
        spm_write_vol(v0,meanperf);
        
%         v0.fname = fullfile(aslpath,['meanCBF_' subid '.nii']);
%          % v0.dt = [16 0];  %  DT_FLOAT   16   /* float (32 bits/voxel) */ from NIFTI header
%         spm_write_vol(v0,meanCBF);
%        
        
        v0.fname = fullfile(aslpath,['Variance_Perf_' subid '.nii']);
        spm_write_vol(v0,var_perf);
        
        v0.fname = fullfile(aslpath,['Weights_' subid '.nii']);   
        spm_write_vol(v0,weight);
        
        v0.fname = fullfile(aslpath,['STD_Perf_' subid '.nii']);
        spm_write_vol(v0,std_perf);
        
        
        M0 = 10 .* M0 ;  % Returns the origina M0, without downscaling
        
        v0.fname = fullfile(aslpath,['meanM0_' subid '.nii']);
        spm_write_vol(v0,M0);

%          
%         v0.fname = fullfile(aslpath,['Variance_M0_' subid '.nii']);
%         spm_write_vol(v0,var_M0);
%         
%         v0.fname = fullfile(aslpath,['STD_M0_' subid '.nii']);
%         spm_write_vol(v0,std_M0);

        
end



