%%%   Dicom 2 nifti converter. Convers one PLD series a time.
%%% 
% 
% clear;
% clc;
% 
% codepath = pwd;
% cd ../
% dicompath = fullfile(pwd,'Dicom','SUB01','20191006.Techdev_Reprod_Oct6_Set3.19.10.06_14_00_21_DST_1.3.12.2.1107.5.2.43.167024/');
% 
% P = spm_select('FPlist',dicompath,'.*dcm$');
% 
% cd(codepath);
% Niftipath = fullfile(pwd,'../Nifti','SUB01');
% 
% cd(codepath);
% 
% filenames = [];
% dicomnames = [];
% 
% preseriesnumber = 10000;
% for s = 1:size(P,1)
%    info = dicominfo(P(s,:)) ;
%    flag = findstr(info.SeriesDescription,'ASL');
%       
%    if isempty(flag)       
%        continue;
%    end
%    
%    if info.SeriesNumber == preseriesnumber
%         filenames = strvcat(filenames,P(s,:));
%    else
%        hdr = spm_dicom_headers(filenames);
%        spm_dicom_convert(hdr,'all','series','nii',Niftipath);
%        
%        filenames = P(s,:);
%        preseriesnumber = info.SeriesNumber;
%    end
%    
% end





clear
codepath = pwd;
Niftipath = fullfile(pwd,'../Nifti');

cd(codepath);

A = dir(Niftipath); 
A = {A.name};
A(ismember(A,{'.','..','.DS_Store'})) = [];


 LD=1500;
 pld_msec=[100,300,500,700,1000,1300,1700,2200,2700,3200,3500];
 pump = [200, 300, 400];
%  pump = 400 ; 



for sb = 1:length(A)    % Loop over the number of Subjects
  
    subid = A{sb}; 
    for ipp = 1:length(pump)            
       for i =1:length(pld_msec)
           
           aslsesdir = dir(fullfile(Niftipath,A{sb},strcat('3DSp_ubpCASL_TR6s_LD',num2str(LD),'_PLD',num2str(pld_msec(i)),'_BS0_PUMP',num2str(pump(ipp)),'_ASL*')));
           m0sesdir = dir(fullfile(Niftipath,A{sb},strcat('3DSp_ubpCASL_TR6s_NoLabel_BS0_PUMP',num2str(pump(ipp)),'_ASL*')));   % using Unlabeled ASL images as M0 for each pump rate.
                     
   
%%% Experiments
% 
%  		  aslsesdir = dir(fullfile(Niftipath,A{sb},strcat('3DSp_ubpCASL_TR6s_LD',num2str(LD),'_PLD',num2str(pld_msec(i)),'_FA20_PUMP',num2str(pump),'_ASL*'))); 
%         m0sesdir = dir(fullfile(Niftipath,A{sb},strcat('3DSp_ubpCASL_TR6s_LD',num2str(LD),'_PLD',num2str(pld_msec(i)),'_FA20_PUMP',num2str(pump),'_M0*')));
%        
% 
%  		  aslsesdir = dir(fullfile(Niftipath,A{sb},strcat('3DSp_ubpCASL_TR6s_LD',num2str(LD),'_PLD',num2str(pld_msec(i)),'_GmaxGavg4_PUMP',num2str(pump),'_ASL*'))); 
%         m0sesdir = dir(fullfile(Niftipath,A{sb},strcat('3DSp_ubpCASL_TR6s_LD',num2str(LD),'_PLD',num2str(pld_msec(i)),'_GmaxGavg4_PUMP',num2str(pump),'_M0*')));
%            
% 
%  		  aslsesdir = dir(fullfile(Niftipath,A{sb},strcat('3DSp_ubpCASL_TR6s_LD',num2str(LD),'_PLD',num2str(pld_msec(i)),'_ASYM_PUMP',num2str(pump),'_ASL*'))); 
%          m0sesdir = dir(fullfile(Niftipath,A{sb},strcat('3DSp_ubpCASL_TR6s_LD',num2str(LD),'_PLD',num2str(pld_msec(i)),'_ASYM_PUMP',num2str(pump),'_M0*')));
%              
         
         
         
         
        aslprefix = strcat('ASL_1shot_LD',num2str(LD),'_PLD',num2str(pld_msec(i)),'msec');
         M0prefix = strcat('M0_1shot_LD',num2str(LD),'_PLD',num2str(pld_msec(i)),'msec'); 
            
        origaslprefix = 's';
        origM0prefix = 's';
      

        aslpath = fullfile(aslsesdir.folder,aslsesdir.name)
         m0path = fullfile(m0sesdir.folder,m0sesdir.name)
        
    

       files = spm_select('FPlist',aslpath,['^' origaslprefix '.*nii$']);
      
        if ~isempty(files)
            spm_file_merge(files,[aslprefix '_' subid '.nii']);  % merges all repeatition (10 files in my case), and creates a new files with aslprefix names containing 4D volumes. 
        end 
       
        
       files = spm_select('FPlist',m0path,['^' origM0prefix '.*nii$']);
        
        if ~isempty(files)
            spm_file_merge(files,[M0prefix '_' subid '.nii']);        
        end    
        
     
          
        P = spm_select('ExtFPlist',m0path,['^' M0prefix '.*nii']);
        M0 = mean(spm_read_vols(spm_vol(P)),4);
        
        v = spm_vol(P);
        v0 = v(1);
        fields={'pinfo'};    % Removing this field to extract correct info from memory when writing info as 64-bit float.  
        v0=rmfield(v0,fields);
        v0.dt = [64 0];  %  DT_FLOAT   16   /* float (32 bits/voxel) */ from NIFTI header
         
        v0.fname = fullfile(aslpath,['meanM0_' subid '.nii']);
        spm_write_vol(v0,M0);
        
      end   %   length(pld_msec)

    end     %   length(pump) 
end         %   length(A)




