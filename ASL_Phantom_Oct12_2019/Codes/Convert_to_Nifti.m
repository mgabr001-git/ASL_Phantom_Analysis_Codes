%%%%%   This Macro Converts dicom files to nifti

clear;
close all;

codepath = pwd;
cd ../Dicom
Dicompath = pwd;
cd ../Nifti
Niftipath = pwd;

cd(codepath);


A = dir(Dicompath);
A = {A.name};
A(ismember(A,{'.','..','.DS_Store'})) = [];


 
for sb = 1:length(A)     % Loop over the number of Subjects
  %  dicomdir = my_genpath(fullfile(Dicompath,A{sb}));
    dicomdir = my_genpath(fullfile(Dicompath,A{sb},'TFL_TI'));
    for s = 1:size(dicomdir,1)
        P = spm_select('FPlist',dicomdir(s,:),'.*.dcm');
        if isempty(P)
             continue;
        end
        hdr = spm_dicom_headers(P);
        niftidir = fullfile(Niftipath,A{sb});
        if ~exist(niftidir,'dir') 
            mkdir(niftidir);
        end
        cd(niftidir);
        spm_dicom_convert(hdr,'all','series','nii');
        cd(codepath);
    end
end

% 
% 
%  LD=1800;
%  pld_msec=[1800];  
%  
% for sb = 1:length(A)    % Loop over the number of Subjects
%   
%     subid = A{sb};
% %       PAR.subject(sb).anatdir = ......
% 
%     for i =1:length(pld_msec)
%         
%  
%         aslsesdir = dir(fullfile(Niftipath,A{sb},strcat('3DSp_ubpCASL_BS90_LD',num2str(LD),'_PLD',num2str(pld_msec(i)),'_ASL*'))); 
%          M0sesdir = dir(fullfile(Niftipath,A{sb},strcat('3DSp_ubpCASL_BS90_LD',num2str(LD),'_PLD',num2str(pld_msec(i)),'_M0*')));
%           
%         aslprefix = strcat('ASL_1shot_LD',num2str(LD),'_PLD',num2str(pld_msec(i)),'msec');
%          M0prefix = strcat('M0_1shot_LD',num2str(LD),'_PLD',num2str(pld_msec(i)),'msec'); 
%                 
%         origaslprefix = 's';
%         origM0prefix = 's';
%       
% 
%         aslpath = fullfile(aslsesdir.folder,aslsesdir.name);
%          m0path = fullfile(M0sesdir.folder,M0sesdir.name);
%         
%     
% 
%        files = spm_select('FPlist',aslpath,['^' origaslprefix '.*nii$']);
%       
%         if ~isempty(files)
%             spm_file_merge(files,[aslprefix '_' subid '.nii']);  % merges all repeatition (10 files in my case), and creates a new files with aslprefix names containing 4D volumes. 
%         end 
%        
%         
%        files = spm_select('FPlist',m0path,['^' origM0prefix '.*nii$']);
%         
%         if ~isempty(files)
%             spm_file_merge(files,[M0prefix '_' subid '.nii']);        
%         end    
%         
% 
% 
%     end
% end
% 

