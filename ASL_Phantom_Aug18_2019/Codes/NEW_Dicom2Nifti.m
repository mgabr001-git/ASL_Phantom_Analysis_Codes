% new dicom 2 nifti

clear
codepath = pwd;
cd ../
dicompath = fullfile(pwd,'Dicom','SUB01','20190818.Techdev_New3mm_Reproduce_set1.19.08.18_12_08_58_DST_1.3.12.2.1107.5.2.43.167024');

P = spm_select('FPlist',dicompath,'.*dcm$');

cd(codepath);
niftipath = fullfile(pwd,'../Nifti','SUB01');

cd(codepath);

filenames = [];
dicomnames = [];

preseriesnumber = 10000;
for s = 1:size(P,1)
   info = dicominfo(P(s,:)) ;
   flag = findstr(info.SeriesDescription,'ASL');
   
   
   if isempty(flag)       
       continue;
   end
   
   if info.SeriesNumber == preseriesnumber
        filenames = strvcat(filenames,P(s,:));
   else
       hdr = spm_dicom_headers(filenames);
       spm_dicom_convert(hdr,'all','series','nii',niftipath);
       
       filenames = P(s,:);
       preseriesnumber = info.SeriesNumber;
   end
   
end




