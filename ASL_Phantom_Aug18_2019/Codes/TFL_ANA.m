%%%%   Fitting T1 
%%%%


clear; clc;
close all;

codepath = pwd;

cd ../Nifti
Niftipath = pwd;

cd(codepath);

A = dir(Niftipath); 
A = {A.name};
A(ismember(A,{'.','..','.DS_Store'})) = [];


mycol=['g','b','r','k','m'];
TI= [500 1000 1250 1500 2000];
sig = zeros(length(TI),1);

direction= [0 0 1];
time=0:0.01:2.5;


   P=spm_select('ExtFPlist', fullfile(Niftipath,A{1}),'^TFL_PHANTOM_MASK.nii');  % Pump rate ana 
   GeomMask = spm_read_vols(spm_vol(P));
   
   GeomMask = imrotate(GeomMask,90);
   

% figure(1)
%    imshow(GeomMask(:,:,1),[], 'InitialMagnification', 1000)
%   
   
figure(2)

for i = 1:length(TI)

tfl = dir(fullfile(Niftipath,A{1},strcat('tfl_TI',num2str(TI(i)),'ms_*')));
 impath = fullfile(tfl.folder,tfl.name);

     P=spm_select('ExtFPlist',impath,'s.*nii');
     v = spm_vol(P);
     Y = spm_read_vols(v);
      

     Y = imrotate(Y,90);   
  
    % imshow(Y(:,:,1),[], 'InitialMagnification', 1000);
      
     IM = Y .* GeomMask;
     
   % imshow(IM(:,:,1),[], 'InitialMagnification', 1000);
     
     msk = IM ~= 0;
     sig(i) = mean(IM(msk));
     
 

 
end

[~, indx] = min(sig);
for i = 1:(indx)
 sig(i) = -sig(i);

end 
 
% h = plot( TI./1000 ,sig,'+','Color',mycol(i)); 
%      
% set(h,'LineWidth',3,'MarkerSize',12);
% 
% xlim([0 2.5]);
% %ylim([-1 1.1]);
% 
% hold on 
t = (TI./1000);

ft = fittype( @(p1, p2, x) p2*(1 - 2*exp(-x/p1))) ; 
[ff, gof]= fit(t',(sig),ft,'StartPoint', [1.5, 350])
%  gof
%  ff.p1 


  h = plot(ff, t', sig);
 set(h,'LineWidth',3,'MarkerSize',32);
  

  xlim([0 2.5]);
  xlabel('Time (s)');
  ylabel('Intensity');
  set(gca,'FontSize',24);   % axis font size
 
 
xl = xlim;
yl = ylim;
xt = 0.05 * (xl(2)-xl(1)) + xl(1);
yt = 0.90 * (yl(2)-yl(1)) + yl(1);
cap1 = sprintf('%f *(1 - 2*exp(-x/%f))', ff.p2, ff.p1 );
cap2 = sprintf( 'R2 = %4.2f', gof.rsquare);

text(xt, yt, cap1, 'FontSize', 20, 'Color', 'b', 'FontWeight', 'bold'); 
text(xt, yt -20, cap2, 'FontSize', 20, 'Color', 'b', 'FontWeight', 'bold') 
 
 
 
 
hold off


