close all;

R2cut = 0.86;
nbins = 20;



figure(1)  
set(gcf, 'Position',  [50, 450, 750, 500]);

att = niftiread('../NoSmooth_2M0_PUMP200_BS0/ATT_map_GeomMask_NoSmooth_PUMP200.nii');
R2  = niftiread('../NoSmooth_2M0_PUMP200_BS0/GoodnessOfFit_R2_GeomMask_NoSmooth_PUMP200.nii');

   

msk2 =  R2(:,:,:) > R2cut ;

Y = att .* msk2;
        
msk = Y(:)>0;
        
size(Y(msk))
histogram(Y(msk),nbins)
xlim([0 5.]);  
xlabel("ATT (s)");
set(gca,'FontSize',24);   % axis font size       
text(0.25, 80, 'Pump Rate','Color','k','FontSize',24);
text(0.25, 70, '200 mL/min','Color','k','FontSize',24);

 print('ATT_Hist_Pump200','-djpeg');  
  



figure(2)  
set(gcf, 'Position',  [850, 450, 750, 500]);

att = niftiread('../NoSmooth_2M0_PUMP300_BS0/ATT_map_GeomMask_NoSmooth_PUMP300.nii');
R2  = niftiread('../NoSmooth_2M0_PUMP300_BS0/GoodnessOfFit_R2_GeomMask_NoSmooth_PUMP300.nii');

msk2 =  R2(:,:,:) > R2cut ;

Y = att .* msk2;
        
msk = Y(:)>0;
        
size(Y(msk))
histogram(Y(msk),nbins)
xlim([0 5.]);        
xlabel("ATT (s)");
set(gca,'FontSize',24);   % axis font size       
text(0.25, 50, 'Pump Rate','Color','k','FontSize',24);
text(0.25, 44, '300 mL/min','Color','k','FontSize',24);

print('ATT_Hist_Pump300','-djpeg');  



figure(3) 
set(gcf, 'Position',  [1650, 450, 750, 500]);

att = niftiread('../NoSmooth_2M0_PUMP400_BS0/ATT_map_GeomMask_NoSmooth_PUMP400.nii');
R2  = niftiread('../NoSmooth_2M0_PUMP400_BS0/GoodnessOfFit_R2_GeomMask_NoSmooth_PUMP400.nii');


msk2 =  R2(:,:,:) > R2cut;

Y = att .* msk2;
        
msk = Y(:)>0;
        
size(Y(msk))
histogram(Y(msk),nbins)
xlim([0 5.]);
xlabel("ATT (s)");
set(gca,'FontSize',24);   % axis font size       
text(0.25, 25, 'Pump Rate','Color','k','FontSize',24);
text(0.25, 22, '400 mL/min','Color','k','FontSize',24);

print('ATT_Hist_Pump400','-djpeg');  


        