%%%%  
close all
codepath = pwd;

% % % % % 
% % % % % 
% % % % % Phantom 3mm.  
% % % % % sigma 3
% % % % % contrast -0.02 - +0.02
% % % % % LABELED, BS0
% % % % % 
% % % % % 
% % % % % PUMP 400. 
% % % % % threshold  = 0.03 (un-smoothed), 0.015 (smoothed)
% % % % % volume 5  - PLD = 1000 ms
% % % % %  SAME ROI is segmentation is applied to PUMP 200, 300 mL/min
% % % % 
% % % % asl.pld_msec=[100, 300, 500, 700, 1000, 1300, 1700, 2200, 2700, 3200, 3500];
% % % %            vol  1.   2.   3.   4.   5.    6.    7.    8.    9.   10.   11. 



%% Oct 12 Pump 400, 2019 Data 

%%%%%%%%%%%%%   2 M0   %%%%%%%%%%%%%%%%%%%%

 
%%%   Sphere Segment ROI   %%%%%%%%%%%%
%%%  R = 16,   s > 23, Zmax =33

%%%% NO SMOOTHING 


sig3 = [ 0.0222    0.0264    0.0316    0.0333    0.0322    0.0287    0.0237    0.0170    0.0139    0.0106    0.0105];





%% EXPERIMENTS   Aug 31 , 2019 Data 
%%%%% PR 400 mL/min , BS = 0 


%%%%%%   SET 1  %%%%%%

%%%  Spehrical Segment ROI   %%%%%%%%%%%%

%%%  R = 16, s > 23, Zmax =33 

 %%%% NO SMOOTHING 

%%% ASYMMETRIC  Labeling , 50 meas
sig5_set1 = [0.0139    0.0146    0.0158    0.0162    0.0146    0.0122    0.0114    0.0087    0.0069    0.0070    0.0059];

%%%  GmaxGavg4 , 50 meas
sig6_set1 = [0.0168    0.0203    0.0230    0.0248    0.0231    0.0222    0.0173    0.0144    0.0095    0.0091    0.0081];
 
%%% FA20 , 50 meas
sig7_set1 = [0.0152    0.0193    0.0215    0.0239    0.0233    0.0214    0.0172    0.0132    0.0107    0.0083    0.0070];







%%%%%%   SET 2  %%%%%%


%%%   Connical base + Spehrical top ROI   %%%%%%%%%%%%

%%%  R = 16, s > 23, Zmax =33 

%%%% NO SMOOTHING 

% % %%% ASYMMETRIC  Labeling , 50 meas
sig5_set2 = [0.0124    0.0151    0.0157    0.0154    0.0142    0.0135    0.0117    0.0089    0.0080    0.0075    0.0056];

% %%%  GmaxGavg4 , 50 meas
sig6_set2 = [0.0174    0.0210    0.0231    0.0253    0.0244    0.0217    0.0187    0.0145    0.0114    0.0081    0.0090];
 
%%% FA20 , 50 meas
sig7_set2 = [0.0180    0.0212    0.0254    0.0269    0.0250    0.0219    0.0200    0.0151    0.0105    0.0085    0.0080];







asl.pld=[0.1,0.3,0.5,0.7,1.0,1.3,1.7,2.2,2.7,3.2,3.5];
asl.LD = 1500;        % Label Duration LD in msec
asl.TR =6;            % sec
asl.lambda = 1;
asl.T1blood =2.1;
asl.alpha=0.85;
asl.tau = asl.LD/1000;        %  tau = Label Duration LD in sec 


unit_factor = 6000;   %%%  1/Vvoxel = 1/(53 * 10-3) =  1000 /53  (Flow unit conversion from mL/Vvol/sec into mL/mL/s  or fractional exchange per sec)
smooth_pld=0.1:0.02:5.2;
mycol=['r','b','g','k','m'];


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% dM vs PLD for each Pump Rate: OVERLAY
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
% figure;
% set(gcf, 'Position',  [50, 900, 900, 550]);
% 
%  
%    ft = fittype('DeltaM(x,lambda,tau,T1blood,p1,p2)', 'problem', {'lambda','tau','T1blood'},'coefficients',{'p1','p2'}); 
%    [ff, gof]= fit((asl.pld)',(sig1)',ft,'StartPoint',[0.25,0.05],'problem',{asl.lambda,asl.tau,asl.T1blood});
%    ph =  'Ph 3mm, 200 mL/min'  
% %    dt = ff.p1
% %    flow = ff.p2
%   dt = sprintf('%4.3f',ff.p1);
%   flow = sprintf('%4.3f',ff.p2);
%    
%    %%%% Use Variance Wweighted fit parameters to plot Buxton Model. 
%    DM1 =DeltaM(smooth_pld,asl.lambda,asl.tau,asl.T1blood,ff.p1,ff.p2);
%   % h=errorbar(asl.pld+asl.tau,sig1,err1,err1,'g+');
%    h=plot(asl.pld+asl.tau,sig1,'g+');  
%    hold on 
%    plot(smooth_pld+asl.tau, DM1,'LineWidth',2,'LineStyle', '--', 'Color', 'green'); % plot best fit curve, shifted x-axis values
%   
%   set(gca,'FontSize',14);   % axis font size
%   set(h,'LineWidth',2,'MarkerSize',12);
%   str = {strcat('dt=', num2str(dt)),strcat('flow=', num2str(flow))};
%   text(0.5, 0.02, str, 'Color','green','FontSize',14);
%    
%   
%   
%   
%   hold on 
%   
%  
%    ft = fittype('DeltaM(x,lambda,tau,T1blood,p1,p2)', 'problem', {'lambda','tau','T1blood'},'coefficients',{'p1','p2'}); 
%    [ff, gof]= fit((asl.pld)',(sig2)',ft,'StartPoint',[0.25,0.05],'problem',{asl.lambda,asl.tau,asl.T1blood});
%    
%   ph =  'Ph 3mm, 300 mL/min' 
% %   dt = ff.p1
% %   flow = ff.p2
%   dt = sprintf('%4.3f',ff.p1);
%   flow = sprintf('%4.3f',ff.p2);
%    
% 
%   DM2 =DeltaM(smooth_pld,asl.lambda,asl.tau,asl.T1blood,ff.p1,ff.p2);
%   h=plot(asl.pld+asl.tau,sig2,'b+');
%   hold on 
%   plot(smooth_pld+asl.tau, DM2,'LineWidth',2,'LineStyle', '--', 'Color', 'blue'); % plot best fit curve, shifted x-axis values
% 
%   set(gca,'FontSize',14);   % axis font size
%   set(h,'LineWidth',2,'MarkerSize',12);
%   str = {strcat('dt=', num2str(dt)),strcat('flow=', num2str(flow))};
%   text(0.5, 0.03, str, 'Color','blue','FontSize',14);
%    
%   
%   
%   hold on 
%   
%   
%    ft = fittype('DeltaM(x,lambda,tau,T1blood,p1,p2)', 'problem', {'lambda','tau','T1blood'},'coefficients',{'p1','p2'}); 
%    [ff, gof]= fit((asl.pld)',(sig3)',ft,'StartPoint',[0.25,0.05],'problem',{asl.lambda,asl.tau,asl.T1blood});
%    
%   ph =  'Ph 3mm, 400 mL/min' 
% %   dt = ff.p1
% %   flow = ff.p2
%   dt = sprintf('%4.3f',ff.p1);
%   flow = sprintf('%4.3f',ff.p2);
%    
%    str = {strcat('dt=', num2str(dt)),strcat('flow=', num2str(flow))};
%   
%   
% %  h=plot(ff,asl.pld,sig2,'go'); 
%   
%   DM3 =DeltaM(smooth_pld,asl.lambda,asl.tau,asl.T1blood,ff.p1,ff.p2);
%   h=plot(asl.pld+asl.tau,sig3,'r+');  % plot data points
%   hold on 
%   plot(smooth_pld+asl.tau, DM3,'LineWidth',2,'LineStyle', '--', 'Color', 'red'); % plot best fit curve, shifted x-axis values
% 
%  
%   set(gca,'FontSize',14);   % axis font size
%   set(h,'LineWidth',2,'MarkerSize',12);
%   str = {strcat('dt=', num2str(dt)),strcat('flow=', num2str(flow))};
%   text(0.5, 0.04, str, 'Color','red','FontSize',14); 
%   
% %   hold on 
% %   
% % 
% %    ft = fittype('DeltaM_VarWgtFit(x,lambda,tau,T1blood,wgt,p1,p2)', 'problem', {'lambda','tau','T1blood','wgt'},'coefficients',{'p1','p2'}); 
% %    [ff, gof]= fit((asl.pld)',(wgt4.*sig4)',ft,'StartPoint',[-0.25,0.5],'problem',{asl.lambda,asl.tau,asl.T1blood,wgt4});
% %    
% %   ph =  'Ph 3mm, 600 mL/min' 
% %   dt = ff.p1
% %   flow = ff.p2
% % 
% % %  h=plot(ff,asl.pld,sig2,'go'); 
% %   
% %   DM4 =DeltaM(smooth_pld,asl.lambda,asl.tau,asl.T1blood,ff.p1,ff.p2);
% %   h=errorbar(asl.pld+asl.tau,sig4,err4,err4,'k+');  % plot data points
% %   hold on 
% %   plot(smooth_pld+asl.tau, DM4,'LineWidth',2,'LineStyle', '--', 'Color', 'black'); % plot best fit curve, shifted x-axis values
%  
%   set(gca,'FontSize',14);   % axis font size
%   set(h,'LineWidth',2,'MarkerSize',12);
%   
%    
%   
%    set(gca,'FontSize',14);   % axis font size
%    hTitle = title('3 mm Phantom, BS = 0');
%   
%    set(hTitle,'FontSize',16);  
%    set(h,'LineWidth',2,'MarkerSize',12);
%    xlabel('PLD (s)');
%    ylabel('\Delta M /M0');
%    xlim([0 5.5]);
%    ylim([0 0.1]);
%  
%    y1=get(gca,'ylim');
%    LL=line([1.5 1.5], y1);
%    set(LL,'LineWidth',2);
%       
%    
%  legend('3 mm, PR = 200 mL/min','BM: PR 200 mL/min',...
%         '3 mm, PR = 300 mL/min','BM: PR 300 mL/min',...
%         '3 mm, PR = 400 mL/min','BM: PR 400 mL/min');
%  
% % 
% %  legend('3 mm, PR = 300 mL/min','BM: PR 300 mL/min',...
% %         '3 mm, PR = 400 mL/min','BM: PR 400 mL/min');
% %  
%     hold off
%   
%    fname=fullfile(codepath,'par_maps/results/FitResults/SingleROI_NoSmoothing_pump400_th003_March30_100M0/OVERLAY_Perf_vs_PLD_NoSmoothNoWeight_3pumps');
%     print(strcat(fname),'-djpeg');   
%    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %% Plot dM vs Pump Rates
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
%  pld_msec=[100,300,500,700,1000,1300,1700,2200,2700,3200,3500];
%  pump=[200, 300, 400];
%  
% %  
% %  figure; 
% % 
% %  for i=1:11
% % 
% % pump_pld = [sig1(i), sig2(i), sig3(i)];
% % 
% %   h=plot(pump,pump_pld,'b+');
% %   
% %    xlabel('PUMP rate (mL/min)'); 
% %    ylabel('\Delta M /M0');         
% %    set(gca,'FontSize',14);   % axis font size  
% %    set(h,'LineWidth',2,'MarkerSize',12);
% %    xlim([150 450]);  
% %    hTitle = title(strcat('BS = 0, PLD = ',num2str(pld_msec(i)),' ms')); 
% %   
% %  fname=fullfile(codepath,['par_maps/results/FitResults/IndividualSegment_thresh004/dM_vs_PumpRate_PLD',num2str(pld_msec(i))]);
% %  print(strcat(fname),'-djpeg');  
% % 
% %  end
%   
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Plot dM vs PLD for each pump rate, Individual plots      
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
%  figure;
%  set(gcf, 'Position',  [50, 900, 900, 550]);
%    
% 
%  for i=1:3
% 
% sig=eval(strcat('sig',num2str(i)));
% erry=eval(strcat('err',num2str(i)));
% wgt = eval(strcat('wgt',num2str(i)));
% 
% 
% 
%    ft = fittype('DeltaM(x,lambda,tau,T1blood,p1,p2)', 'problem', {'lambda','tau','T1blood'},'coefficients',{'p1','p2'}); 
%    [ff, gof]= fit(asl.pld',(sig)',ft,'StartPoint',[0.25,0.05],'problem',{asl.lambda,asl.tau,asl.T1blood});   
% 
%    ph =  strcat('Ph 3mm, Pump=',num2str(pump(i)),' mL/min') ;
% %    dt = ff.p1;
% %    flow = ff.p2;
%    
%    %%%% Use Variance Wweighted fit parameters to plot Buxton Model. 
%    DM =DeltaM(smooth_pld,asl.lambda,asl.tau,asl.T1blood,ff.p1,ff.p2);
%  
%    h=plot(asl.pld+asl.tau,sig,'b+');
%     
%    hold on
%    plot(smooth_pld+asl.tau, DM,'LineWidth',2, 'Color', 'red'); % plot best fit curve, shifted x-axis values
%     
%    hold off 
%    
%    
%    set(gca,'FontSize',14);   % axis font size
%    hTitle = title('3 mm Phantom, BS = 0');
%   
%    set(hTitle,'FontSize',16);  
%    set(h,'LineWidth',2,'MarkerSize',12);
%    xlabel('PLD (s)');
%    ylabel('\Delta M /M0');
%    xlim([0 5.5] );
%    ylim([0 0.1] );
%    y1=get(gca,'ylim');
%    LL=line([1.5 1.5], y1);
%    set(LL,'LineWidth',2);
% 
%    legend(strcat('Ph 3mm, Pump=',num2str(pump(i)),' mL/min'),strcat('BM, Pump=',num2str(pump(i)),' mL/min'));
%      
%    dt = sprintf('%4.3f',ff.p1);
%    flow = sprintf('%4.3f',ff.p2)
%    
%    str = {strcat('dt=', num2str(dt)),strcat('flow=', num2str(flow))};
%    text(0.5, 0.04,str, 'Color','red','FontSize',14);
%    
%    
%    
%    fname=fullfile(codepath,['par_maps/results/FitResults/SingleROI_NoSmoothing_pump400_th003_March30_100M0/TotalPerf_PR',num2str(pump(i))]);
%     print(strcat(fname),'-djpeg');  
%      
%  end
%   
%  
%  
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%  EXPERIMENTS
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 figure;
 set(gcf, 'Position',  [50, 900, 900, 550]);
   

 for i=5:7

sig_set1=eval(strcat('sig',num2str(i),'_set1'));
% erry_set1=eval(strcat('err',num2str(i),'_set1'));
% wgt_set1 = eval(strcat('wgt',num2str(i),'_set1'));


sig_set2=eval(strcat('sig',num2str(i),'_set2'));
% erry_set2=eval(strcat('err',num2str(i),'_set2'));
% wgt_set2= eval(strcat('wgt',num2str(i),'_set2'));


   ft = fittype('DeltaM(x,lambda,tau,T1blood,p1,p2)', 'problem', {'lambda','tau','T1blood'},'coefficients',{'p1','p2'}); 
   [ff, gof]= fit((asl.pld+asl.tau)',(sig_set1)',ft,'StartPoint',[1.0,0.02],'problem',{asl.lambda,asl.tau,asl.T1blood});

     
   %%%% Use Variance Wweighted fit parameters to plot Buxton Model. 
   DM =DeltaM(smooth_pld,asl.lambda,asl.tau,asl.T1blood,ff.p1,ff.p2);
   h=plot(asl.pld+asl.tau,sig_set1,'+','Color','g');  
   set(h,'LineWidth',2,'MarkerSize',12);   
   
   hold on
   hf = plot(smooth_pld, DM,'LineWidth',2, 'Color', 'g'); % plot best fit curve, shifted x-axis values
   hf.Annotation.LegendInformation.IconDisplayStyle = 'off'; 

   dt = sprintf('%4.2f',ff.p1+asl.tau);
  flow = sprintf('%4.2f',ff.p2*unit_factor); 
%   str = {strcat('dt=', num2str(dt)),strcat('flow=', num2str(flow))};
%   text(0.25, 0.015, str, 'Color','red','FontSize',14);  
 
 
  hold on
  
  
   ft = fittype('DeltaM(x,lambda,tau,T1blood,p1,p2)', 'problem', {'lambda','tau','T1blood'},'coefficients',{'p1','p2'}); 
   [ff, gof]= fit((asl.pld+asl.tau)',(sig_set2)',ft,'StartPoint',[1.0,0.02],'problem',{asl.lambda,asl.tau,asl.T1blood});

     
   %%%% Use Variance Wweighted fit parameters to plot Buxton Model. 
   DM =DeltaM(smooth_pld,asl.lambda,asl.tau,asl.T1blood,ff.p1,ff.p2);
   h=plot(asl.pld+asl.tau,sig_set2,'+','Color','r');  
   set(h,'LineWidth',2,'MarkerSize',12);   
   
   hold on
   hf = plot(smooth_pld, DM,'LineWidth',2, 'Color', 'r'); % plot best fit curve, shifted x-axis values
    hf.Annotation.LegendInformation.IconDisplayStyle = 'off';

    dt = sprintf('%4.2f',ff.p1+asl.tau);
  flow = sprintf('%4.2f',ff.p2*unit_factor); 
%   str = {strcat('dt=', num2str(dt)),strcat('flow=', num2str(flow))};
%   text(0.25, 0.025, str, 'Color','green','FontSize',14);  
  
  
  

   %%%%%%%  Overlay NOMINAL    
   %%%%%%%  BS0, Pump 400 ml/Min, Symmetric, Gmax/Gavg = 8, FA =28.7

   hold on 
    ft = fittype('DeltaM(x,lambda,tau,T1blood,p1,p2)', 'problem', {'lambda','tau','T1blood'},'coefficients',{'p1','p2'});   
    [ff, gof]= fit((asl.pld+asl.tau)',(sig3)',ft,'StartPoint',[1.0,0.02],'problem',{asl.lambda,asl.tau,asl.T1blood});

   
  
   DM =DeltaM(smooth_pld,asl.lambda,asl.tau,asl.T1blood,ff.p1,ff.p2);
   h=plot(asl.pld+asl.tau,sig3,'b+');
   
   hold on 
   hf = plot(smooth_pld, DM,'LineWidth',2,'LineStyle', '--', 'Color', 'blue'); % plot best fit curve, shifted x-axis values
   hf.Annotation.LegendInformation.IconDisplayStyle = 'off';
  
   dt = sprintf('%4.2f',ff.p1+asl.tau);
   flow = sprintf('%4.2f',ff.p2*unit_factor); 
  
%   str = {strcat('dt=', num2str(dt)),strcat('flow=', num2str(flow))};
%   text(0.25, 0.035, str, 'Color','blue','FontSize',14);
  ylim([0 0.05] ); 
    
  
  
  
   xlabel('Time (s)');
   ylabel('\Delta M /M0');
   xlim([0 5.5]);  
   %ylim([0 0.06] );   
   y1=get(gca,'ylim');
   LL=line([1.5 1.5], y1);
   set(LL,'LineWidth',2);
% LL.Annotation.LegendInformation.IconDisplayStyle = 'off';
   
   
switch i
    case 5
       %  hTitle = title('ASYMMETRIC vs SYMMETRIC Labeling, PR = 400 mL/min');
         legend('Asymmetric, SET 1','Asymmetric, SET 2','REFERENCE: Symmetric');
           
    case 6
      %   hTitle = title('Gmax/Gavg = 4 vs Gmax/Gavg = 8 ,  PR = 400 mL/min');
         legend('GmaxGavg = 4, SET 1','GmaxGavg = 4, SET 2','REFERENCE: GmaxGavg = 8');
              
           
    case 7
     %    hTitle = title('FA = 20 vs FA = 28.7,  PR = 400 mL/min');
         legend('FA = 20, SET 1','FA = 20, SET 2','REFERENCE, FA = 28.7');

end
   
   
   set(gca,'FontSize',20);   % axis font size 
 %  set(hTitle,'FontSize',20);  
   set(h,'LineWidth',2,'MarkerSize',12);
 
   hold off   
   
%     
%     
% switch i
%     case 5
%     fname=fullfile(codepath,['par_maps/results/FitResults/Experiments/ConicalBase_ROI_overlayed/Smoothed/TotalPerf_PR400_ASYMMETRIC']);
%     print(strcat(fname),'-djpeg');  
%      
%     case 6
%     fname=fullfile(codepath,['par_maps/results/FitResults/Experiments/ConicalBase_ROI_overlayed/Smoothed/TotalPerf_PR400_GmaxGavg4']);
%     print(strcat(fname),'-djpeg');  
%             
%    
%     case 7
%     fname=fullfile(codepath,['par_maps/results/FitResults/Experiments/ConicalBase_ROI_overlayed/Smoothed/TotalPerf_PR400_FA20']);
%     print(strcat(fname),'-djpeg'); 
%     
% end
  
switch i
    case 5
    fname=fullfile(codepath,['par_maps/results/FitResults/Experiments/SphericalSegmentROI/Perf_PR400_ASYMMETRIC']);
    print(strcat(fname),'-djpeg');  
     
    case 6
    fname=fullfile(codepath,['par_maps/results/FitResults/Experiments/SphericalSegmentROI/Perf_PR400_GmaxGavg4']);
    print(strcat(fname),'-djpeg');  
            
   
    case 7
    fname=fullfile(codepath,['par_maps/results/FitResults/Experiments/SphericalSegmentROI/Perf_PR400_FA20']);
    print(strcat(fname),'-djpeg'); 
    
end
    


    
 end
  
 
 