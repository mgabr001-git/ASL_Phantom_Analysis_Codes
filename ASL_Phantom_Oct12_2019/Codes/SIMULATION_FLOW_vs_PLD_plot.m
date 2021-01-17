 close all
 clear
 clc
 
%%% Column 1 = time, Column 2 = Full phantom + Inflow and Outflow tubes,  Column 3 = Full phantom, Column 6 = Spherical Segment ROI

T = importdata('SIMULATIONS_T1_2.1s/markerAvrg_200.dat');
t =  T.data(:,1) ;
sim_full_InOut_pump200 =T.data(:,2) ;  % full phantom including Inflow /Outflow tubes
sim_full_pump200 =T.data(:,3) ;        % full phantom without Inflow /Outflow tubes
sim_roi_pump200 =T.data(:,6) ;


T = importdata('SIMULATIONS_T1_2.1s/markerAvrg_300.dat');
sim_full_InOut_pump300 =T.data(:,2) ;
sim_full_pump300 =T.data(:,3) ;
sim_roi_pump300 =T.data(:,6) ;

T = importdata('SIMULATIONS_T1_2.1s/markerAvrg_400.dat');
sim_full_InOut_pump400 =T.data(:,2) ;
sim_full_pump400 =T.data(:,3) ;
sim_roi_pump400 =T.data(:,6) ;


%% OCT 12 DATA

%%% Full phantom
Oct12_sig_full_pump200 = [0.0024    0.0026    0.0030    0.0034    0.0035    0.0032    0.0027    0.0022    0.0018    0.0014    0.0014];
Oct12_sig_full_pump300 = [0.0045    0.0053    0.0054    0.0054    0.0051    0.0048    0.0040    0.0033    0.0028    0.0022    0.0020];
Oct12_sig_full_pump400 = [0.0066    0.0071    0.0076    0.0074    0.0069    0.0064    0.0055    0.0042    0.0034    0.0027    0.0026];

%%% Spherical Segment
Oct12_sig_roi_pump200 = [ 0.0034    0.0050    0.0067    0.0084    0.0118    0.0127    0.0117    0.0100    0.0078    0.0061    0.0053];
Oct12_sig_roi_pump300 = [ 0.0117    0.0159    0.0179    0.0208    0.0224    0.0219    0.0185    0.0145    0.0119    0.0093    0.0085];
Oct12_sig_roi_pump400 = [ 0.0222    0.0264    0.0316    0.0333    0.0322    0.0287    0.0237    0.0170    0.0139    0.0106    0.0105];


 pump = [200, 300, 400] ;
 pld=[0.1,0.3,0.5,0.7,1.0,1.3,1.7,2.2,2.7,3.2,3.5];
 smooth_pld=0.1:0.02:5.2;
 
 mycol=['g','b','r','k','m'];
 txt = {};  % Used for Legend generation
 tau = 1.5;  % Label Duration
 lambda =1;
 T1blood =2.1;


figure(1)
set(gcf, 'Position',  [200, 600, 1100, 650]);
 

txt = {};  % Used for Legend generation
%set(gcf, 'Position',  [200, 600, 1100, 650]);

for i=1:length(pump)
    
  data1 = eval(sprintf('Oct12_sig_roi_pump%d', pump(i)));
  sim1 = eval(sprintf('sim_roi_pump%d', pump(i)));
    
  
   h1=plot(t,sim1, '--','Color', mycol(i));
   set(h1,'LineWidth',2,'MarkerSize',13);
   title('Spherical Segment ROI');
   hold on 
   
   h=plot(pld+tau, data1, '+','Color', mycol(i));
   set(h,'LineWidth',2,'MarkerSize',13);
   h.Annotation.LegendInformation.IconDisplayStyle = 'off';
  
        
   tt = sprintf('Pump Rate = %d mL/min',pump(i));
   txt = [txt; tt];
   
   hold on 
end
   
   hold off
   grid on

  xlabel('Time (s)');
  ylabel('\Delta M /M0');
  xlim([0. 5.5]);
  ylim([0 0.05]);
  set(gca,'FontSize',24);   % axis font size

  
  y1=get(gca,'ylim');
  LL=line([1.5 1.5], y1);
  set(LL,'LineWidth',2);
   
  legend(txt);
  
 fname=fullfile('OVERLAY_SIM_DATA_Flow_vs_Time_ROI');
        print(strcat(fname),'-djpeg');  
  
     
    



figure(3)
txt = {};  % Used for Legend generation
% set(gcf, 'Position',  [200, 600, 1100, 650]);
set(gcf, 'Position',  [50, 900, 900, 550]);

for i=1:length(pump)
    
  data1 = eval(sprintf('Oct12_sig_roi_pump%d', pump(i)));
  sim1 = eval(sprintf('sim_roi_pump%d', pump(i)));
   
   
  % Plot simulations
  
   h1=plot(t,sim1, ':','Color', mycol(i));
   set(h1,'LineWidth',2.5,'MarkerSize',13);
%    title('Spherical Segment ROI');
   h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
   
   hold on 
   
  % Plot data 
  
   h=plot(pld+tau, data1, '+','Color', mycol(i));
   set(h,'LineWidth',3,'MarkerSize',13);
  % h.Annotation.LegendInformation.IconDisplayStyle = 'off';
  
   % Fit Buxton Model
   ft = fittype('DeltaM(x,lambda,tau,T1blood,p1,p2)', 'problem', {'lambda','tau','T1blood'},'coefficients',{'p1','p2'}); 
   [ff, gof]= fit((pld+tau)',(data1)',ft,'StartPoint',[0.3,0.05],'problem',{lambda,tau,T1blood});
  
   
   %%%% Use Fit parameters to plot Buxton Model. 
   
   DM =DeltaM(smooth_pld,lambda,tau,T1blood,ff.p1,ff.p2);
   
   hold on 
   hf=plot(smooth_pld, DM,'LineWidth',2.5,'LineStyle', '--', 'Color', mycol(i)); % plot best fit curve, shifted x-axis values
   hf.Annotation.LegendInformation.IconDisplayStyle = 'off';
   
        
   tt = sprintf('Pump Rate = %d mL/min',pump(i));
   txt = [txt; tt];
   
   hold on 
end
   
   hold off
   grid on

  xlabel('Time (s)');
  ylabel('\Delta M /M0');
  xlim([0. 5.5]);
  ylim([0 0.05]);
  set(gca,'FontSize',24);   % axis font size

  
  y1=get(gca,'ylim');
  LL=line([1.5 1.5], y1);
  set(LL,'LineWidth',2);
   
  legend(txt);
  
 fname=fullfile('OVERLAY_SIM_DATA_GKM_Flow_vs_Time_ROI');
        print(strcat(fname),'-djpeg');  
  
     
        

