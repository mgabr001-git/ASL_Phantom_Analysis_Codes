 close all
 clear
 clc
 
%REF	 = [22.87	27.15	31.93	34.08	36.62	39.2	40.86	39.22	39.95	40.47	45.85];  % OCT 12 dataset

%%%% 4 dataset avg REF and STDERR between 4 datasets
REF  = [23.2836	27.81035	31.046225	33.586125	34.962	37.2379	38.840375	39.54065	42.115125	43.913	47.179875 ];
STDERR =[1.004206625	0.953823241	1.261226642	1.079465322	1.483642273	1.517240239	1.554616731	1.59141783	2.672587869	2.78954313	2.956834865];
   


ASYM = [12.99	14.68	15.13	15.88	15.15	18.36	20.2	19.6	19.54	30.68	25.26];
ASYM_STDERR = [1.0488	0.34135	1.0625	1.2751	1.77295	0.74815	0.68325	0.46225	2.03715	0.56625	1.21285];

GmaxGavg4 = [13.95	17.13	20.59	20.9	21.25	24.64	24.84	27.32	26.13	31.83	36.37];
GmaxGavg4_STDERR = [1.269	0.59475	0.38615	0.084	0.7974	0.3551	0.8644	0.85165	1.95865	3.41245	3.5818];

FA20 = [15.55	20.01	24.31	25.26	26.07	27.02	29.62	32.95	31.76	34.83	34.85];
FA20_STDERR = [ 1.38285	0.56335	1.8746	2.3078	0.6471	0.3818	2.166	1.82925	1.31395	2.8589	1.0104];


PCASL_exp = ["ASYM", "GmaxGavg4", "FA20"];

 pump = 400 ;
 pld=[0.1,0.3,0.5,0.7,1.0,1.3,1.7,2.2,2.7,3.2,3.5];
 mycol=['g','b','r','k','m'];
 tau = 1.5;  % Label Duration
 
 
 A = REF./ASYM;
 B = REF./GmaxGavg4;
 C = REF./FA20;
 
% mean(A(1:9))
% mean(B(1:9))
% mean(C(1:9))
 
 
 
A = (REF-ASYM) ./ REF ;
B = (REF-GmaxGavg4) ./ REF ;
C = (REF-FA20) ./ REF  ;


mean(A(1:9))*100
mean(B(1:9))*100
mean(C(1:9))*100
%  
% sprintf("###############################################################################")
% 

% 
e_ASYM = (1-mean(A(1:9)))*85
e_GmaxGavg = (1-mean(B(1:9)))*85
e_FA20 = (1-mean(C(1:9)))*85

 
figure(1);


set(gcf, 'Position',  [1200, 600, 1100, 650]);

for i=1:length(PCASL_exp)

  txt = {};  % Used for Legend generation
  
  sig = eval(sprintf(PCASL_exp(i)));
  err = eval(sprintf(strcat(PCASL_exp(i),'_STDERR')));
% Plot Reference
   h=errorbar(pld+tau,REF,STDERR, STDERR,'.-','Color', 'r');
   set(h,'LineWidth',2,'MarkerSize',13);
      
   hold on 
% Plot PCASL Experiment data   
   h=errorbar(pld+tau,sig ,err,err,' .-','Color', 'b');
     
   tt = sprintf(PCASL_exp(i));
   set(h,'LineWidth',2,'MarkerSize',13);
   txt = ['REF', tt];
   
  grid on

  xlabel('Time (s)');
  ylabel('Flow (mL/100g/min)');
  xlim([0. 5.5]);
  ylim([0 60]);
  set(gca,'FontSize',24);   % axis font size

  
   y1=get(gca,'ylim');
   LL=line([1.5 1.5], y1);
   set(LL,'LineWidth',2);
  
   
  legend(txt);
  
   hold off
  
   fname=fullfile(strcat(PCASL_exp(i),'_Flow_vs_Time_4DataSet_avg_stderr'));
   print(strcat(fname),'-djpeg');  
  
  
   
end





 