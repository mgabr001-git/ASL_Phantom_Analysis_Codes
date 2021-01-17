 close all
 clear
 clc

 
 
 ROIflag =1;   %%%  Sphere-segment ROI = 1 , Full Phantom = 2
 

 
if ~ismember(ROIflag,[1 2])
    sprintf("Please select correct ROIflag: 1 = Sphere-segment ROI, 2 = Full Phantom !")
    return
elseif ROIflag == 1  
 %%%  SPHERE-SEGMENT ROI
 %%%  Mean Flow from 4 datasets at all PLDs 
 %%%  eerr here is the STDERR = STDEV/sqrt(4)   
  sig_pump200 = [11.415925    17.788925     26.9547      38.53625       57.950475       75.182125      88.652   	 96.6308      98.550225	    97.2817        96.142975];
  err_pump200 = [1.435445289	2.429397551  3.037805033   4.667409324    5.376725304     4.891375188	 4.23425488  4.391862484   3.938418363	 3.853709072   3.953456549];
   
  sig_pump300 = [43.438775	  59.547075	   76.8185	     95.1134	  115.211025	128.504325	  136.061525	     135.59675	    138.2555	137.638525      145.44405 ];
  err_pump300 = [4.329176928   4.882370949	6.297690278   5.804742364	5.13959589	  5.438299063	4.565573238	   3.570515928	  4.902958442	  5.863645644	  6.800181185 ];
   
   
  sig_pump400 = [78.1575	  103.637625      129.269        148.768075      163.690825     169.029775 	  168.9227   	162.02805	  169.492975    164.7001	  177.140725 ];
  err_pump400 = [5.440338333	6.946840197	    7.745073544    5.725058408	   6.387022674	  5.223607556	5.370197771   5.969465993	9.950375569	  7.823567207	6.589582781 ];

  fname=fullfile('ROI1_SphSeg_Flow_vs_Time_4DataSet_avg_stderr');
  
 elseif ROIflag == 2

 %%%  FULL PHANTOM
 %%%  Mean Flow from 4 datasets at all PLDs 
 %%%  eerr here is the STDERR = STDEV/sqrt(4) 
  sig_pump200 = [8.481675       10.069325	12.649925	15.360275	16.866425	18.36365	19.3735     20.4932     21.240125	21.637475	23.702675];
  err_pump200 = [0.542847169	0.688200629	0.738449041	1.123789476	1.357003881	0.892134207	0.565589162	1.222066377	1.262767078	1.04415708	1.165486866];
   
  sig_pump300 = [16.36915       19.55015	22.656475	24.55315	25.9385     27.701925	29.216325	31.40915	33.669575	33.0358     36.705825];
  err_pump300 = [ 0.848805225	0.898493397	1.107157452	0.897001056	0.80782588	1.240125327	1.240492094	0.712272319	1.419988383	0.712347567	1.530338354];
   
  
  sig_pump400 = [23.2836        27.81035	31.046225	33.586125	34.962      37.2379     38.840375	39.54065	42.115125	43.913  	47.179875];
  err_pump400 = [1.004206625	0.953823241	1.261226642	1.079465322	1.483642273	1.517240239	1.554616731	1.59141783	2.672587869	2.78954313	2.956834865];
   
  ylim([0 105]);
   fname=fullfile('ROI2_FullPhantom_Flow_vs_Time_4DataSet_avg_stderr');
   
 end        
 
  
 
  
 pump = [200, 300, 400] ;
 pld=[0.1,0.3,0.5,0.7,1.0,1.3,1.7,2.2,2.7,3.2,3.5];
 mycol=['g','b','r','k','m'];
 txt = {};  % Used for Legend generation
 tau = 1.5;  % Label Duration

 
 
  
figure(1);
set(gcf, 'Position',  [1200, 600, 1100, 650]);

 

for i=1:length(pump)
    
  sig = eval(sprintf('sig_pump%d', pump(i)));
  err = eval(sprintf('err_pump%d', pump(i)));
  stdev = err *2;
  
  % h=errorbar(pld+tau,sig ,stdev, stdev,'+-','Color', mycol(i));
   
   h=errorbar(pld+tau,sig ,err,err,'+-','Color', mycol(i));
     
   tt = sprintf('Pump Rate = %d mL/min',pump(i));
   set(h,'LineWidth',2,'MarkerSize',13);
   txt = [txt; tt];
   
   hold on 
  
end
   
   grid on

  xlabel('Time (s)');
  ylabel('Flow (mL/100g/min)');
  xlim([0. 5.5]);
%  ylim([0 65]);
  set(gca,'FontSize',24);   % axis font size

  
  y1=get(gca,'ylim');
  LL=line([1.5 1.5], y1);
  set(LL,'LineWidth',2);
   
 
  if ROIflag == 1
    ylim([0 230]);
  
  elseif ROIflag == 2  
      xl=get(gca,'xlim');
     L20=line(xl, [20. 20.]);
     L30=line(xl, [30. 30.]);
     L40=line(xl, [40. 40.]);
     set(L20,'LineWidth',2, 'Color', 'g', 'LineStyle','--');  
     set(L30,'LineWidth',2, 'Color', 'b', 'LineStyle','--');  
     set(L40,'LineWidth',2, 'Color', 'r', 'LineStyle','--');  
  
      %%% Do not include line in the legend
    %  L20.Annotation.LegendInformation.IconDisplayStyle = 'off';  
      ylim([0 65]);
      
  end



   legend(txt);
 
  
 %  text(0.25, 42, 'Nominal Flow','Color','k','FontSize',24);
   
    print(strcat(fname),'-djpeg');  
  
