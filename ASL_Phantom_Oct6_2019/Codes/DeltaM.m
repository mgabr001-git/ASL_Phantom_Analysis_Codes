function dM=DeltaM(t,lambda,tau,T1blood,p1,p2)

% Perfusion Kinetic Model (Buxton et al 1998 )
% asl struct defines the model parameter values, p1 = transit time, p2 =
% flow 


dM = zeros(size(t));

    T1prime= T1blood*lambda/(lambda+T1blood*p2);    %  1/T1prime = 1/T1blood+f/lambda

  for i = 1:length(t)   

     if 0 < t(i) && t(i) < p1 
         dM(i) = 0;
     elseif p1 <= t(i) && t(i) < tau+p1
         dM(i) = 2*p2*0.85*T1prime*exp(-p1/T1blood)*(1-exp(-(t(i)-p1)/T1prime));
    
     elseif tau+p1 <= t(i)
         dM(i) = 2*p2*0.85*T1prime*exp(-p1/T1blood)*exp(-(t(i)-tau-p1)/T1prime)*(1-exp(-tau/T1prime)); 

     end
   end  
     
 
end






