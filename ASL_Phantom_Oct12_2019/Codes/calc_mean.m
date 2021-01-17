clear; clc;


%%% Flow and ATT from Full phantom Fits

a= [19.72 21.12 16.72 17.20];
b = [29.29 29.61 25.80 27.06];
c = [37.62 39.93 33.25 34.61];


m = mean(a)
stderr = std(a)/2


m = mean(b)
stderr = std(b)/2


m = mean(c)
stderr = std(c)/2


t1 = [0.90 0.82 0.86 0.87];
t2 = [0.74 0.60 0.72 0.63];
t3 = [0.62 0.59 0.59 0.48];


m = mean(t1)
stderr = std(t1)/2


m = mean(t2)
stderr = std(t2)/2


m = mean(t3)
stderr = std(t3)/2




