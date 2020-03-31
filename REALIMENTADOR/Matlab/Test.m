clear all
close all
clc



%% 

s = tf('s');

AD817 = db2mag(73.9)*(s/60E6/2/pi+1)*(s/100E6/2/pi+1)/(s/120E6/2/pi+1)^2/(s/160E6/2/pi+1)/(s/8.8E3/2/pi+1);

figure
h = bodeplot(AD817);
setoptions(h,'FreqUnits','Hz');
hold on
h = bodeplot(feedback(AD817,1/100));
setoptions(h,'FreqUnits','Hz');


R = 8;
C = 330E-9;
L = 33E-6;

