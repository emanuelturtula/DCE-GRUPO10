close all
clear all
clc
set(0,'defaulttextinterpreter','latex');
set(0,'defaultlegendinterpreter','latex');

%%

vo = dlmread('VO.CSV',',',10,0);

x=vo(:,4);
x = x + abs(x(1)); %desplazo el tiempo negativo a 0
y = vo(:,5);  

figure
plot(x*1E6,y)
grid minor
xlim([10 25])
ylim([2 14])
xlabel("Tiempo [$\mu$s]")
ylabel("Tensi\'on [V]")

print('med_tension.eps', '-depsc', '-tiff');


% figure
% plot(abs())
% grid minor
% xlim([10 25])
% ylim([2 14])
% xlabel("Tiempo [$\mu$s]")
% ylabel("Tensi\'on [V]")
% 
% print('med_tension.eps', '-depsc', '-tiff');

%%
io = dlmread('IO2.CSV',',',10,0);

x= io(:,4);
x = x + abs(x(1)); %desplazo el tiempo negativo a 0
y = (io(:,5)-0.149)/10;   

figure
plot(x*1E6,y*1E3)
grid minor
xlim([10 25])
% ylim([-10 8])
xlabel("Tiempo [$\mu$s]")
ylabel("Corriente [mA]")

print('med_corriente.eps', '-depsc', '-tiff');

