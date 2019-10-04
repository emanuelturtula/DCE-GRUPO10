close all
clear all
clc
set(0,'defaulttextinterpreter','latex');
set(0,'defaultlegendinterpreter','latex');

%%

T = dlmread('Sim.CSV', ';', 1, 0);
x=T(:,1);

%%
y=T(:,2);
figure
plot(x*1E6,y*1E3)
grid minor
xlim([10 25]);
xlabel("Tiempo [$\mu$s]")
ylabel("Corriente [mA]")
print('sim_corriente.eps', '-depsc', '-tiff');

%%
y=T(:,5);
figure
plot(x*1E6,y)
grid minor
xlim([10 25]);
xlabel("Tiempo [$\mu$s]")
ylabel("Tensi\'on [V]")
print('sim_tension.eps', '-depsc', '-tiff');