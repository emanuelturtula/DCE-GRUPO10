close all
clear all
clc
set(0,'defaulttextinterpreter','latex');
set(0,'defaultlegendinterpreter','latex');

%%
vo = readtable('VO.CSV');
x=table2array(vo(:,4));
x = x + abs(x(1));
y=(table2array(vo(:,5))+6.08)/2;  %6.08 por el desplazamiento vertical y /2 por la escala

figure
plot(x*1E6,y)
grid minor
xlim([0 80])
ylim([0 15])
xlabel("Tiempo [$\mu$s]")
ylabel("Tensi\'on [V]")

Fs = 1.5E6;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = length(y);             % Length of signal
t = (0:L-1)*T;        % Time vector

Y = fft(y);
Y = abs(Y/L);
Y = Y(1:L);
Y(2:end-1) = 2*Y(2:end-1);
f = (0:Fs)/L;

figure
plot(f*1E-3,20*log(Y)) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (kHz)')
ylabel('|P1(f)|')

%%
% io = readtable('IO.CSV');                     ESTE DA MUCHO MAS FEO QUE IO2
% x=table2array(io(:,4));                       POR ESO LO COMENTO
% x = x + abs(x(1));
% y=(table2array(io(:,5)))/0.05;
% 
% figure 
% plot(x*1E6,y)
% grid minor
% xlim([0 80])
% xlabel("Tiempo [$\mu$s]")
% ylabel("Corriente [mA]")

%%
io2 = readtable('IO2.CSV');
x=table2array(io2(:,4));
x = x + abs(x(1));
y=(table2array(io2(:,5))+0.14900)/0.05;

figure 
plot(x*1E6,y)
grid minor
xlim([0 20])                        %HASTA 80us le puedo dar
xlabel("Tiempo [$\mu$s]")
ylabel("Corriente [mA]")

