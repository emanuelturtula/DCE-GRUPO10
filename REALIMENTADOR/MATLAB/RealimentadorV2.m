clear all
close all
clc

R = 8;
L = 33E-6;
C = 300E-9;

s = tf('s');

Filtro = (1/L/C)/(s^2 + s*1/R/C + 1/L/C);   %TRANSFERENCIA DEL FILTRO
a = 20*Filtro;
plot_bode(a, 1)

%Compenso primero y después ajusto las ganancias:
phi = 2*pi*55/360;
wmax = 2*pi*2.27E5;
alpha = (1+sin(phi))/(1-sin(phi));
tau = 1/(wmax*sqrt(alpha));

lead = (s*alpha*tau+1)/(s*tau+1);

plot_bode(a*lead, 1);
plot_bode(lead, 0);



%Ganancia de lazo
Kerror = 100;
Kfeed = 1/20;
feed = lead*Kfeed;

af = a*feed*Kerror;

plot_bode(af, 1);


return




















%% SIMULINK
sim_output = sim('Realimentador_SIM.slx', ...
            'SimulationMode','normal', 'AbsTol','1e-5',...
            'SaveState','on','StateSaveName','xoutNew',...
            'SaveOutput','on','OutputSaveName','youtNew');
        
vo_ol = sim_output.Vo_ol.signals(1).values;
t = sim_output.Vo_ol.time;

T = t(2)-t(1);        % Sampling period
Fs = 1/T;             % Sampling frequency                     
L = length(vo_ol);       % Length of signal
t = (0:L-1)*T;        % Time vector



VO = abs(fft(vo_ol)/L);
VO = VO(1:L/2+1);
VO(2:end-1) = 2*VO(2:end-1);


f = Fs*(0:(L/2))/L;
figure
semilogx(f, mag2db(VO)) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
grid minor
xlim([1E4 1E7])

idx = buscar_armonicos(f, 20E3, 3);
THD_ol = rssq(VO(idx(2:end)))/VO(idx(1));









function plot_bode(func, create_new)
    w = logspace(0,10,50E3);
    if (create_new)
        figure
    else 
        hold on
    end    
    P = bodeoptions; P.FreqUnits = 'Hz';
    bodeplot(func, w, P);
    grid minor
end

function indexes = buscar_armonicos(f, fundamental, n)
    epsilon = 0.5;
    for i=1:n
        indexes(i)=find(abs(f-i*fundamental)<epsilon);
    end
end

