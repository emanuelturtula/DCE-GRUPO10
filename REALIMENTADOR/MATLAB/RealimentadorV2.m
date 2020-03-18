clear all
close all
clc



s = tf('s');
p_tl082 = 1/(s/4E6/2/pi + 1);

%200 ns de delay 
%POR EL MOMENTO NO LO USO
delay = exp(-s*200E-9);              

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Filtro LC
R = 8;
L = 33E-6;
C = 300E-9;
Filtro = (1/L/C)/(s^2 + s*1/R/C + 1/L/C);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ganancia a lazo abierto
a = 20*Filtro;                

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Restador con ganancia = 100
K_error = 82*p_tl082;              

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Realimentador:

%Etapa 1: bajar la tension de salida por un factor de 30. Requiere un
%operacional

feed = 1/20*p_tl082;

%Compensador
grados = 45;
phi = 2*pi*grados/360;
freq_max = 5E5;

alpha = (sin(phi)+1)/(1-sin(phi));
tau = 1/(freq_max*2*pi*sqrt(alpha));

comp = 1/alpha*(s*alpha*tau+1)/(s*tau+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ganancia de lazo
af(1) = K_error*a*feed; %Sin compensar
af(2) = af(1)*comp; %Con compensador

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAFICO

%SIN COMPENSAR
figure                                            
[Gm, Pm, Fcg, Fcp] = plot_bode(af(1));
grid minor
txt1 = sprintf('Margen de ganancia = %.2f  (%.1i Hz)', Gm, Fcg);
txt2 = sprintf('Margen de fase = %.2f ^{\\circ} (%.1i Hz)', Pm, Fcp);
title({txt1, txt2}, 'Fontsize', 12);

%CON COMPENSADOR
figure                                            
[Gm, Pm, Fcg, Fcp] = plot_bode(af(2));
hold on
plot_bode(comp);
grid minor
txt1 = sprintf('Margen de ganancia = %.2f  (%.1i Hz)', Gm, Fcg);
txt2 = sprintf('Margen de fase = %.2f ^{\\circ} (%.1i Hz)', Pm, Fcp);
title({txt1, txt2}, 'Fontsize', 12);




return
%%
syms R1 R2 C;
eqn1 = alpha == (R1+R2)/R2;
eqn2 = tau == R1/alpha*C;
eqn3 = C == 100E-12;

sol = solve([eqn1, eqn2, eqn3], [R1 R2 C]);

vpa(sol.R1)
vpa(sol.R2)
vpa(sol.C)


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









function [Gm, Pm, Fcg, Fcp] = plot_bode(func)
    w = logspace(0,10,50E3);
    [g, p] = bode(func, w);
    [Gm, Pm, Wcg, Wcp] = margin(g, p, w);
    Gm = mag2db(Gm);
    Fcg = Wcg/2/pi;
    Fcp = Wcp/2/pi;
    
    %Plot
    P = bodeoptions; P.FreqUnits = 'Hz';
    bodeplot(func, w, P);    
end

function indexes = buscar_armonicos(f, fundamental, n)
    epsilon = 0.5;
    for i=1:n
        indexes(i)=find(abs(f-i*fundamental)<epsilon);
    end
end

