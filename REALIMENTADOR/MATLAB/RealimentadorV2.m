close all
clear all
clc

s = tf('s');          

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
%ATENCION: SI SE MODIFICA LA GANANCIA, HAY QUE RECALCULAR EL POLO
%Restador
fc = 387E3;
g = 100;
Restador = g/(s/fc/2/pi+1);              

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Amplificador de error 
%ATENCION: SI SE MODIFICA LA GANANCIA, HAY QUE RECALCULAR EL POLO
fc = 4.2E6;
g = 5.1;
G_error = g/(s/fc/2/pi+1);

% figure                                            
% plot_bode(G_error);
% grid minor
% title("Bode de amplificador de error", 'Fontsize', 12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Realimentador:
feed = 1/20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ganancia de lazo sin compensar
af(1) = Restador*G_error*a*feed;

figure                                            
[Gm, Pm, Fcg, Fcp] = plot_bode(af(1));
grid minor
txt1 = "Sin compensar";
txt2 = sprintf('Margen de fase = %.2f ^{\\circ} (%.1i Hz)', Pm, Fcp);
title({txt1, txt2}, 'Fontsize', 12);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Ganancia de lazo con compensador

%Compensador (tiene un operacional con ganancia = alpha, polo en )
grados = 45;
phi = 2*pi*grados/360;
freq_max = 7.6E5;
alpha = (sin(phi)+1)/(1-sin(phi));
tau = 1/(freq_max*2*pi*sqrt(alpha));
comp = 1/alpha*(s*alpha*tau+1)/(s*tau+1);

%Ganancia
af(2) = af(1)*comp; 

figure                                            
[Gm, Pm, Fcg, Fcp] = plot_bode(af(2));
hold on
plot_bode(comp);
grid minor
txt1 = sprintf('Margen de ganancia = %.2f  (%.1i Hz)', Gm, Fcg);
txt2 = sprintf('Margen de fase = %.2f ^{\\circ} (%.1i Hz)', Pm, Fcp);
title({txt1, txt2}, 'Fontsize', 12);

%% COMPONENTES
syms R1 R2 C;
eqn1 = alpha == (R1+R2)/R2;
eqn2 = tau == R1/alpha*C;
eqn3 = C == 100E-12;

sol = solve([eqn1, eqn2, eqn3], [R1 R2 C]);

R1 = vpa(sol.R1);
R2 = vpa(sol.R2);


%% Funciones
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

