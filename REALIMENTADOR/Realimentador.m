clear all
close all
clc
set(0,'defaulttextinterpreter','latex');
set(0,'defaultlegendinterpreter','latex');

R = 8;
L = 33E-6;
C = 300E-9;

Vtriang = 1.5;
Vsource = 30;               %VALORES DEL AMPLIFICADOR

s = tf('s');
w = logspace(0,10,50E3);

%% Transferencias lazo abierto


td_ol = 87E-9;                  %Delay introducido el amplificador, 
                                %desde entrada de PWM a salida de PWM
Delay = exp(-s*td_ol);
Kpwm = Vsource/Vtriang;         %Ganancia (linealización de PWM)

Filtro = (1/L/C)/(s^2 + s*1/R/C + 1/L/C);   %TRANSFERENCIA DEL FILTRO
[numF, denF] = tfdata(Filtro, 'v');  %NUMERADOR Y DENOMINADOR DEL FILTRO

Planta = Kpwm*Delay*Filtro;     %Planta


[mag_ol, fase_ol] = bode(Planta, w);

figure
semilogx(w/2/pi, mag2db(squeeze(mag_ol)));
title ("Transferencia a lazo abierto");
grid minor
ylabel("Amplitud (dB)")
xlabel("Frecuencia [Hz]")
xlim([1E3 10E6])

%% Ganancia de lazo

                       
p_tl082 = 1/(1+s/(4E6*2*pi)); %Hay dos operacionales: uno es el amplificador 
                            %de error y otro es el seguidor en el
                            %realimentador

%f = 1/Kpwm;    %Realimento con constante 1/a_ol. Agrego el polo del
                       %operacional (seguidor)
Kerror = 15;   %Pruebo con amplificador de error de ganancia 26dB 
                       %(20 en veces) y pongo el polo del operacional
                       %(4 MHz)
%G_error = Kerror*p_tl082; 
G_error = Kerror*p_tl082; 
                       
% GH1 = Planta*f*G_error;  %Ganancia de lazo
% 
% [mag_GH1, fase_GH1] = bode(GH1, w);
% [Gm, Pm, Wcg, Wcp] = margin(mag_GH1, fase_GH1, w);
% freq_pm1 = Wcp/2/pi;
% 
% figure
% yyaxis left
% semilogx(w/2/pi, mag2db(squeeze(mag_GH1)));
% ylabel("Amplitud (dB)")
% hold on
% 
% yyaxis right
% semilogx(w/2/pi, squeeze(fase_GH1));
% plot([Wcp/2/pi Wcp/2/pi], [-360 180])
% plot([400E3 400E3], [-360 180])
% ylim([-360 180])
% txt1 = "Ganancia de lazo sin compensar";
% txt2 = sprintf('Margen de ganancia = %.2f dB', mag2db(Gm));
% txt3 = sprintf('Margen de fase = %.2f $^{\\circ}$', Pm);
% title ({txt1, txt2, txt3});
% grid minor
% xlabel("Frecuencia [Hz]")
% xlim([1E3 10E6])

%% Compensador

% RED DE ADELANTO
syms alpha tau;
% wmax = Wcp;
% grados = 45;
% phimax = 2*pi*grados/360;
% 
% eqn = sin(phimax) == (1-alpha)/(alpha+1);
% alpha = double(solve(eqn, alpha));
% 
% eqn = 1/tau/sqrt(alpha) == wmax;
% tau = double(solve(eqn, tau));
% Kgc = alpha*5;

z = 3E5*2*pi;
p = 4E6*2*pi;
Kgc = 1/20;

Gc = Kgc*(1+s/z)/(1+s/p)*p_tl082^2;
[mag_Gc, fase_Gc] = bode(Gc, w);


GH2 = Planta*G_error*Gc;

[mag_GH2, fase_GH2] = bode(GH2, w);
[Gm, Pm, Wcg, Wcp] = margin(mag_GH2, fase_GH2, w);         %BODE
 

figure
yyaxis left
semilogx(w/2/pi, mag2db(squeeze(mag_GH2)), '-b');
hold on
semilogx(w/2/pi, mag2db(squeeze(mag_Gc)), '-k');
ylabel("Amplitud (dB)")

yyaxis right
semilogx(w/2/pi, squeeze(fase_GH2), '-r');
semilogx(w/2/pi, squeeze(fase_Gc), '-m');
ylim([-360 180])
txt1 = "Ganancia de lazo con compensador";
txt2 = sprintf('Margen de ganancia = %.2f dB', mag2db(Gm));
txt3 = sprintf('Margen de fase = %.2f $^{\\circ}$', Pm);
title ({txt1, txt2, txt3});
grid minor
xlabel("Frecuencia [Hz]")
xlim([1E3 10E6])
% 
% a = 10E3;   %R1 + R2
% R2 = a/alpha;
% R1 = a-R2;
% C1 = tau*alpha/R1;
% 
% 
% fprintf("CON RED DE ADELANTO\n");
% fprintf("In    R1=%d                  Out   \n", R1);
% fprintf("----/\\/\\/\\-------------------\n");  
% fprintf("  |        |           |\n");
% fprintf("  |---||---|           / \n");
% fprintf("      C1=%d  \\ R2=%d\n", C1, R2);
% fprintf("                       / \n");
% fprintf("                       \\\n");
% fprintf("                       | \n");
% fprintf("-----------------------------\n\n\n");
% 
% fprintf("Margen de fase: %d\n", Pm);
% fprintf("Margen de ganancia: %d\n\n\n", Gm);
% 
% 
% R1 = 5.1E3;
% R2 = 4.7E3;
% C1 = 220E-12;
% 
% tau = R1*R2/(R1+R2)*C1;         %RED CON VALORES COMERCIALES
% alpha = (R1+R2)/R2;
% 
% Kh = db2mag(-26); %Ganancia del realimentador
% Kgc = 1/alpha;
% 
% Gc = Kgc*(1+s*alpha*tau)/(1+s*tau);
% 
% Kh = Kh*Gc;
% 
% [mag, phase] = bode(Kerror*Planta*Kh, w);
% [Gm, Pm, Wcg, Wcp] = margin(mag, phase, w);         %BODE
% 
% figure
% yyaxis left
% semilogx(w/2/pi, mag2db(squeeze(mag)));
% ylim([-50 50])
% hold on 
% yyaxis right
% semilogx(w/2/pi, squeeze(phase));
% semilogx([Wcp/2/pi Wcp/2/pi], [-500 500]);
% ylim([-360 180])
% xlim([1 10^7/2/pi])
% grid minor
% 
% fprintf("CON VALORES COMERCIALES: \n")
% fprintf("Margen de fase: %d\n", Pm);
% fprintf("Margen de ganancia: %d\n\n\n", Gm);



