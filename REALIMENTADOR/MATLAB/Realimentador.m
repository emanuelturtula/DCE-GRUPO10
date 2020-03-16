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

epsilon = 1E2;
f_200k = find(abs(w/2/pi-200E3)<epsilon);
%% Transferencias lazo abierto

td_ol = 87E-9;                  %Delay introducido el amplificador, 
                                %desde entrada de PWM a salida de PWM
Delay = exp(-s*td_ol);
Kpwm = Vsource/Vtriang;         %Ganancia (linealización de PWM)

Filtro = (1/L/C)/(s^2 + s*1/R/C + 1/L/C);   %TRANSFERENCIA DEL FILTRO
[numF, denF] = tfdata(Filtro, 'v');  %NUMERADOR Y DENOMINADOR DEL FILTRO

Planta = Kpwm*Delay*Filtro;     %Planta
[mag_P, fase_P] = bode(Planta, w);

figure
yyaxis left
semilogx(w/2/pi, mag2db(squeeze(mag_P)));
ylabel("Amplitud [dB]")
hold on

yyaxis right
semilogx(w/2/pi, squeeze(fase_P));
ylim([-360 180])
txt1 = "Planta";
title (txt1);
grid minor
ylabel("Fase")
xlabel("Frecuencia [Hz]")
xlim([20 10E6])

%% Ganancia de lazo

                       
p_tl082 = 1/(1+s/(4E6*2*pi)); %POLO DEL TL082
                
A = 1;
f = 130E3;
Q = 0.5771;
k = 1.2754;
C1 = 100E-12;
filtro = get_filter(A, f, Q, k, C1);    %OBTENGO FILTRO PASABAJOS CON 
                                        %FRECUENCIA DE CORTE DE 130 kHz

% feed = filtro*1/Kpwm*p_tl082;          %Realimentador: filtro+divisor de 
                                       %tension.       
feed = 1/Kpwm*p_tl082;
[mag_R, fase_R] = bode(feed, w);     

%Tension a 200k a la salida del realimentador:
v_vtriang = Vtriang*mag_P(f_200k(1))*mag_R(f_200k(1));

%Planta es el amplificador a lazo abierto. 1/Kpwm es un divisor de tensión 
%colocado justo después de la salida, en el realimentador. El filtro es un 
%filtro pasabajos para eliminar lo mas posible las frecuencias de 130 kHz 
%para arriba. Kerror es el restador con ganancia para obtener la señal de
%error. El filtro tiene un operacional, al igual que el restador. Además, 
%el filtro invierte la fase 180 grados, así que se agrega otro inversor (un
%operacional más (3 operacionales en total).

Kerror = 82*p_tl082;

GH1 = Planta*feed*Kerror;   %Ganancia de lazo sin compensar  

[mag_GH1, fase_GH1] = bode(GH1, w);
[Gm, Pm, Wcg, Wcp] = margin(mag_GH1, fase_GH1, w);

figure
yyaxis left
semilogx(w/2/pi, mag2db(squeeze(mag_GH1)));
ylabel("Amplitud [dB]")
hold on

yyaxis right
semilogx(w/2/pi, squeeze(fase_GH1));
plot([Wcp/2/pi Wcp/2/pi], [-360 180])
ylim([-360 180])
txt1 = "Ganancia de lazo sin compensar";
txt2 = sprintf('Margen de ganancia = %.2f dB', mag2db(Gm));
txt3 = sprintf('Margen de fase = %.2f $^{\\circ}$', Pm);
title ({txt1, txt2, txt3});
grid minor
ylabel("Fase")
xlabel("Frecuencia [Hz]")
xlim([20 10E6])

%Compensando:
grados = 40;
phi = 2*pi*grados/360;
freq_max = Wcp/2/pi;
gain = 1;
[Gc, alpha, tau] = get_lead(phi, freq_max, gain);
GH2 = GH1*Gc;

[mag_GH2, fase_GH2] = bode(GH2, w);
[Gm, Pm, Wcg, Wcp] = margin(mag_GH2, fase_GH2, w);
[mag_Gc, fase_Gc] = bode(Gc, w);

figure
yyaxis left
semilogx(w/2/pi, mag2db(squeeze(mag_GH2)), '-b');
hold on
semilogx(w/2/pi, mag2db(squeeze(mag_Gc)), '--k');
ylabel("Amplitud [dB]")


yyaxis right
semilogx(w/2/pi, squeeze(fase_GH2), '-r');
semilogx(w/2/pi, squeeze(fase_Gc), '--g');
plot([Wcp/2/pi Wcp/2/pi], [-360 360])
ylim([-360 180])
txt1 = "Ganancia de lazo con compensador";
txt2 = sprintf('Margen de ganancia = %.2f dB', mag2db(Gm));
txt3 = sprintf('Margen de fase = %.2f $^{\\circ}$', Pm);
title ({txt1, txt2, txt3});
grid minor
ylabel("Fase")
xlabel("Frecuencia [Hz]")
xlim([20 10E6])

%% FUNCIONES
function H = get_filter(A, fc, Q, k, C1)
s=tf('s');
R3 = Q*(A+1)/(pi*k*fc*C1);
C2 = 1/(4*pi*k*fc*Q*R3);
R1 = R3/A;
R2 = R3/(A+1);

a1 = 1/C1/C2/R1/R2;
b1 = 1/C1*(1/R1+1/R2+1/R3);
b2 = 1/C1/C2/R2/R3;

H = -a1/(s^2+s*b1+b2);
end

function [comp, alpha, tau] = get_lead(phi, freq_max, gain)
s = tf('s');

alpha = (sin(phi)+1)/(1-sin(phi));
tau = 1/(freq_max*2*pi*sqrt(alpha));

comp = gain/alpha*(1+s*alpha*tau)/(1+s*tau);
end

