close all
clear all
clc
set(0,'defaulttextinterpreter','latex');
set(0,'defaultlegendinterpreter','latex');

%% Intermodulación

% Primero leer los archivos. Adentro estan las señales obtenidas para 1W y 
% 30W. Solo se simuló un período porque la simulación tarda muhco si no.
% Así que hay que repetir la señal con el MATLAB para poder hacer una mejor
% FFT.

signal_1W = dlmread("Intermod_1W.txt", "\t", 2, 0);
signal_30W = dlmread("Intermod_30W.txt", "\t", 2, 0);

vo_1W = signal_1W(:,2); %Solo me interesa Vo
vo_30W = signal_30W(:,2);

% Considero que 16 períodos es suficiente. Algo muy importante es que al
% mirar el vector de tiempo de los archivos, se ve que los puntos no estan
% espaciados de forma equidistante, sino que va variando según vaya
% resolviendo la simulación. Esto no nos sirve, así que vamos a reacomodar
% el vector de forma que el tiempo de muestreo sea fijo, sabiendo que el
% período para una señal de 100Hz es de 10ms.

N = 16; %Cantidad de períodos
T = 10E-3; %Período

                                                
for i=1:4
    vo_1W = vertcat(vo_1W, vo_1W);              %Esta función pega las 
                                                %señales una a continuación
                                                %de la otra. En la primer
                                                %iteración (i=1), se
                                                %obtendrán 2 períodos, en
                                                %la segunda (i=2), 4
                                                %períodos, y en (i=4), 16
                                                %períodos
    vo_30W = vertcat(vo_30W, vo_30W);
end

time_1W = 0:T*N/(length(vo_1W)-1):T*N;          %15 períodos de 10ms
time_30W = 0:T*N/(length(vo_30W)-1):T*N;        %Si bien los tiempos de 
                                                %simulación van a ser 
                                                %iguales, la cantidad de 
                                                %puntos obtenidos durante 
                                                %la simulación no tienen 
                                                %por qué ser iguales.
                                                
figure
plot(time_1W*1E3, vo_1W, '-k', 'DisplayName', 'Potencia: 1W');
hold on
plot(time_30W*1E3, vo_30W, '-r', 'DisplayName', 'Potencia: 30W');

title("\begin{tabular}{c} Lazo abierto \\ Distorsi\'on por intermodulaci\'on \end{tabular}", 'Fontsize', 14)
legend('location', 'Southeast');
xlabel("Tiempo [ms]");
ylabel("Tensi\'on [V]");
xlim([0 80]);
grid minor

%Ahora la FFT. El intervalo de muestreo nuevamente es distinto para ambas
%señales porque la cantidad de puntos entre ellas difiere. La frecuencia de
%muestreo es la inversa, entonces:

Fs_1W = 1/(time_1W(end)-time_1W(end-1));
%Fs_30W = 1/(time_30W(end)-time_30W(end-1));

L_1W = length(vo_1W);
f_1W = 0: Fs_1W/L_1W : Fs_1W/2;  %Vector de frecuencia para 1W
Y_1W = fft(vo_1W)/L_1W;
Y_1W = 2*abs(Y_1W( 1 : L_1W/2+1));  %Vector de amplitud

figure 
semilogx(f, 20*log10(Y_1W), '-b');




