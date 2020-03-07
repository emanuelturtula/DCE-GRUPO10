close all
clear all
clc
set(0,'defaulttextinterpreter','latex');
set(0,'defaultlegendinterpreter','latex');

%% Intermodulaci�n

% Primero leer los archivos. Adentro estan las se�ales obtenidas para 1W y 
% 30W. Solo se simul� un per�odo porque la simulaci�n tarda muhco si no.
% As� que hay que repetir la se�al con el MATLAB para poder hacer una mejor
% FFT.

signal_1W = dlmread("Intermod_1W.txt", "\t", 2, 0);
signal_30W = dlmread("Intermod_30W.txt", "\t", 2, 0);

vo_1W = signal_1W(:,2); %Solo me interesa Vo
vo_30W = signal_30W(:,2);

% Considero que 16 per�odos es suficiente. Algo muy importante es que al
% mirar el vector de tiempo de los archivos, se ve que los puntos no estan
% espaciados de forma equidistante, sino que va variando seg�n vaya
% resolviendo la simulaci�n. Esto no nos sirve, as� que vamos a reacomodar
% el vector de forma que el tiempo de muestreo sea fijo, sabiendo que el
% per�odo para una se�al de 100Hz es de 10ms.

N = 16; %Cantidad de per�odos
T = 10E-3; %Per�odo

                                                
for i=1:4
    vo_1W = vertcat(vo_1W, vo_1W);              %Esta funci�n pega las 
                                                %se�ales una a continuaci�n
                                                %de la otra. En la primer
                                                %iteraci�n (i=1), se
                                                %obtendr�n 2 per�odos, en
                                                %la segunda (i=2), 4
                                                %per�odos, y en (i=4), 16
                                                %per�odos
    vo_30W = vertcat(vo_30W, vo_30W);
end

time_1W = 0:T*N/(length(vo_1W)-1):T*N;          %15 per�odos de 10ms
time_30W = 0:T*N/(length(vo_30W)-1):T*N;        %Si bien los tiempos de 
                                                %simulaci�n van a ser 
                                                %iguales, la cantidad de 
                                                %puntos obtenidos durante 
                                                %la simulaci�n no tienen 
                                                %por qu� ser iguales.
                                                
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
%se�ales porque la cantidad de puntos entre ellas difiere. La frecuencia de
%muestreo es la inversa, entonces:

Fs_1W = 1/(time_1W(end)-time_1W(end-1));
%Fs_30W = 1/(time_30W(end)-time_30W(end-1));

L_1W = length(vo_1W);
f_1W = 0: Fs_1W/L_1W : Fs_1W/2;  %Vector de frecuencia para 1W
Y_1W = fft(vo_1W)/L_1W;
Y_1W = 2*abs(Y_1W( 1 : L_1W/2+1));  %Vector de amplitud

figure 
semilogx(f, 20*log10(Y_1W), '-b');




