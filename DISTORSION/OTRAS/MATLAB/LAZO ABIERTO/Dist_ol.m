close all
clear all
clc
set(0,'defaulttextinterpreter','latex');
set(0,'defaultlegendinterpreter','latex');

%% Distorsión a 1k (1W)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tiempo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim_1k = dlmread(fullfile("signals","Distorsion_1k_1W_ol.txt"), "\t", 1, 0);  
time = sim_1k(:,1)*1E3; %Tiempo (en ms)                                                                  
vo = sim_1k(:,2);       %Salida
vo_rms = rms(vo);       %Valor RMS de la salida
vo_avg = mean(vo);      %Offset de la salida

figure
plot(time, vo, '-k', 'DisplayName', 'Salida');
hold on
plot([time(1) time(end)], [vo_rms vo_rms], '--b', 'linewidth', 2, ...
    'DisplayName', sprintf('Salida (RMS) = %.2f V', vo_rms));
plot([time(1) time(end)], [vo_avg vo_avg], '--r', 'linewidth', 2, ...
    'DisplayName', sprintf('DC Offset = %.2f V', vo_avg));

title("\begin{tabular}{c} Lazo abierto \\ Simulaci\'on a 1 kHz @ 1 W \end{tabular}", 'Fontsize', 14)
legend('location', 'Southeast');
xlabel("Tiempo [ms]");
ylabel("Tensi\'on [V]");
xlim([0 10])
grid minor

print("fig/Dist_1k_1W_ol.eps", '-depsc', '-r0', '-tiff');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fft_1k = dlmread(fullfile("signals","FFT_1k_1W_ol.txt"), "|", 1, 0);
f = fft_1k(:,1);                %Frecuencia
VO_db = fft_1k(:,2);            %Amplitud (dB RMS)
VO = db2mag(VO_db)*sqrt(2);     %Amplitud (lineal pico)
fundamental = 1E3;              %Frecuencia fundamental
N_arm = 9;                      %Cantidad de armónicos

indexes = find_freq_idx(f, fundamental, N_arm); %Busco los indices de las 
                                                %frecuencias
                                                
dist_1k_1W = rssq(VO(indexes(2:end)))/VO(indexes(1));   %Calculo la 
                                                        %distorsión

figure
semilogx(f/1E3, VO_db,'DisplayName', "FFT");
hold on
plot(f(indexes(1))/1E3, VO_db(indexes(1)), '*', 'Markersize', 4, 'DisplayName', ...
    sprintf('Fundamental = %.2f V', VO(indexes(1))));
for i=2:N_arm
    plot(f(indexes(i))/1E3, VO_db(indexes(i)), '*', 'Markersize', 4, 'DisplayName', ...
        sprintf('Armonico %i = %.2f mV', i-1, VO(indexes(i))*1E3));
end

txt1 = "\begin{tabular}{c} Lazo abierto \\ FFT de la simulaci\'on a 1 kHz @ 1 W \end{tabular}";
txt2 = sprintf("Distorsi\\'on = %.3f \\%%", dist_1k_1W*100);
title({txt1, txt2}, 'Fontsize', 12);
ylabel('Amplitud [$dB_{RMS}$]')
xlabel('Frecuencia [kHz]')
xlim([0.5 10E3])
hLeg=legend("location", 'Northeast');
hLeg.NumColumns=2;
grid minor

print("fig/FFT_1k_1W_ol.eps", '-depsc', '-r0', '-tiff');

%% Distorsión a 1k (30W)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tiempo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim_1k = dlmread(fullfile("signals","Distorsion_1k_30W_ol.txt"), "\t", 1, 0);  
time = sim_1k(:,1)*1E3; %Tiempo (en ms)                                                                  
vo = sim_1k(:,2);       %Salida
vo_rms = rms(vo);       %Valor RMS de la salida
vo_avg = mean(vo);      %Offset de la salida

figure
plot(time, vo, '-k', 'DisplayName', 'Salida');
hold on
plot([time(1) time(end)], [vo_rms vo_rms], '--b', 'linewidth', 2, ...
    'DisplayName', sprintf('Salida (RMS) = %.2f V', vo_rms));
plot([time(1) time(end)], [vo_avg vo_avg], '--r', 'linewidth', 2, ...
    'DisplayName', sprintf('DC Offset = %.2f V', vo_avg));

title("\begin{tabular}{c} Lazo abierto \\ Simulaci\'on a 1 kHz @ 30 W \end{tabular}", 'Fontsize', 14)
legend('location', 'Southeast');
xlabel("Tiempo [ms]");
ylabel("Tensi\'on [V]");
xlim([0 10])
grid minor

print("fig/Dist_1k_30W_ol.eps", '-depsc', '-r0', '-tiff');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fft_1k = dlmread(fullfile("signals","FFT_1k_30W_ol.txt"), "|", 1, 0);
f = fft_1k(:,1);                %Frecuencia
VO_db = fft_1k(:,2);            %Amplitud (dB RMS)
VO = db2mag(VO_db)*sqrt(2);     %Amplitud (lineal pico)
fundamental = 1E3;              %Frecuencia fundamental
N_arm = 9;                      %Cantidad de armónicos

indexes = find_freq_idx(f, fundamental, N_arm); %Busco los indices de las 
                                                %frecuencias
                                                
dist_1k_30W = rssq(VO(indexes(2:end)))/VO(indexes(1));   %Calculo la 
                                                        %distorsión

figure
semilogx(f/1E3, VO_db,'DisplayName', "FFT");
hold on
plot(f(indexes(1))/1E3, VO_db(indexes(1)), '*', 'Markersize', 4, 'DisplayName', ...
    sprintf('Fundamental = %.2f V', VO(indexes(1))));
for i=2:N_arm
    plot(f(indexes(i))/1E3, VO_db(indexes(i)), '*', 'Markersize', 4, 'DisplayName', ...
        sprintf('Armonico %i = %.2f mV', i-1, VO(indexes(i))*1E3));
end

txt1 = "\begin{tabular}{c} Lazo abierto \\ FFT de la simulaci\'on a 1 kHz @ 30 W \end{tabular}";
txt2 = sprintf("Distorsi\\'on = %.3f \\%%", dist_1k_30W*100);
title({txt1, txt2}, 'Fontsize', 12);
ylabel('Amplitud [$dB_{RMS}$]')
xlabel('Frecuencia [kHz]')
xlim([0.5 10E3])
hLeg=legend("location", 'Northeast');
hLeg.NumColumns=2;
grid minor

print("fig/FFT_1k_30W_ol.eps", '-depsc', '-r0', '-tiff');

%% Distorsión a 10k (30W)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tiempo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim_10k = dlmread(fullfile("signals","Distorsion_10k_30W_ol.txt"), "\t", 1, 0);  
time = sim_10k(:,1)*1E3; %Tiempo (en ms)                                                                  
vo = sim_10k(:,2);       %Salida
vo_rms = rms(vo);       %Valor RMS de la salida
vo_avg = mean(vo);      %Offset de la salida

figure
plot(time, vo, '-k', 'DisplayName', 'Salida');
hold on
plot([time(1) time(end)], [vo_rms vo_rms], '--b', 'linewidth', 2, ...
    'DisplayName', sprintf('Salida (RMS) = %.2f V', vo_rms));
plot([time(1) time(end)], [vo_avg vo_avg], '--r', 'linewidth', 2, ...
    'DisplayName', sprintf('DC Offset = %.2f V', vo_avg));

title("\begin{tabular}{c} Lazo abierto \\ Simulaci\'on a 10 kHz @ 30 W \end{tabular}", 'Fontsize', 14)
legend('location', 'Southeast');
xlabel("Tiempo [ms]");
ylabel("Tensi\'on [V]");
xlim([0 1])
grid minor

print("fig/Dist_10k_30W_ol.eps", '-depsc', '-r0', '-tiff');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fft_10k = dlmread(fullfile("signals","FFT_10k_30W_ol.txt"), "|", 1, 0);
f = fft_10k(:,1);                %Frecuencia
VO_db = fft_10k(:,2);            %Amplitud (dB RMS)
VO = db2mag(VO_db)*sqrt(2);     %Amplitud (lineal pico)
fundamental = 10E3;              %Frecuencia fundamental
N_arm = 6;                      %Cantidad de armónicos

indexes = find_freq_idx(f, fundamental, N_arm); %Busco los indices de las 
                                                %frecuencias
                                                
dist_10k_30W = rssq(VO(indexes(2:end)))/VO(indexes(1));   %Calculo la 
                                                        %distorsión

figure
semilogx(f/1E3, VO_db,'DisplayName', "FFT");
hold on
plot(f(indexes(1))/1E3, VO_db(indexes(1)), '*', 'Markersize', 4, 'DisplayName', ...
    sprintf('Fundamental = %.2f V', VO(indexes(1))));
for i=2:N_arm
    plot(f(indexes(i))/1E3, VO_db(indexes(i)), '*', 'Markersize', 4, 'DisplayName', ...
        sprintf('Armonico %i = %.2f mV', i-1, VO(indexes(i))*1E3));
end

txt1 = "\begin{tabular}{c} Lazo abierto \\ FFT de la simulaci\'on a 10 kHz @ 30 W \end{tabular}";
txt2 = sprintf("Distorsi\\'on = %.3f \\%%", dist_10k_30W*100);
title({txt1, txt2}, 'Fontsize', 12);
ylabel('Amplitud [$dB_{RMS}$]')
xlabel('Frecuencia [kHz]')
xlim([1 10E3])
hLeg=legend("location", 'Northeast');
hLeg.NumColumns=1;
grid minor

print("fig/FFT_10k_30W_ol.eps", '-depsc', '-r0', '-tiff');

%% Distorsión a 20k (30W)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tiempo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim_20k = dlmread(fullfile("signals","Distorsion_20k_30W_ol.txt"), "\t", 1, 0);  
time = sim_20k(:,1)*1E6; %Tiempo (en us)                                                                  
vo = sim_20k(:,2);       %Salida
vo_rms = rms(vo);       %Valor RMS de la salida
vo_avg = mean(vo);      %Offset de la salida

figure
plot(time, vo, '-k', 'DisplayName', 'Salida');
hold on
plot([time(1) time(end)], [vo_rms vo_rms], '--b', 'linewidth', 2, ...
    'DisplayName', sprintf('Salida (RMS) = %.2f V', vo_rms));
plot([time(1) time(end)], [vo_avg vo_avg], '--r', 'linewidth', 2, ...
    'DisplayName', sprintf('DC Offset = %.2f V', vo_avg));

title("\begin{tabular}{c} Lazo abierto \\ Simulaci\'on a 20 kHz @ 30 W \end{tabular}", 'Fontsize', 14)
legend('location', 'Southeast');
xlabel("Tiempo [$\mu s$]");
ylabel("Tensi\'on [V]");
xlim([0 500])
grid minor

print("fig/Dist_20k_30W_ol.eps", '-depsc', '-r0', '-tiff');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fft_20k = dlmread(fullfile("signals","FFT_20k_30W_ol.txt"), "|", 1, 0);
f = fft_20k(:,1);                %Frecuencia
VO_db = fft_20k(:,2);            %Amplitud (dB RMS)
VO = db2mag(VO_db)*sqrt(2);     %Amplitud (lineal pico)
fundamental = 20E3;              %Frecuencia fundamental
N_arm = 3;                      %Cantidad de armónicos

indexes = find_freq_idx(f, fundamental, N_arm); %Busco los indices de las 
                                                %frecuencias
                                                
dist_20k_30W = rssq(VO(indexes(2:end)))/VO(indexes(1));   %Calculo la 
                                                        %distorsión

figure
semilogx(f/1E3, VO_db,'DisplayName', "FFT");
hold on
plot(f(indexes(1))/1E3, VO_db(indexes(1)), '*', 'Markersize', 4, 'DisplayName', ...
    sprintf('Fundamental = %.2f V', VO(indexes(1))));
for i=2:N_arm
    plot(f(indexes(i))/1E3, VO_db(indexes(i)), '*', 'Markersize', 4, 'DisplayName', ...
        sprintf('Armonico %i = %.2f mV', i-1, VO(indexes(i))*1E3));
end

txt1 = "\begin{tabular}{c} Lazo abierto \\ FFT de la simulaci\'on a 20 kHz @ 30 W \end{tabular}";
txt2 = sprintf("Distorsi\\'on = %.3f \\%%", dist_20k_30W*100);
title({txt1, txt2}, 'Fontsize', 12);
ylabel('Amplitud [$dB_{RMS}$]')
xlabel('Frecuencia [kHz]')
xlim([2 10E3])
hLeg=legend("location", 'Northeast');
hLeg.NumColumns=1;
grid minor

print("fig/FFT_20k_30W_ol.eps", '-depsc', '-r0', '-tiff');

%% Funciones
function idx = find_freq_idx(freq,freq_fund, n)
    idx = 0;
    for i=1:n
        idx(i) = find(freq == i*freq_fund);
    end
end