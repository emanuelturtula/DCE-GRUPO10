close all
clear all
clc
set(0,'defaulttextinterpreter','latex');
set(0,'defaultlegendinterpreter','latex');

%% Distorsión a 1k

sim_1k = dlmread(fullfile("30W","Distorsion_1k.txt"), "\t", 1, 0);     %leo el archivo
fft_1k = dlmread(fullfile("30W","FFT_1k.txt"), "\t", 1, 0);

time = sim_1k(:,1)*1E3;                                    %separo las columnas
vo = sim_1k(:,2);
vin = sim_1k(:,3);
f = fft_1k(:,1);
VO_db = fft_1k(:,2);
VO = db2mag(VO_db); 


%genero un vector de tiempo
vo_rms = rms(vo);
vin_rms = rms(vin);

figure
plot(time, vo, '-k', 'DisplayName', 'Salida');
hold on
plot(time, vin, '-r', 'DisplayName', 'Entrada');
plot([time(1) time(end)], [vo_rms vo_rms], '-- b', 'linewidth', 2, ...
    'DisplayName', sprintf('Salida (RMS) = %.2f V', vo_rms));
plot([time(1) time(end)], [vin_rms vin_rms], '-- g', 'linewidth', 2, ...
    'DisplayName', sprintf('Entrada (RMS) = %.2f V', vin_rms));

title("\begin{tabular}{c} Lazo abierto \\ Simulaci\'on a 1 kHz @ 30 W \end{tabular}", 'Fontsize', 14)
legend('location', 'Southeast');
xlabel("Tiempo [ms]");
ylabel("Tensi\'on [V]");
xlim([0 4])
grid minor

print("Dist_ol_1k_30W.eps", '-depsc', '-r0', '-tiff');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFT

figure
semilogx(f/1E3, VO_db,'DisplayName', "FFT");
hold on
indexes(1) = find(f == 1E3);
plot(f(indexes(1))/1E3, VO_db(indexes(1)), '*', 'Markersize', 4, 'DisplayName', ...
    sprintf('Fundamental = %.2f V', VO(indexes(1))));

dist1k = 0;
for i=2:9
    indexes(i) = find(f == i*1E3);
    plot(f(indexes(i))/1E3, VO_db(indexes(i)), '*', 'Markersize', 4, 'DisplayName', ...
        sprintf('Armonico %i = %.2f mV', i-1, VO(indexes(i))*1E3));
    dist1k = dist1k + VO(indexes(i))^2;
end

dist1k = sqrt(dist1k)/VO(indexes(1));

title("\begin{tabular}{c} Lazo abierto \\ FFT de la simulaci\'on a 1 kHz @ 30 W \end{tabular}", 'Fontsize', 14)
ylabel('Amplitud [dB]')
xlabel('Frecuencia [kHz]')
xlim([0.5 10E3])
hLeg=legend("location", 'Northeast');
hLeg.NumColumns=2;
grid minor

print("Dist_ol_1k_FFT_30W.eps", '-depsc', '-r0', '-tiff');

%% Distorsión a 10k

sim_10k = dlmread(fullfile("30W","Distorsion_10k.txt"), "\t", 1, 0);     %leo el archivo
fft_10k = dlmread(fullfile("30W","FFT_10k.txt"), "\t", 1, 0);

time = sim_10k(:,1)*1E3;                                    %separo las columnas
vo = sim_10k(:,2);
vin = sim_10k(:,3);
f = fft_10k(:,1);
VO_db = fft_10k(:,2);
VO = db2mag(VO_db); 


%genero un vector de tiempo
vo_rms = rms(vo);
vin_rms = rms(vin);

figure
plot(time, vo, '-k', 'DisplayName', 'Salida');
hold on
plot(time, vin, '-r', 'DisplayName', 'Entrada');
plot([time(1) time(end)], [vo_rms vo_rms], '-- b', 'linewidth', 2, ...
    'DisplayName', sprintf('Salida (RMS) = %.2f V', vo_rms));
plot([time(1) time(end)], [vin_rms vin_rms], '-- g', 'linewidth', 2, ...
    'DisplayName', sprintf('Entrada (RMS) = %.2f V', vin_rms));

title("\begin{tabular}{c} Lazo abierto \\ Simulaci\'on a 10 kHz @ 30 W \end{tabular}", 'Fontsize', 14)
legend('location', 'Southeast');
xlabel("Tiempo [ms]");
ylabel("Tensi\'on [V]");
xlim([0 1])
grid minor

print("Dist_ol_10k_30W.eps", '-depsc', '-r0', '-tiff');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFT

figure
semilogx(f/1E3, VO_db,'DisplayName', "FFT");
hold on
indexes(1) = find(f == 10E3);
plot(f(indexes(1))/1E3, VO_db(indexes(1)), '*', 'Markersize', 4, 'DisplayName', ...
    sprintf('Fundamental = %.2f V', VO(indexes(1))));

dist10k = 0;
for i=2:9
    indexes(i) = find(f == i*10E3);
    plot(f(indexes(i))/1E3, VO_db(indexes(i)), '*', 'Markersize', 4, 'DisplayName', ...
        sprintf('Armonico %i = %.2f mV', i-1, VO(indexes(i))*1E3));
    dist10k = dist10k + VO(indexes(i))^2;
end

dist10k = sqrt(dist10k)/VO(indexes(1));

title("\begin{tabular}{c} Lazo abierto \\ FFT de la simulaci\'on a 10 kHz @ 30 W \end{tabular}", 'Fontsize', 14)
ylabel('Amplitud [dB]')
xlabel('Frecuencia [kHz]')
xlim([0.5 100E3])
hLeg=legend("location", 'Northeast');
hLeg.NumColumns=1;
grid minor

print("Dist_ol_10k_FFT_30W.eps", '-depsc', '-r0', '-tiff');

%% Distorsión a 20k

sim_20k = dlmread(fullfile("30W","Distorsion_20k.txt"), "\t", 1, 0);     %leo el archivo
fft_20k = dlmread(fullfile("30W","FFT_20k.txt"), "\t", 1, 0);

time = sim_20k(:,1)*1E3;                                    %separo las columnas
vo = sim_20k(:,2);
vin = sim_20k(:,3);
f = fft_20k(:,1);
VO_db = fft_20k(:,2);
VO = db2mag(VO_db); 


%genero un vector de tiempo
vo_rms = rms(vo);
vin_rms = rms(vin);

figure
plot(time, vo, '-k', 'DisplayName', 'Salida');
hold on
plot(time, vin, '-r', 'DisplayName', 'Entrada');
plot([time(1) time(end)], [vo_rms vo_rms], '-- b', 'linewidth', 2, ...
    'DisplayName', sprintf('Salida (RMS) = %.2f V', vo_rms));
plot([time(1) time(end)], [vin_rms vin_rms], '-- g', 'linewidth', 2, ...
    'DisplayName', sprintf('Entrada (RMS) = %.2f V', vin_rms));

title("\begin{tabular}{c} Lazo abierto \\ Simulaci\'on a 20 kHz @ 30 W \end{tabular}", 'Fontsize', 14)
legend('location', 'Southeast');
xlabel("Tiempo [ms]");
ylabel("Tensi\'on [V]");
xlim([0 0.5])
grid minor

print("Dist_ol_20k_30W.eps", '-depsc', '-r0', '-tiff');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFT

figure
semilogx(f/1E3, VO_db,'DisplayName', "FFT");
hold on
indexes(1) = find(f == 20E3);
plot(f(indexes(1))/1E3, VO_db(indexes(1)), '*', 'Markersize', 4, 'DisplayName', ...
    sprintf('Fundamental = %.2f V', VO(indexes(1))));

dist20k = 0;
for i=2:4
    indexes(i) = find(f == i*20E3);
    plot(f(indexes(i))/1E3, VO_db(indexes(i)), '*', 'Markersize', 4, 'DisplayName', ...
        sprintf('Armonico %i = %.2f mV', i-1, VO(indexes(i))*1E3));
    dist20k = dist20k + VO(indexes(i))^2;
end

dist20k = sqrt(dist20k)/VO(indexes(1));

title("\begin{tabular}{c} Lazo abierto \\ FFT de la simulaci\'on a 20 kHz @ 30 W \end{tabular}", 'Fontsize', 14)
ylabel('Amplitud [dB]')
xlabel('Frecuencia [kHz]')
xlim([1 100E3])
hLeg=legend("location", 'Northeast');
hLeg.NumColumns=1;
grid minor

print("Dist_ol_20k_FFT_30W.eps", '-depsc', '-r0', '-tiff');

%% Distorsion 1k 1W
clear all
close all

sim_1k = dlmread(fullfile("1W","Distorsion_1k_1W.txt"), "\t", 1, 0);     %leo el archivo
fft_1k = dlmread(fullfile("1W","FFT_1k_1W.txt"), "\t", 1, 0);

time = sim_1k(:,1)*1E3;                                    %separo las columnas
vo = sim_1k(:,2);
vin = sim_1k(:,3);
f = fft_1k(:,1);
VO_db = fft_1k(:,2);
VO = db2mag(VO_db); 


%genero un vector de tiempo
vo_rms = rms(vo);
vin_rms = rms(vin);

figure
plot(time, vo, '-k', 'DisplayName', 'Salida');
hold on
plot(time, vin, '-r', 'DisplayName', 'Entrada');
plot([time(1) time(end)], [vo_rms vo_rms], '-- b', 'linewidth', 2, ...
    'DisplayName', sprintf('Salida (RMS) = %.2f V', vo_rms));
plot([time(1) time(end)], [vin_rms vin_rms], '-- g', 'linewidth', 2, ...
    'DisplayName', sprintf('Entrada (RMS) = %.2f V', vin_rms));

title("\begin{tabular}{c} Lazo abierto \\ Simulaci\'on a 1 kHz @ 1 W \end{tabular}", 'Fontsize', 14)
legend('location', 'Southeast');
xlabel("Tiempo [ms]");
ylabel("Tensi\'on [V]");
xlim([0 4])
grid minor

print("Dist_ol_1k_1W.eps", '-depsc', '-r0', '-tiff');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFT

figure
semilogx(f/1E3, VO_db,'DisplayName', "FFT");
hold on
indexes(1) = find(f == 1E3);
plot(f(indexes(1))/1E3, VO_db(indexes(1)), '*', 'Markersize', 4, 'DisplayName', ...
    sprintf('Fundamental = %.2f V', VO(indexes(1))));

dist1k = 0;
for i=2:9
    indexes(i) = find(f == i*1E3);
    plot(f(indexes(i))/1E3, VO_db(indexes(i)), '*', 'Markersize', 4, 'DisplayName', ...
        sprintf('Armonico %i = %.2f mV', i-1, VO(indexes(i))*1E3));
    dist1k = dist1k + VO(indexes(i))^2;
end

dist1k = sqrt(dist1k)/VO(indexes(1));

title("\begin{tabular}{c} Lazo abierto \\ FFT de la simulaci\'on a 1 kHz @ 1 W \end{tabular}", 'Fontsize', 14)
ylabel('Amplitud [dB]')
xlabel('Frecuencia [kHz]')
xlim([0.5 10E3])
hLeg=legend("location", 'Northeast');
hLeg.NumColumns=2;
grid minor

print("Dist_ol_1k_FFT_1W.eps", '-depsc', '-r0', '-tiff');