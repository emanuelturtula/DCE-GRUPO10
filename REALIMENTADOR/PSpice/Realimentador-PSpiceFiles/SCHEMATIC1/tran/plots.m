clear all;
close all;
%% Lazo cerrado
cl_1k = csvread('1k_cl.csv',1,0);
cl_2k = csvread('2k_cl.csv',1,0);
cl_3k = csvread('3k_cl.csv',1,0);
cl_5k = csvread('5k_cl.csv',1,0);
cl_7k = csvread('7k_cl.csv',1,0);
cl_10k = csvread('10k_cl.csv',1,0);
cl_15k = csvread('15k_cl.csv',1,0);
cl_20k = csvread('20k_cl.csv',1,0);
cl_50k = csvread('50k_cl.csv',1,0);
cl_100k = csvread('100k_cl.csv',1,0);


cl = {cl_1k,cl_2k,cl_3k,cl_5k,cl_7k,cl_10k,cl_15k,cl_20k,cl_50k,cl_100k};

F = [1 2 3 5 7 10 15 20 50 100];
A = [];
for i = 1:length(cl)
    A = [A (max(cl{i}(:,2))-min(cl{i}(:,2)))/2];
end

A = A/1.1;
A_db = mag2db(A);

fig=figure;
plot(F,A);
grid on;
xlabel('Frecuencia [kHz]');
ylabel('Amplitud [dB]');
print('Rta_cl.png','-dpng')

Dist = [1.9032E+00, 1.3251E-01, 2.2351E-01, 9.7933E-02, 1.0521E-01, 9.5484E-02, 2.1812E-01, 1.4366E-01];

fig=figure;
plot(F(1:end-2),Dist);
grid on;
xlabel('Frecuencia [kHz]');
ylabel('Distorsión (%)');
print('Dist_cl.png','-dpng')

%% Lazo abierto
clear A;

ol_1k = csvread('1k_ol.csv',1,0);
ol_2k = csvread('2k_ol.csv',1,0);
ol_3k = csvread('3k_ol.csv',1,0);
ol_5k = csvread('5k_ol.csv',1,0);
ol_7k = csvread('7k_ol.csv',1,0);
ol_10k = csvread('10k_ol.csv',1,0);
ol_15k = csvread('15k_ol.csv',1,0);
ol_20k = csvread('20k_ol.csv',1,0);
ol_50k = csvread('50k_ol.csv',1,0);
ol_100k = csvread('100k_ol.csv',1,0);


ol = {ol_1k,ol_2k,ol_3k,ol_5k,ol_7k,ol_10k,ol_15k,ol_20k,ol_50k,ol_100k};

A = [];
for i = 1:length(ol)
    A = [A (max(ol{i}(:,2))-min(ol{i}(:,2)))/2];
end

A = A/1.1;
A_db = mag2db(A);

fig=figure;
plot(F,A);
grid on;
xlabel('Frecuencia [kHz]');
ylabel('Amplitud [dB]');
print('Rta_ol.png','-dpng')

Dist = [3.0777e0, 3.2019e0, 1.0676e00, 9.3421e-01, 6.6152e-01, 9.6653e-01, 4.3759e-01, 1.1885e00 ];

fig=figure;
plot(F(1:end-2),Dist);
grid on;
xlabel('Frecuencia [kHz]');
ylabel('Distorsión (%)');
print('Dist_ol.png','-dpng')
