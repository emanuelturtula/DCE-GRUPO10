close all
clear all
clc
set(0,'defaulttextinterpreter','latex');
set(0,'defaultlegendinterpreter','latex');

%%

eficiencia = dlmread("Eficiencia_ol.txt", '\t', 1, 0);

y = eficiencia(:,2);
x = eficiencia(:,3);

figure
plot(x, y, 'DisplayName', 'Curva de eficiencia');

title("\begin{tabular}{c} Lazo abierto \\ Simulaci\'on de eficiencia \end{tabular}", 'Fontsize', 14)
legend('location', 'Southeast');
xlabel("Potencia [W]");
ylabel("Eficiencia [\%]");
xlim([0 50])
ylim([55 100])
grid minor

print("Eficiencia_ol.eps", '-depsc', '-r0', '-tiff');