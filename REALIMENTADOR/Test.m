clear all 
close all
clc
set(0,'defaulttextinterpreter','latex');
set(0,'defaultlegendinterpreter','latex');

s = tf('s');
w = logspace(0,10,50E3);

GH1 = 1/(1+s/(2*pi*4E6));

[mag_GH1, fase_GH1] = bode(GH1, w);
[Gm, Pm, Wcg, Wcp] = margin(mag_GH1, fase_GH1, w);
freq_pm1 = Wcp/2/pi;

figure
yyaxis left
semilogx(w/2/pi, mag2db(squeeze(mag_GH1)));
ylabel("Amplitud (dB)")
hold on

yyaxis right
semilogx(w/2/pi, squeeze(fase_GH1));
plot([400E3 400E3], [-90 0])
grid minor
xlabel("Frecuencia [Hz]")
xlim([1 1E10])