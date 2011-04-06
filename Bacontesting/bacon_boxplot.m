clc
clear all


% CREATING TIME ARRAY

% SKOGSBACON
t = length(2*(70-25)/5);

j=0;
for i=25:5:70
    j = j+1;
    t(2*j) = i;
    t(2*j-1) = i;
end


% DATA SHEET FOR 'SKOGSBACON', FROM TIME t = 25s -> 40s
sb_vekt_for = [12.65 12.51 14.25 11.87 11.55 12.02 13.17 12.76];
sb_vekt_etter = [6.46 6.3 7.38 5.43 4.61 4.85 5.31 4.85 ];

tpapir_for = [2.53 2.03 2.60 2.43 2.51 2.57 2.37 2.65];
tpapir_etter = [5.05 4.80 5.05 5.07 5.2 5.6 4.8 5.45];

diff_sbvekt = sb_vekt_for-sb_vekt_etter;
diff_tpapir = tpapir_etter-tpapir_for;

volum_for = [95*34*3.4 98*36*3.3 103.5*42.4*2.7 99*39*3.4 98*35*3.8 89*38*3.5 97*40*3.4 94*45*3.3];
volum_etter = [60*31*2 62*30*2.2 59*34*2.4 59*34*2.4 59*26*2 62*24*2.2 60*75*2.5 67*34*2.2];
volumtap = volum_for - volum_etter; 

sprohet = [0.6 0.6 0.5 0.6 0.8 0.8 0.9 0.85];
N = length(diff_tpapir);

% DATA SHEET FOR REGULAR BACON, FROM TIME t = 40s -> 70s
fb_vekt_for = [17.23 15.7 20.23 17.64 17.8 31.73 18.6 18.6 17 16.06 15.81 15.73 16.21 16.5];
fb_vekt_etter = [6.14 5.32 8.46 7.2 6.35 15.89 5.95 5.95 4.57 4.2 4.29 4.16 4.35 4.67];

fb_tpapir_for = [2.64 2.60 2.95 2.88 4.64 4.79 3.8 3.8 2.8 2.73 2.43 2.49 2.52 2.51];
fb_tpapir_etter = [5.95 5.77 6.33 6.58 8.57 10.23 7.74 7.74 7.28 7.1 6.72 6.81 6.99 6.45];
diff_fb_tpapir = fb_tpapir_etter-fb_tpapir_for;

fb_volum_for = [192*43*1.9 175*47*2.8 220*50*2.4 198*51*2.6 267*33*2.5 275*42*3.6 212*40*1.8 212*40*1.8 225*41*1 214*41.2*1.0 215*40.5*1.0 214*41*1 222*40*1 212*40*1];
fb_volum_etter = [88*33*1.6 86*34*1.5 220*50*2.4 95*44*1.5 132*25*2.0 157*33.6*2.2 96*34*1.0 96*34*1.0 95*23.5*1.0 214*41.2*1.0 105*26*1 103*22*1 222*40*1 80*22*1];

fb_sprohet = [0.7 0.7 0.8 0.8 0.9 0.7 0.95 0.95 0.98 0.98 1.0 1.0 1.1 1.1];

M = length(fb_vekt_for);

for k=1:(N/2)
    sb_relative(2*k-1) = (sb_vekt_for(2*k-1) - sb_vekt_etter(2*k-1))/(sb_vekt_for(2*k-1));
    sb_relative(2*k) = (sb_vekt_for(2*k) - sb_vekt_etter(2*k))/(sb_vekt_for(2*k));
end
for m=1:(M/2)
    fb_relative(2*m-1) = (fb_vekt_for(2*m-1)-fb_vekt_etter(2*m-1))/(fb_vekt_for(2*m-1));
    fb_relative(2*m) = (fb_vekt_for(2*m)-fb_vekt_etter(2*m))/(fb_vekt_for(2*m));
end

figure(1) % Plot for thick bacon
title('Thick bacon')
subplot(2,2,1)
boxplot(diff_sbvekt,t(1:N))
ylabel('Weight (g)','InterPreter','LaTeX')
xlabel('Time (s)','InterPreter','LaTeX')
title('Weight loss','InterPreter','LaTeX')

subplot(2,2,2)
boxplot(sb_relative,t(1:N))
ylabel('Weight loss ($\%$)','InterPreter','LaTeX')
xlabel('Time (s)','InterPreter','LaTeX')
title('Relative weight','InterPreter','LaTeX')

subplot(2,2,3)
boxplot(diff_sbvekt-diff_tpapir,t(1:N))
title('Water loss','InterPreter','LaTeX')
ylabel('Weight (g)','InterPreter','LaTeX')
xlabel('Time (s)','InterPreter','LaTeX')

subplot(2,2,4)
boxplot(diff_tpapir,t(1:N))
title('Fat loss','InterPreter','LaTeX')
ylabel('Weight (g)','InterPreter','LaTeX')
xlabel('Time (s)','InterPreter','LaTeX')

figure(2) % Plot for regular bacon
title('Regular bacon')

subplot(2,2,1)
boxplot(fb_vekt_for-fb_vekt_etter,t(7:end))
title('Weight loss','InterPreter','LaTeX')
ylabel('Weight (g)','InterPreter','LaTeX')
xlabel('Time (s)','InterPreter','LaTeX')

subplot(2,2,2)
boxplot(fb_relative,t(7:end))
title('Relative weight loss','InterPreter','LaTeX')
ylabel('Weight loss ($\%$)','InterPreter','LaTeX')
xlabel('Time (s)','InterPreter','LaTeX')

subplot(2,2,3)
boxplot(diff_fb_tpapir,t(7:end))
title('Fat loss','InterPreter','LaTeX')
xlabel('Time (s)','InterPreter','LaTeX')
ylabel('Weight (g)','InterPreter','LaTeX')

subplot(2,2,4)
boxplot(diff_fb_tpapir-fb_vekt_for-fb_vekt_etter,t(7:end))
title('Water loss','InterPreter','LaTeX')
xlabel('Time (s)','InterPreter','LaTeX')
ylabel('Weight (g)','InterPreter','LaTeX')

sb_vekt_for = [12.65 12.51 14.25 11.87 11.55 12.02 13.17 12.76];
sb_vekt_etter = [6.46 6.3 7.38 5.43 4.61 4.85 5.31 4.85 ];

% Starter på 40s for FirstPrice-bacon
% fb_vekt_for = [17.23 15.7 20.23 17.64 17.8 31.73 18.6 18.6 17 16.06 15.81 15.73 16.21 16.5]; 
% fb_vekt_etter = [6.14 5.32 8.46 7.2 6.35 15.89 5.95 5.95 4.57 4.2 4.29 4.16 4.35 4.67];
% 
% sprohet = [0.6 0.6 0.5 0.6 0.8 0.8 0.9 0.85];
% 
% tpapir_for = [2.53 2.03 2.60 2.43 2.51 2.57 2.37 2.65];
% tpapir_etter = [5.05 4.80 5.05 5.07 5.2 5.6 4.8 5.45];
% diff_sbvekt = sb_vekt_for-sb_vekt_etter;
% 
% diff_tpapir = tpapir_etter-tpapir_for;
% N = length(diff_tpapir);

% subplot(2,2,1)
% boxplot(diff_sbvekt,t(1:N))
% subplot(2,2,2)
% boxplot(diff_tpapir,t(1:N))
% subplot(2,2,3)
% boxplot(diff_sbvekt-diff_tpapir,t(1:N))
% subplot(2,2,4)
% boxplot(fb_vekt_for-fb_vekt_etter,t(7:end))
