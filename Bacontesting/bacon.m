clc
clear all

% SKOGSBACON
t = [25:5:60];
sb_vekt_for = [12.65 12.51 14.25 11.87 11.55 12.02 13.17 12.76];
sb_vekt_etter = [6.46 6.3 7.38 5.43 4.61 4.85 5.31 4.85 ];

sprohet = [0.6 0.6 0.5 0.6 0.8 0.8 0.9 0.85];

tpapir_for = [2.53 2.03 2.60 2.43 2.51 2.57 2.37 2.65];
tpapir_etter = [5.05 4.80 5.05 5.07 5.2 5.6 4.8 5.45];
diff_sbvekt = sb_vekt_for-sb_vekt_etter;

diff_tpapir = tpapir_etter-tpapir_for;
N = length(diff_tpapir);

for i=1:1:(N/2)
subplot(2,2,1)
plot(t(i),diff_sbvekt(2*i-1),'.',t(i),diff_sbvekt(2*i),'.',t(i),mean([diff_sbvekt(2*i-1) diff_sbvekt(2*i)]),'d');
hold on
plot([t(i) t(i)],[diff_sbvekt(2*i-1) diff_sbvekt(2*i)]);
end

title('Massetap bacon');
xlabel('Tid (s)')
ylabel('Massetap (g)')

for i=1:1:(N/2)
subplot(2,2,2)
plot(t(i),diff_tpapir(2*i-1),'.',t(i),diff_tpapir(2*i),'.',t(i),mean([diff_tpapir(2*i-1) diff_tpapir(2*i)]),'d');
hold on
plot([t(i) t(i)],[diff_tpapir(2*i-1) diff_tpapir(2*i)]);
end

title('Fettabsorbering tørkepapir');
xlabel('Tid (s)')
ylabel('Masse (g)')

vanntap = diff_sbvekt - diff_tpapir;

for i=1:1:(N/2)
subplot(2,2,3)
plot(t(i),vanntap(2*i-1),'.',t(i),vanntap(2*i),'.',t(i),mean([vanntap(2*i-1) vanntap(2*i)]),'d');
hold on
plot([t(i) t(i)],[vanntap(2*i-1) vanntap(2*i)]);
end

title('Vanntap');
xlabel('Tid (s)')
ylabel('Masse (g)')

volum_for = [95*34*3.4 98*36*3.3 103.5*42.4*2.7 99*39*3.4 98*35*3.8 89*38*3.5 97*40*3.4 94*45*3.3];
volum_etter = [60*31*2 62*30*2.2 59*34*2.4 59*34*2.4 59*26*2 62*24*2.2 60*75*2.5 67*34*2.2];
volumtap = volum_for - volum_etter;

for i=1:1:(N/2)
subplot(2,2,4)
plot(t(i),volumtap(2*i-1),'.',t(i),volumtap(2*i),'.',t(i),mean([volumtap(2*i-1) volumtap(2*i)]),'d');
hold on
plot([t(i) t(i)],[volumtap(2*i-1) volumtap(2*i)]);
end

title('Volumtap');
xlabel('Tid (s)')
ylabel('Volum (mm^3)')

figure(2)
for i=1:1:(N/2)
plot(t(i),sprohet(2*i-1),'x',t(i),sprohet(2*i),'x');
hold on
end