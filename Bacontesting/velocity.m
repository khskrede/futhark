t = [25:5:60];

tpapir_for = [2.53 2.03 2.60 2.43 2.51 2.57 2.37 2.65];
tpapir_etter = [5.05 4.80 5.05 5.07 5.2 5.6 4.8 5.45];
diff_tpapir = tpapir_etter-tpapir_for;

areal_for = [95*34 98*36 103.5*42.4 99*39 98*35 89*38 97*40 94*45];
areal_etter = [60*31 62*30 59*34 59*34 59*26 62*24 60*75 67*34];
areal_avg = 0.5*(areal_for + areal_etter);

rho = 0.9;
diff_vol = rho*diff_tpapir;
flow = diff_vol./t;

% I m/s
velocity = flow./areal_avg

