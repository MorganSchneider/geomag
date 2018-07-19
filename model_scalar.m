%
load('./EEJ_Data/Swarm_Scalar.mat')
%%
rad = pi/180; %radians

t_scalar = [scalar(1).time, scalar(2).time];
yearEEJ = decimalYear(t_EEJ);
yearScalar = [scalar(1).year, scalar(2).year];

timeStart = 2015.0;
timeEnd = 2015.25;

inds = find(yearScalar >= timeStart & yearScalar < timeEnd);
t = t_scalar(inds);
mjd = datenum(datetime(t,'ConvertFrom','posixtime')) - datenum(2000,1,1,0,0,0);

r = [scalar(1).rad, scalar(2).rad]; r = r(inds);
theta = (90 - [scalar(1).geolat, scalar(2).geolat]) * rad; theta = theta(inds);
phi = [scalar(1).lon, scalar(2).lon] * rad; phi = phi(inds);
N = 10;

NSH = (N+1)^2-1;

F_synth = rot90([scalar(1).F, scalar(2).F], 3); F_synth = F_synth(inds);

%%

g_init = ones(NSH,1);

for x = 1:50
    gprev = g_init;
    beta = find_beta(F_synth, r, theta, phi, g_init, N);
    [~,J] = find_J(r, theta, phi, g_init, N);
    delta = J\beta;
    g_init = g_init + gamma*delta;
    error = norm(g_init - gprev) / norm(g_init)
end
%final_error = norm(g_synth - g_init) / norm(g_synth)

%%

B_scalar_model = find_B(r, theta, phi, g_init, N);
F_scalar_model = find_F(r, theta, phi, g_init, N);

% compare to CHAOS
B_chaos = synth_values(r, theta, phi, pp_N, mjd);
F_chaos = find_F(r, theta, phi, pp_N, mjd);



%%

Br = B_scalar_model(:,1);
Bt = B_scalar_model(:,2);
Bp = B_scalar_model(:,3);
deltaBr = B_scalar_model(:,1) - B_chaos(:,1);
deltaBt = B_scalar_model(:,2) - B_chaos(:,2);
deltaBp = B_scalar_model(:,3) - B_chaos(:,3);
deltaF = F_scalar_model - F_chaos;
lat_deg = 90 - (theta_s./rad);
lon_deg = phi_s./rad;
lon_deg_r = phi_r./rad;

figure(1)
subplot(2,1,1)
plot(lon_deg, B_scalar_model(:,1) - B_model(:,1), '*')
xlabel('Longitude (deg)')
ylabel('\alpha (nT)')
title('B_{r} Model Residuals')
subplot(2,1,2)
plot(qd_s, F_scalar_model - F_model, '*')
xlabel('Quasi-Dipole Latitude (deg)')
ylabel('\beta (nT)')
title('F Model Residuals')

lat_rng = [-90 90];
lon_rng = [-180 180];

lat_lin = linspace(min(lat_deg), max(lat_deg), 2000);
lon_lin = linspace(min(lon_deg), max(lon_deg), 2000);
[lat_grid, lon_grid] = meshgrid(lat_lin, lon_lin);
Br_grid = griddata(lat_deg, lon_deg, Br, lat_grid, lon_grid, 'cubic');
Bt_grid = griddata(lat_deg, lon_deg, Bt, lat_grid, lon_grid, 'cubic');
Bp_grid = griddata(lat_deg, lon_deg, Bp, lat_grid, lon_grid, 'cubic');
F_grid = griddata(lat_deg, lon_deg, F_scalar_model, lat_grid, lon_grid, 'cubic');
dBr_grid = griddata(lat_deg, lon_deg, deltaBr, lat_grid, lon_grid, 'cubic');
dBt_grid = griddata(lat_deg, lon_deg, deltaBt, lat_grid, lon_grid, 'cubic');
dBp_grid = griddata(lat_deg, lon_deg, deltaBp, lat_grid, lon_grid, 'cubic');
dF_grid = griddata(lat_deg, lon_deg, deltaF, lat_grid, lon_grid, 'cubic');

% why
figure(2)
load coastlines

subplot(2,2,1)

ax = worldmap(lat_rng, lon_rng);
contourfm(lat_grid, lon_grid, Br_grid)
plotm(coastlat, coastlon, 'Color', 'black')
title('B_r (nT)')
contourcbar('southoutside')

subplot(2,2,2)

worldmap(lat_rng, lon_rng);
contourfm(lat_grid, lon_grid, Bt_grid)
plotm(coastlat, coastlon, 'Color', 'black')
title('B_{\theta} (nT)')
contourcbar('southoutside')

subplot(2,2,3)

worldmap(lat_rng, lon_rng);
contourfm(lat_grid, lon_grid, Bp_grid)
plotm(coastlat, coastlon, 'Color', 'black')
title('B_{\phi} (nT)')
contourcbar('southoutside')

subplot(2,2,4)

worldmap(lat_rng, lon_rng);
contourfm(lat_grid, lon_grid, F_grid)
plotm(coastlat, coastlon, 'Color', 'black')
title('F (nT)')
contourcbar('southoutside')


figure(3)

subplot(2,2,1)

worldmap(lat_rng, lon_rng);
contourfm(lat_grid, lon_grid, dBr_grid)
plotm(coastlat, coastlon, 'Color', 'black')
title('\delta B_r (nT)')
contourcbar('southoutside')

subplot(2,2,2)

worldmap(lat_rng, lon_rng);
contourfm(lat_grid, lon_grid, dBt_grid)
plotm(coastlat, coastlon, 'Color', 'black')
title('\delta B_{\theta} (nT)')
contourcbar('southoutside')

subplot(2,2,3)

worldmap(lat_rng, lon_rng);
contourfm(lat_grid, lon_grid, dBp_grid)
plotm(coastlat, coastlon, 'Color', 'black')
title('\delta B_{\phi} (nT)')
contourcbar('southoutside')

subplot(2,2,4)

worldmap(lat_rng, lon_rng);
contourfm(lat_grid, lon_grid, dF_grid)
plotm(coastlat, coastlon, 'Color', 'black')
title('\delta F (nT)')
contourcbar('southoutside')




