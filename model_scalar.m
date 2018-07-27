%
load('./EEJ_Data/Swarm_Scalar.mat')
% load CHAOS model
load('./CHAOS-6_FWD/CHAOS-6-x5.mat')

%%
rad = pi/180; %radians

t_scalar = [scalar(1).time, scalar(2).time];
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

pp_N = pp;
pp_N.dim = N*(N+2);
coefs_tmp = reshape(pp.coefs, [], pp.pieces, pp.order);
pp_N.coefs = reshape(coefs_tmp(1:N*(N+2),:,:), [], pp.order);

qd = [scalar(1).qdlat, scalar(2).qdlat]; qd = qd(inds);

F_synth = rot90([scalar(1).F, scalar(2).F], 3); F_synth = F_synth(inds);

%%

% Build weighting function
w = zeros(length(r), length(r));
for i = 1:length(r)
    w(i,i) = sin(theta(i));
end

g_init = zeros(NSH,1);
g_init(index(1,0)) = -30000;

iter = 50;
gam = 0.5;

error = zeros(1,iter);
chisq = zeros(1,iter);
for x = 1:iter
    gprev = g_init;
    beta = find_beta(F_synth, r, theta, phi, g_init, N);
    [~,J] = find_J(r, theta, phi, g_init, N);
    delta = (J'*w*J)\(J'*w*beta);
    g_init = g_init + gam*delta;
    error = norm(g_init - gprev) / norm(g_init)
    chisq = norm(sqrt(w)*beta)
end


%%

F_scalar_model = find_F(r, theta, phi, g_init, N);
lon_deg = phi./rad;

close all

figure(1)
plot(qd, F_synth - F_scalar_model, '*')
xlabel('Quasi-Dipole Latitude (deg)')
ylabel('F Residuals (nT)')
title('F_{Swarm} - F_{Model}, 2015.0 - 2015.25')

%% Plot difference matrix

g_chaos = reindex(pp_N, N, mean(mjd));

dg = difference_mat(g_init, g_chaos);
Sm = sensitivity_mat(g_init, g_chaos);
Ps = powerspec(g_init);


figure(5)

subplot(1,3,1)
pcolor(-N:N, 1:N, dg)
caxis([-25 25])
set(gca, 'Ydir', 'reverse')
colormap(cool)
colorbar('h')
title('Difference Matrix')
xlabel('Spherical harmonic order')
ylabel('Spherical harmonic degree')

subplot(1,3,2)
pcolor(-N:N, 1:N, Sm)
caxis([-100 100])
set(gca, 'Ydir', 'reverse')
colormap(cool)
colorbar('h')
title('Sensitivity Matrix')
xlabel('Spherical harmonic order')
ylabel('Spherical harmonic degree')


subplot(1,3,3)
plot(1:N, Ps)
title('Power Spectrum')
xlabel('Spherical harmonic degree')
ylabel('Power spectrum')

%%

lat_rng = [-90 90];
lon_rng = [-180 180];

lat_lin = linspace(-90, 90, 200);
lon_lin = linspace(-180, 180, 400);
[lat_grid, lon_grid] = meshgrid(lat_lin, lon_lin);
lat_vector = 90 - reshape(lat_grid, [], 1);
lon_vector = reshape(lon_grid, [], 1);
r_const = ones(1, 80000) * 6371.2;
t_const = ones(1, 80000) * mean(mjd);

s = size(lat_grid);

B_scalar_model = find_B(r_const, lat_vector*rad, lon_vector*rad, g_init, N);
F_scalar_model = find_F(r_const, lat_vector*rad, lon_vector*rad, g_init, N);
B_chaos = synth_values(r_const, lat_vector, lon_vector, pp_N, t_const);
F_chaos = find_F(r_const, lat_vector, lon_vector, pp_N, t_const);
dB_scalar = B_scalar_model - B_chaos;
dF_scalar = F_scalar_model - F_chaos;

Br_grid = reshape(B_scalar_model(:,1), s);
Bt_grid = reshape(B_scalar_model(:,2), s);
Bp_grid = reshape(B_scalar_model(:,3), s);
F_grid = reshape(F_scalar_model, s);
dBr_grid = reshape(dB_scalar(:,1), s);
dBt_grid = reshape(dB_scalar(:,2), s);
dBp_grid = reshape(dB_scalar(:,3), s);
dF_grid = reshape(dF_scalar, s);
Bchaos_grid = reshape(B_chaos(:,1), s);
Fchaos_grid = reshape(F_chaos, s);

%%


load coastlines

figure(2)


subplot(2,2,1)

worldmap(lat_rng, lon_rng);
pcolorm(lat_grid, lon_grid, Br_grid)
plotm(coastlat, coastlon, 'Color', 'black')
title('B_r (nT)')
contourcbar('southoutside')

subplot(2,2,2)

worldmap(lat_rng, lon_rng);
pcolorm(lat_grid, lon_grid, Bt_grid)
plotm(coastlat, coastlon, 'Color', 'black')
title('B_{\theta} (nT)')
contourcbar('southoutside')

subplot(2,2,3)

worldmap(lat_rng, lon_rng);
pcolorm(lat_grid, lon_grid, Bp_grid)
plotm(coastlat, coastlon, 'Color', 'black')
title('B_{\phi} (nT)')
contourcbar('southoutside')

subplot(2,2,4)

worldmap(lat_rng, lon_rng);
pcolorm(lat_grid, lon_grid, F_grid)
plotm(coastlat, coastlon, 'Color', 'black')
title('F (nT)')
contourcbar('southoutside')



figure(3)

subplot(2,2,1)

worldmap(lat_rng, lon_rng);
pcolorm(lat_grid, lon_grid, dBr_grid)
plotm(coastlat, coastlon, 'Color', 'black')
title('B_r - B_r_{CHAOS}')
contourcbar('southoutside')

subplot(2,2,2)

worldmap(lat_rng, lon_rng);
pcolorm(lat_grid, lon_grid, dBt_grid)
plotm(coastlat, coastlon, 'Color', 'black')
title('B_{\theta} - B_{\theta}_{CHAOS}')
contourcbar('southoutside')

subplot(2,2,3)

worldmap(lat_rng, lon_rng);
pcolorm(lat_grid, lon_grid, dBp_grid)
plotm(coastlat, coastlon, 'Color', 'black')
title('B_{\phi} - B_{\phi}_{CHAOS}')
contourcbar('southoutside')

subplot(2,2,4)

worldmap(lat_rng, lon_rng);
pcolorm(lat_grid, lon_grid, dF_grid)
plotm(coastlat, coastlon, 'Color', 'black')
title('F - F_{CHAOS}')
contourcbar('southoutside')




figure(4)

subplot(2,3,1)

worldmap(lat_rng, lon_rng);
pcolorm(lat_grid, lon_grid, Br_grid)
plotm(coastlat, coastlon, 'Color', 'black')
title('B_r (nT)')
contourcbar('southoutside')

subplot(2,3,2)

worldmap(lat_rng, lon_rng);
pcolorm(lat_grid, lon_grid, Bchaos_grid)
plotm(coastlat, coastlon, 'Color', 'black')
title('B_r_{CHAOS} (nT)')
contourcbar('southoutside')

subplot(2,3,3)

worldmap(lat_rng, lon_rng);
pcolorm(lat_grid, lon_grid, dBr_grid)
plotm(coastlat, coastlon, 'Color', 'black')
title('B_r - B_r_{CHAOS}')
contourcbar('southoutside')

subplot(2,3,4)

worldmap(lat_rng, lon_rng);
pcolorm(lat_grid, lon_grid, F_grid)
plotm(coastlat, coastlon, 'Color', 'black')
title('F (nT)')
contourcbar('southoutside')

subplot(2,3,5)

worldmap(lat_rng, lon_rng);
pcolorm(lat_grid, lon_grid, Fchaos_grid)
plotm(coastlat, coastlon, 'Color', 'black')
title('F_{CHAOS} (nT)')
contourcbar('southoutside')

subplot(2,3,6)

worldmap(lat_rng, lon_rng);
pcolorm(lat_grid, lon_grid, dF_grid)
plotm(coastlat, coastlon, 'Color', 'black')
title('F - F_{CHAOS}')
contourcbar('southoutside')




