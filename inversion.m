%% Load data

load('./EEJ_Data/Swarm_1HzData.mat')
load('./EEJ_Data/Swarm_Scalar.mat')
load('./CHAOS-6_FWD/CHAOS-6-x5.mat')

%% Set CHAOS coefficients

N = 20;
NSH = (N+1)^2-1;
pp_N = pp;
pp_N.dim = N*(N+2);
coefs_tmp = reshape(pp.coefs, [], pp.pieces, pp.order);
pp_N.coefs = reshape(coefs_tmp(1:N*(N+2),:,:), [], pp.order);

%% Retrieve EEJ peaks

[t1, lat1, lon1, r1, qd1, ~, ~] = find_EEJ(swarm, 1);
[t2, lat2, lon2, r2, qd2, ~, ~] = find_EEJ(swarm, 2);
[t3, lat3, lon3, r3, qd3, ~, ~] = find_EEJ(swarm, 3);

%% Separate times

t_EEJ = [t1, t3];
t_scalar = [scalar(1).time, scalar(2).time];
yearEEJ = decimalYear(t_EEJ);
yearScalar = [scalar(1).year, scalar(2).year];

timeStart = 2015.0;
timeEnd = 2015.25;

inds_r = find(yearEEJ >= timeStart & yearEEJ < timeEnd);
inds_s = find(yearScalar >= timeStart & yearScalar < timeEnd);
t_r = t_EEJ(inds_r);
t_s = t_scalar(inds_s);

mjd_r = datenum(datetime(t_r,'ConvertFrom','posixtime')) - datenum(2000,1,1,0,0,0);
mjd_s = datenum(datetime(t_s,'ConvertFrom','posixtime')) - datenum(2000,1,1,0,0,0);

%% Create vectors for radial and scalar data

% length of g = # of columns in J = (N+1)^2-1
rad = pi/180; %radians

r_r = [r1, r3]; r_r = r_r(inds_r);
theta_r = (90 - [lat1, lat3]) * rad; theta_r = theta_r(inds_r);
phi_r = [lon1, lon3] * rad; phi_r = phi_r(inds_r);
r_s = [scalar(1).rad, scalar(2).rad]; r_s = r_s(inds_s);
theta_s = (90 - [scalar(1).geolat, scalar(2).geolat]) * rad; theta_s = theta_s(inds_s);
phi_s = [scalar(1).lon, scalar(2).lon] * rad; phi_s = phi_s(inds_s);

qd_r = [qd1, qd3]; qd_r = qd_r(inds_r);
qd_s = [scalar(1).qdlat, scalar(2).qdlat]; qd_s = qd_s(inds_s);

%Br_EEJ = zeros(length(r_r), 1);
B_EEJ = synth_values(r_r, theta_r./rad, phi_r./rad, pp_N, mjd_r);
Br_EEJ = B_EEJ(:,1);
F_swarm = rot90([scalar(1).F, scalar(2).F], 3); F_swarm = F_swarm(inds_s);

%% Invert

g_init = zeros(NSH,1);
g_init(index(1,0)) = -30000; % Earth's dipole moment is roughly 30000 nT
gam = 0.5;
[J_alpha,~] = find_J(r_r, theta_r, phi_r, g_init, N);

% Build weighting functions
w_r = diag(ones(1,length(r_r)) * 10);
w_s = zeros(length(r_s), length(r_s));
for i = 1:length(r_s)
    w_s(i,i) = sin(theta_s(i));
end

w_robust_r = ones(length(r_r), 1);
w_robust_s = ones(length(r_s), 1);
for x = 1:50
    gprev = g_init;
    alpha = find_alpha(Br_EEJ, r_r, theta_r, phi_r, g_init, N);
    beta = find_beta(F_swarm, r_s, theta_s, phi_s, g_init, N);
    [~,J_beta] = find_J(r_s, theta_s, phi_s, g_init, N);
    
    % robust statistical removal of outliers from EEJ dataset
    if x > 1
        rms_r = sqrt(alpha' * alpha / length(alpha));
        rms_s = sqrt(beta' * beta / length(beta));
        w_robust_r = min(1.5 ./ abs(alpha./rms_r), 1);
        w_robust_s = min(1.5 ./ abs(beta./rms_s), 1);
    end
    
    W_r = diag(w_r * w_robust_r);
    W_s = diag(w_s * w_robust_s);
    
    JT_W_J_alpha = J_alpha' * W_r * J_alpha; %left side of radial inversion
    JT_W_alpha = J_alpha' * W_r * alpha; %right side of radial inversion
    JT_W_J_beta = J_beta' * W_s * J_beta; %left side of scalar inversion
    JT_W_beta = J_beta' * W_s * beta; % right side of scalar inversion
    
    JT_W_J = JT_W_J_alpha + JT_W_J_beta; %left side of total inversion
    JT_W_R = JT_W_alpha + JT_W_beta; %right side of total inversion
    
    delta = (JT_W_J)\(JT_W_R);
    g_init = g_init + gam*delta;
    error = norm(g_init - gprev) / norm(g_init)
    chisq = norm(sqrt(W_r) * alpha) + norm(sqrt(W_s) * beta)
end
g_model = g_init; % final vector of Gauss coeffs

%%

B_model_r = find_B(r_r, theta_r, phi_r, g_model, N);
F_model_s = find_F(r_s, theta_s, phi_s, g_model, N);
lon_deg_r = phi_r./rad;

close all

figure(1)
subplot(2,1,1)
plot(lon_deg_r, Br_EEJ - B_model_r(:,1), '*')
xlabel('Longitude (deg)')
ylabel('B_r Residuals (nT)')
title('B_r_{Swarm} - B_r_{Model}, 2015.0 - 2015.25')
subplot(2,1,2)
plot(qd_s, F_swarm - F_model_s, '*')
xlabel('Quasi-Dipole Latitude (deg)')
ylabel('F Residuals (nT)')
title('F_{Swarm} - F_{Model}, 2015.0 - 2015.25')


g_chaos = reindex(pp_N, N, mean(mjd_s));

dg = difference_mat(g_model, g_chaos);
Sm = sensitivity_mat(g_model, g_chaos);
Ps = powerspec(g_model);


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
t_const = ones(1, 80000) * mean(mjd_s);

s = size(lat_grid);

B_model = find_B(r_const, lat_vector*rad, lon_vector*rad, g_model, N);
F_model = find_F(r_const, lat_vector*rad, lon_vector*rad, g_model, N);
B_chaos = synth_values(r_const, lat_vector, lon_vector, pp_N, t_const);
F_chaos = find_F(r_const, lat_vector, lon_vector, pp_N, t_const);
dB = B_model - B_chaos;
dF = F_model - F_chaos;

Br_grid = reshape(B_model(:,1), s);
Bt_grid = reshape(B_model(:,2), s);
Bp_grid = reshape(B_model(:,3), s);
F_grid = reshape(F_model, s);
dBr_grid = reshape(dB(:,1), s);
dBt_grid = reshape(dB(:,2), s);
dBp_grid = reshape(dB(:,3), s);
dF_grid = reshape(dF, s);
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

%% Statistics



theta_A = (90 - scalar(1).geolat(inds_A)) * rad;
phi_A = scalar(1).lon(inds_A) * rad;
theta_B = (90 - scalar(2).geolat(inds_B)) * rad;
phi_B = scalar(2).lon(inds_B) * rad;
r_A = ones(1, length(theta_A)) * 6371.2;
r_B = ones(1, length(theta_B)) * 6371.2;
t_A = ones(1, length(theta_A)) * mean(mjd_s);
t_B = ones(1, length(theta_B)) * mean(mjd_s);

r = [r_A r_B];
theta = [theta_A theta_B];
phi = [phi_A phi_B];
t = [t_A t_B];

B_model = find_B(r, theta, phi, g_model, N);
B_chaos = synth_values(r, theta./rad, phi./rad, pp_N, t);
dB = B_model - B_chaos;
dBr = dB(:,1);

hilat_Br = find(abs(90 - theta./rad) > 55);
lolat_Br = find(abs(90 - theta./rad) <= 55);
dBr_hilat = dBr(hilat_Br);
dBr_lolat = dBr(lolat_Br);

mu_Br_hilat = mean(dBr_hilat);
mu_Br_lolat = mean(dBr_lolat);
sigma_Br_hilat = std(dBr_hilat);
sigma_Br_lolat = std(dBr_lolat);
rms_Br_hilat = sqrt(dBr_hilat' * dBr_hilat / length(dBr_hilat));
rms_Br_lolat = sqrt(dBr_lolat' * dBr_lolat / length(dBr_lolat));

statsBr.hilat = struct('mu', mu_Br_hilat, 'sigma', sigma_Br_hilat, 'rms', rms_Br_hilat);
statsBr.lolat = struct('mu', mu_Br_lolat, 'sigma', sigma_Br_lolat, 'rms', rms_Br_lolat);

F_model_A = find_F(r_A, theta_A, phi_A, g_model, N);
F_chaos_A = find_F(r_A, theta_A./rad, phi_A./rad, pp_N, t_A);
F_model_B = find_F(r_B, theta_B, phi_B, g_model, N);
F_chaos_B = find_F(r_B, theta_B./rad, phi_B./rad, pp_N, t_B);

hilat_A = find(abs(90 - theta_A./rad) > 55);
lolat_A = find(abs(90 - theta_A./rad) <= 55);
hilat_B = find(abs(90 - theta_B./rad) > 55);
lolat_B = find(abs(90 - theta_B./rad) <= 55);

F_model_A_hilat = F_model_A(hilat_A);
F_chaos_A_hilat = F_chaos_A(hilat_A);
F_model_A_lolat = F_model_A(lolat_A);
F_chaos_A_lolat = F_chaos_A(lolat_A);

F_model_B_hilat = F_model_B(hilat_B);
F_chaos_B_hilat = F_chaos_B(hilat_B);
F_model_B_lolat = F_model_B(lolat_B);
F_chaos_B_lolat = F_chaos_B(lolat_B);

dF_A_hilat = F_model_A_hilat - F_chaos_A_hilat;
dF_A_lolat = F_model_A_lolat - F_chaos_A_lolat;
dF_B_hilat = F_model_B_hilat - F_chaos_B_hilat;
dF_B_lolat = F_model_B_lolat - F_chaos_B_lolat;

mu_A_hilat = mean(dF_A_hilat);
mu_A_lolat = mean(dF_A_lolat);
mu_B_hilat = mean(dF_B_hilat);
mu_B_lolat = mean(dF_B_lolat);

sigma_A_hilat = std(dF_A_hilat);
sigma_A_lolat = std(dF_A_lolat);
sigma_B_hilat = std(dF_B_hilat);
sigma_B_lolat = std(dF_B_lolat);

rms_A_hilat = sqrt(dF_A_hilat' * dF_A_hilat / length(dF_A_hilat));
rms_A_lolat = sqrt(dF_A_lolat' * dF_A_lolat / length(dF_A_lolat));
rms_B_hilat = sqrt(dF_B_hilat' * dF_B_hilat / length(dF_B_hilat));
rms_B_lolat = sqrt(dF_B_lolat' * dF_B_lolat / length(dF_B_lolat));

statsF.A.hilat = struct('mu', mu_A_hilat, 'sigma', sigma_A_hilat, 'rms', rms_A_hilat);
statsF.A.lolat = struct('mu', mu_A_lolat, 'sigma', sigma_A_lolat, 'rms', rms_A_lolat);
statsF.B.hilat = struct('mu', mu_B_hilat, 'sigma', sigma_B_hilat, 'rms', rms_B_hilat);
statsF.B.lolat = struct('mu', mu_B_lolat, 'sigma', sigma_B_lolat, 'rms', rms_B_lolat);

save('./stats2015_00_25.mat', 'statsBr', 'statsF', '-v7.3', '-nocompression')


