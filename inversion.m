%% Load data

load('./EEJ_Data/Swarm_1HzData.mat')
load('./EEJ_Data/Swarm_Scalar.mat')
load('./CHAOS-6_FWD/CHAOS-6-x5.mat')

N = 10;
pp_N = pp;
pp_N.dim = N*(N+2);
coefs_tmp = reshape(pp.coefs, [], pp.pieces, pp.order);
pp_N.coefs = reshape(coefs_tmp(1:N*(N+2),:,:), [], pp.order);

%% Retrieve EEJ peaks

[t1, lat1, lon1, r1, qd1, ~, ~] = find_EEJ(swarm, 1);
[t2, lat2, lon2, r2, qd2, ~, ~] = find_EEJ(swarm, 2);
[t3, lat3, lon3, r3, qd3, ~, ~] = find_EEJ(swarm, 3);

%% Separate times
t_EEJ = [t1, t2, t3];
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
NSH = (N+1)^2-1;

r_r = [r1, r2, r3]; r_r = r_r(inds_r);
theta_r = (90 - [lat1, lat2, lat3]) * rad; theta_r = theta_r(inds_r);
phi_r = [lon1, lon2, lon3] * rad; phi_r = phi_r(inds_r);
r_s = [scalar(1).rad, scalar(2).rad]; r_s = r_s(inds_s);
theta_s = (90 - [scalar(1).geolat, scalar(2).geolat]) * rad; theta_s = theta_s(inds_s);
phi_s = [scalar(1).lon, scalar(2).lon] * rad; phi_s = phi_s(inds_s);

qd_r = [qd1, qd2, qd3]; qd_r = qd_r(inds_r);
qd_s = [scalar(1).qdlat, scalar(2).qdlat]; qd_s = qd_s(inds_s);

B_EEJ = synth_values(r_r, theta_r./rad, phi_r./rad, pp_N, mjd_r);
Br_EEJ = B_EEJ(:,1);
F_swarm = rot90([scalar(1).F, scalar(2).F], 3); F_swarm = F_swarm(inds_s);

%% Invert

g_init = zeros(NSH,1);
g_init(index(1,0)) = -30000; % Earth's dipole moment is roughly 30000 nT
gamma = 0.5;
[J_alpha,~] = find_J(r_r, theta_r, phi_r, g_init, N);

% Build weighting functions
W_r = diag(ones(1,length(r_r)) * 10);
W_s = zeros(length(r_s), length(r_s));
for i = 1:length(r_s)
    W_s(i,i) = sin(theta_s(i));
end


for x = 1:50
    gprev = g_init;
    alpha = find_alpha(Br_EEJ, r_r, theta_r, phi_r, g_init, N);
    beta = find_beta(F_swarm, r_s, theta_s, phi_s, g_init, N);
    [~,J_beta] = find_J(r_s, theta_s, phi_s, g_init, N);
    
    JT_W_J_alpha = J_alpha' * W_r * J_alpha; %left side of radial inversion
    JT_W_alpha = J_alpha' * W_r * alpha; %right side of radial inversion
    JT_W_J_beta = J_beta' * W_s * J_beta; %left side of scalar inversion
    JT_W_beta = J_beta' * W_s * beta; % right side of scalar inversion
    
    JT_W_J = JT_W_J_alpha + JT_W_J_beta; %left side of total inversion
    JT_W_R = JT_W_alpha + JT_W_beta; %right side of total inversion
    
    delta = (JT_W_J)\(JT_W_R);
    g_init = g_init + gamma*delta;
    error = norm(g_init - gprev) / norm(g_init)
end
g_model = g_init; % final vector of Gauss coeffs

%%

B_model = find_B(r_r, theta_r, phi_r, g_model, N);
F_model = find_F(r_s, theta_s, phi_s, g_model, N);

% compare to CHAOS
B_chaos = synth_values(r_r, theta_r, phi_r, pp_N, mjd_r);
F_chaos = find_F(r_s, theta_s, phi_s, pp_N, mjd_s);


%% Plot
subplot(2,1,1)
plot(lon_deg_r, Br_EEJ - B_model(:,1), '*')
subplot(2,1,2)
plot(qd_s, F_swarm - F_model, '*')

%%

Br = B_model(:,1);
Bt = B_model(:,2);
Bp = B_model(:,3);
deltaBr = B_model(:,1) - B_chaos(:,1);
deltaBt = B_model(:,2) - B_chaos(:,2);
deltaBp = B_model(:,3) - B_chaos(:,3);
deltaF = F_model - F_chaos;
lat_deg = 90 - (theta_s./rad);
lon_deg = phi_s./rad;
lon_deg_r = phi_r./rad;
%%
figure(1)
subplot(2,1,1)
plot(lon_deg, B_model(:,1) - B_chaos(:,1), '*')
xlabel('Longitude (deg)')
ylabel('\alpha (nT)')
title('B_{r} Model Residuals')
subplot(2,1,2)
plot(qd_s, F_model - F_chaos, '*')
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
F_grid = griddata(lat_deg, lon_deg, F_model, lat_grid, lon_grid, 'cubic');
dBr_grid = griddata(lat_deg, lon_deg, deltaBr, lat_grid, lon_grid, 'cubic');
dBt_grid = griddata(lat_deg, lon_deg, deltaBt, lat_grid, lon_grid, 'cubic');
dBp_grid = griddata(lat_deg, lon_deg, deltaBp, lat_grid, lon_grid, 'cubic');
dF_grid = griddata(lat_deg, lon_deg, deltaF, lat_grid, lon_grid, 'cubic');

% why
figure(2)
load coastlines

subplot(2,2,1)

worldmap(lat_rng, lon_rng);
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




