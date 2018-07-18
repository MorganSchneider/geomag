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

[t1, lat1, lon1, r1, ~, ~, ~] = find_EEJ(swarm, 1);
[t2, lat2, lon2, r2, ~, ~, ~] = find_EEJ(swarm, 2);
[t3, lat3, lon3, r3, ~, ~, ~] = find_EEJ(swarm, 3);

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

Br_EEJ = zeros(length(r_r), 1);
F_swarm = rot90([scalar(1).F, scalar(2).F], 3); F_swarm = F_swarm(inds_s);

%% Invert

g_init = ones(NSH,1);
g_init(index(1,0)) = -30000; % Earth's dipole moment is roughly 30000 nT
gamma = 0.5;
[J_alpha,~] = find_J(r_r, theta_r, phi_r, g_init, N);

% Build weighting functions
W_r = diag(ones(1,length(r_r)) * 30);
W_s = zeros(length(r_s), length(r_s));
for i = 1:length(r_s)
    W_s(i,i) = sin(theta_s(i));
end


for x = 1:30
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



