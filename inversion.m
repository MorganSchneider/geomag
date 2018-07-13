% length of g = # of columns in J = (N+1)^2-1
rad = pi/180; %radians

% load('./EEJ_Data/Swarm_1HzData.mat')
% [~, lat1, lon1, r1, ~, ~, ~] = find_EEJ(swarm, 1, 'interp');
% [~, lat2, lon2, r2, ~, ~, ~] = find_EEJ(swarm, 2, 'interp');
% [~, lat3, lon3, r3, ~, ~, ~] = find_EEJ(swarm, 3, 'interp');
% r_r = [r1, r2, r3];
% theta_r = (90 - [lat1, lat2, lat3]) * rad;
% phi_r = [lon1, lon2, lon3] * rad;
% load('./EEJ_Data/Swarm_Scalar.mat')
% r_s = [scalar(1).rad, scalar(2).rad];
% theta_s = (90 - [scalar(1).geolat, scalar(2).geolat]) * rad;
% phi_s = [scalar(1).lon, scalar(2).lon] * rad;

r_r = ones(1,10000) * 6500; %km
theta_r = rand(10000,1) * 180; theta_r = theta_r * rad; %radians
phi_r = rand(10000,1) * 360 - 180; phi_r = phi_r * rad; %radians
r_s = ones(1,1000) * 6500; %km
theta_s = rand(1000,1) * 180; theta_s = theta_s * rad; %radians
phi_s = rand(1000,1) * 360 - 180; phi_s = phi_s * rad; %radians

N = 10;
NSH = (N+1)^2-1;
g_synth = -5 + rand(NSH,1) * 10;
[Br_synth, ~, ~] = find_B(r_r, theta_r, phi_r, g_synth, N);
F_synth = find_F(r_s, theta_s, phi_s, g_synth, N);


% Br_EEJ = zeros(length(r_r), 1);
% F_swarm = rot90([scalar(1).F, scalar(2).F]);

%%

g_init = ones(NSH,1);
% W = weighting function;
gamma = 0.5;
[J_alpha,~] = find_J(r_r, theta_r, phi_r, g_init, N);

for x = 1:50
    gprev = g_init;
    alpha = find_alpha(Br_synth, r_r, theta_r, phi_r, g_init, N);
    beta = find_beta(F_synth, r_s, theta_s, phi_s, g_init, N);
    [~,J_beta] = find_J(r_s, theta_s, phi_s, g_init, N);
    J = vertcat(J_alpha, J_beta);
    R = vertcat(alpha, beta);
    delta = (J'*J)\(J'*R);
    g_init = g_init + gamma*delta;
    error = norm(g_init - gprev) / norm(g_init)
end
final_error = norm(g_synth - g_init) / norm(g_synth)