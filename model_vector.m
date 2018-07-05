% length of g = # of columns in J = (N+1)^2-1

r = ones(1,10000) * 6500; %km
theta = rand(10000,1) * 180; %degrees
phi = rand(10000,1) * 360 - 180; %degrees
N = 10;

NSH = (N+1)^2-1;
g_synth = zeros(NSH,1);
g_synth(index(1,0)) = 1.0;

[Br_synth, Btheta_synth, Bphi_synth] = find_B(r, theta, phi, g_synth, N);


%%
g_init = zeros(NSH, 1);
J = find_J(r, theta, phi, N);
% W = eye(length(r)); % weighting function
gamma = 0.5;

for x = 1:20
    gprev = g_init;
    alpha = find_alpha(Br_synth, r, theta, phi, g_init, N);
    delta = J\alpha; %rcond(J'*J) is too low, close to singular
    g_init = g_init + gamma*delta;
    error = norm(g_init - gprev) / norm(g_init);
    error
end

final_error = norm(g_synth - g_init) / norm(g_synth)