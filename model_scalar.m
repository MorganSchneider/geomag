% it isn't working and I think the problem is in find_J

rad = pi/180; %radians

r = ones(1,1000) * 6500; %km
theta = rand(1000,1) * 180; theta = theta * rad; %radians
phi = rand(1000,1) * 360 - 180; phi = phi * rad; %radians
N = 10;

NSH = (N+1)^2-1;
g_synth = zeros(NSH,1);
g_synth(index(1,0)) = 1.0;

F_synth = find_F(r, theta, phi, g_synth, N);

%%

g_init = ones(NSH, 1) * 0.5;
% W = weighting function;
gamma = 0.2;

for x = 1:50
    gprev = g_init;
    beta = find_beta(F_synth, r, theta, phi, g_init, N);
    [~,J] = find_J(r, theta, phi, g_init, N);
    delta = J\beta;
    g_init = g_init + gamma*delta;
    error = norm(g_init - gprev) / norm(g_init)
    total_error = norm(g_synth - g_init) / norm(g_synth)
end
%final_error = norm(g_synth - g_init) / norm(g_synth)