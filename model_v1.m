% length of g = # of columns in J = (N+1)^2-1

a = 6371.2;
rad = pi/180;

N = 10;
NSH = (N+1)^2-1;

r = 6500; %km
r_n = a./r;
theta = 90; %degrees
theta_n = theta * rad;
phi = 90; %degrees
phi_n = phi * rad;
g_synth = ones(1, NSH) .* rand([1 NSH]);

B_synth = find_B(r_n, theta_n, phi_n, g_synth, N);
Br_synth = B_synth(:,1);

%%

g_init = zeros(1, NSH);








