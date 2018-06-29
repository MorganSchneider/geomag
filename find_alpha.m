function alpha = find_alpha(Br_synth, r, theta, phi, g_init, N)

B_i = find_B(r, theta, phi, g_init, N);
Br_i = B_i(:,1);

alpha = Br_synth - Br_i;


return