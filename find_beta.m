function beta = find_beta(F_obs, r, theta, phi, g, N)
% F_synth: total magnetic flux density in nT
% r: geocentric radius in kilometers
% theta: colatitude in radians
% phi: longitude in radians
% g: Gauss coeffiecients
% N: maximum spherical harmonic degree

F = find_F(r, theta, phi, g, N);
beta = zeros(length(r), 1);
for i = 1:length(r)
    beta(i) = F_obs(i) - F(i);
end

return