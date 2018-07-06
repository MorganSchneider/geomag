function alpha = find_alpha(Br_synth, r, theta, phi, g, N)
% Br_synth: radial magnetic flux density in nT
% r: geocentric radius in kilometers
% theta: colatitude in radians
% phi: longitude in radians
% g: Gauss coeffiecients
% N: maximum spherical harmonic degree

[Br, ~, ~] = find_B(r, theta, phi, g, N);
alpha = zeros(length(r), 1);
for i = 1:length(r)
    alpha(i) = Br_synth(i) - Br(i); % once I start using actual data, Br_synth will be Br_EEJ (= 0).
end

return