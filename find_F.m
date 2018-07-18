function F = find_F(r, theta, phi, g, N)
% r: geocentric radius in kilometers
% theta: colatitude in radians
% phi: longitude in radians
% g: Gauss coeffiecients
% N: maximum spherical harmonic degree (if CHAOS, N is the vector of mjd times)

if length(N) > 1
    t = N;
    B = synth_values(r, theta, phi, g, t);
else
    B = find_B(r, theta, phi, g, N);
end
F = zeros(length(r), 1);
for i = 1:length(r)
    F(i) = norm([B(i,1), B(i,2), B(i,3)]);
end

return