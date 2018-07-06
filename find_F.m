function F = find_F(r, theta, phi, g, N)
% r: geocentric radius in kilometers
% theta: colatitude in radians
% phi: longitude in radians
% g: Gauss coeffiecients
% N: maximum spherical harmonic degree

[Br, Bt, Bp] = find_B(r, theta, phi, g, N);
F = zeros(length(r), 1);
for i = 1:length(r)
    B = [Br(i), Bt(i), Bp(i)];
    F(i) = norm(B);
end

return