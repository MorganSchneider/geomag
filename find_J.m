function [J_alpha, J_beta] = find_J(r, theta, phi, g, N)
% r: geocentric radius in kilometers
% theta: colatitude in radians
% phi: longitude in radians
% g: Gauss coefficients
% N: maximum spherical harmonic degree

a = 6371.2; %km

kr = length(r);
NSH = (N+1)^2 - 1;

J_alpha = zeros(kr, NSH);
J_beta = zeros(kr, NSH);
S = find_Snm(theta, phi, N);
dS_theta = find_dS_theta(theta, phi, N);
dS_phi = find_dS_phi(theta, phi, N);
F = find_F(r, theta, phi, g, N);
[Br, Bt, Bp] = find_B(r, theta, phi, g, N);
for i = 1:kr
    for n = 1:N
        for m = -n:n
            j = index(n,m);
            J_alpha(i,j) = (a/r(i))^(n+2) * (n+1) * S(j,i);
            dBr = J_alpha(i,j);
            dBt = (a/r(i))^(n+2) * -dS_theta(j,i);
            dBp = (a/r(i))^(n+2) * -dS_phi(j,i);
            J_beta(i,j) = F(i)^-1 * ((Br(i)*dBr) + (Bt(i)*dBt) + (Bp(i)*dBp));
        end
    end
end


return