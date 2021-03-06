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
B = find_B(r, theta, phi, g, N);
for i = 1:kr
    ar = a / r(i);
    for n = 1:N
        for m = -n:n
            j = index(n,m);
            J_alpha(i,j) = ar^(n+2) * (n+1) * S(j,i);
            dBr = J_alpha(i,j);
            dBt = ar^(n+2) * -dS_theta(j,i);
            dBp = ar^(n+2) * -dS_phi(j,i);
            J_beta(i,j) = (F(i)^-1) * ((B(i,1)*dBr) + (B(i,2)*dBt) + (B(i,3)*dBp));
        end
    end
end


return