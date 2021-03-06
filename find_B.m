function B = find_B(r, theta, phi, g, N)
%
%   INPUTS:
%   r:          geocentric radius (km)
%   theta:      geocentric colatitude (radians)
%   phi:        geocentric longitude (radians)
%   g:          numerical vector of Gauss coefficients
%   N:          maximum spherical harmonic degree
%
%   OUTPUT:
%   B:          magnetic field components at position (r, theta, phi)

a = 6371.2; %km

if length(r) ~= length(theta) || length(r) ~= length(phi)
    fprintf('Position vectors must be of equal length.')
    return
end

r = a * r.^-1;

B = zeros(length(r), 3);

S = find_Snm(theta, phi, N);
dS_theta = find_dS_theta(theta, phi, N);
dS_phi = find_dS_phi(theta, phi, N);
for i = 1:length(r)
    B_r_i = zeros(1, (N+1)^2-1);
    B_theta_i = zeros(1, (N+1)^2-1);
    B_phi_i = zeros(1, (N+1)^2-1);
    for n = 1:N
        for m = -n:n
            k = index(n,m);
            B_r_i(k) = g(k) * r(i)^(n+2) * (n+1) * S(k,i);
            B_theta_i(k) = g(k) * r(i)^(n+2) * -dS_theta(k,i);
            B_phi_i(k) = g(k) * r(i)^(n+2) * -dS_phi(k,i);
        end
    end
    B(i,1) = sum(B_r_i);
    B(i,2) = sum(B_theta_i);
    B(i,3) = sum(B_phi_i);
end




return