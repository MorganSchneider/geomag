function B = find_B(r, theta, phi, g, N)
%
%   INPUTS:
%   r:          geocentric radius (km)
%   theta:      geocentric colatitude (degrees)
%   phi:        geocentric longitude (degrees)
%   g:          numerical vector of Gauss coefficients (reshaped using index)
%   N:          maximum spherical harmonic degree
%
%   OUTPUT:
%   B:          magnetic field components at position (r, theta, phi)

a = 6371.2; %km
rad = pi/180; %radians

if length(r) ~= length(theta) || length(r) ~= length(phi)
    fprintf('Position vectors must be of equal length.')
    return
end
    

% normalize radius by Earth's radius
% convert colatitude and longitude from degrees to radians
r_n = a * r.^-1;
theta = theta * rad;
phi = phi * rad;

B = zeros(length(r), 3);
for i = 1:length(r)
    S = find_Snm(theta(i), phi(i), N);
    B_r = zeros(1, (N+1)^2-1);
    B_theta = zeros(1, (N+1)^2-1);
    B_phi = zeros(1, (N+1)^2-1);
    for n = 1:N
        for m = -n:n
            k = index(n,m);
            B_r(k) = g(k) * r_n(i)^(n+2) * (n+1) * S(k);
            B_theta(k) = g(k) * r_n(i)^(n+2) * dS; %% how to calculate derivatives of S?
            B_phi(k) = g(k) * r_n(i)^(n+2) * dS;
        end
    end
    B(i,1) = sum(B_r);
    B(i,2) = sum(B_theta);
    B(i,3) = sum(B_phi);
end




return