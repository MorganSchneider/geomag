function J_alpha = find_J(r, theta, phi, N)
% right now, only calculates the Jacobian for the radial component.
a = 6371.2;
rad = pi/180;

kr = length(r);
NSH = (N+1)^2 - 1;

J_alpha = zeros(kr, NSH);
S = find_Snm(theta*rad, phi*rad, N);
for i = 1:kr
    for n = 1:N
        for m = -n:n
            j = index(n,m);
            J_alpha(i,j) = (a/r(i))^(n+2) * (n+1) * S(j,i);
        end
    end
end

return