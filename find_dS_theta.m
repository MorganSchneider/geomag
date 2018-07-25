function dS = find_dS_theta(theta, phi, N)
% theta: colatitude in radians
% phi: longitude in radians
% N: maximum spherical harmonic degree

maxK = index(N,N);
dS = zeros(maxK, length(theta));
h = 1e-10;
theta(theta==0) = nan;
theta(theta==pi) = nan;
for n = 1:N
    P_plus = legendre(n, cos(theta) + h, 'sch');
    P_minus = legendre(n, cos(theta) - h, 'sch');
    dP = (P_plus - P_minus) ./ (2*h);
    for i = 1:length(theta)
        for m = -n:n
            k = index(n,m);
            if m >= 0
                dS(k,i) = cos(m*phi(i)) * dP(m+1,i) * -sin(theta(i));
            elseif m < 0
                dS(k,i) = sin(abs(m)*phi(i)) * dP(abs(m)+1,i) * -sin(theta(i));
            end
        end
    end
end

return