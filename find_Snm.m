function S = find_Snm(theta, phi, N)
% theta: colatitude in radians
% phi: longitude in radians
% N: maximum spherical harmonic degree

maxK = index(N,N);
S = zeros(maxK, length(theta));
for n = 1:N
    P = legendre(n, cos(theta), 'sch');
    for i = 1:length(theta)
        for m = -n:n
            k = index(n,m);
            if m >= 0
                S(k,i) = cos(m*phi(i)) * P(m+1,i);
            elseif m < 0
                S(k,i) = sin(abs(m)*phi(i)) * P(abs(m)+1,i);
            end
        end
    end
end

return