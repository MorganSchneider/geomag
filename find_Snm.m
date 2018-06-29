function S = find_Snm(theta, phi, N)

S = zeros(1, (N+1)^2-1);
for n = 1:N
    for m = -n:n
        k = index(n,m);
        P = legendre(n, cos(theta));
        if m >= 0
            S(k) = cos(m*phi) * P;
        elseif m < 0
            S(k) = sin(abs(m)*phi) * P;
        end
    end
end

return