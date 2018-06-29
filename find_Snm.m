function S = find_Snm(theta, phi, N)

for n = 1:N
    for m = -n:n
        k = index(n,m);
        P = legendre(n, cos(theta));
        if m >= 0
            S(1:length(P),k) = cos(m*phi) * P;
        elseif m < 0
            S(1:length(P),k) = sin(abs(m)*phi) * P;
        end
    end
end

return