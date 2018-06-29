function J = find_J(r, theta, phi, N)

for i = 1:length(r)
    S = find_Snm(theta, phi, N);
    for n = 1:N
        for m = -n:n
            j = index(n,m);
            J(i,j) = r(i)^(n+2) * (n+1) * S(j);
        end
    end
end

return