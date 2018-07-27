function g_new = reindex(pp, N, t)

% Routine to reorder the CHAOS Gauss coefficients to match my indexing
% scheme.

g_old = fnval(t, pp);

order = zeros(1, 1+2*N);
for i = 1:N
    order(2*i) = i;
    order(2*i+1) = -i;
end
NSH = (N+1)^2-1;

l = 1;
g_new = zeros(NSH, 1);
nn = [];
mm = [];
for n = 1:N
    num = 1 + 2 * n;
    nn = [nn, n * ones(1,num)];
    mm = [mm, -n:n];
    for x = 1:num
        k = index(n, order(x));
        g_new(k) = g_old(l);
        l = l + 1;
    end
end

return