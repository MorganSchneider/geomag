function dg_mat = difference_mat(g_1, g_2)

% Routine for evaluating the difference between two sets of Gauss
% coefficients.

if length(g_1) ~= length(g_2)
    error('Models must have the same number of coefficients')
end

NSH = length(g_1);
N = sqrt(NSH + 1) - 1;

dg = g_1 - g_2;
dg_mat = NaN(N, 1+2*N);
for n = 1:N
    for m = -n:n
        dg_mat(n,m+N+1) = dg(index(n,m));
    end
end

return