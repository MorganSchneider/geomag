function Sm = sensitivity_mat(m_1,m_2)

% Sm = sensitivity_mat(m_1,m_2)
%
% Sensitivity matrix between two SH models
% 
% A. Chulliat, 2016-11-23

N_koeff = length(m_1);
if (length(m_2) ~= N_koeff)
    error('the models don''t have the same number of Gauss coefficients')
end

N = sqrt(N_koeff+1)-1;

R_n = powerspec(m_2);
dm  = m_1 - m_2;

Sm = NaN(N,2*N+1);
k = 1;
for n = 1:N
    for m = 0:n
        Sm(n,N+1+m) = 100* dm(k)./sqrt(R_n(n)./(n+1)./(2*n+1));
        k = k+1;
        if m > 0
            Sm(n,N+1-m) = 100*dm(k)./sqrt(R_n(n)./(n+1)./(2*n+1));
            k = k+1;
        end
    end
end




