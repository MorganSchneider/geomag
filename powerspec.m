function R_n = powerspec(m_i, varargin)
% R_n = powerspec(m_i)
% Lowe-Mauersberger spectrum of Gauss coefficients m_i at radius r/a = 1
% R_n = powerspec(m_i, rho) yields spectrum at radius r/a = rho
% 
% May 2000, Nils Olsen, DSRI, Copenhagen

if nargin == 1
   rho = 1;
else
   rho = varargin{1};
end

if nargin > 1
    source = varargin{2};
else
    source = 'int'; % internal sources by default
end
    
maxdeg_MF = sqrt(length(m_i)+1)-1;
R_n = zeros(maxdeg_MF, 1);
i=0;
for n = 1:maxdeg_MF
   for m = 0:n
      i=i+1;
      R_n(n) = R_n(n)+m_i(i).^2;
      if m > 0
         i=i+1;
         R_n(n) = R_n(n)+m_i(i).^2;
      end
   end
   if strcmp(source, 'int')   
       R_n(n) = R_n(n)*(n+1)*rho.^-(2*n+4);
   elseif strcmp(source, 'ext')   
       R_n(n) = R_n(n)*n*rho.^(2*n-4);
   elseif strcmp(source, 'tor')   
       R_n(n) = R_n(n)*n*(n+1)/(2*n+1)*rho.^-2;
   else
       warning 'Unrecognized source: not ''int'', ''ext'', or ''tor'''
   end
end
