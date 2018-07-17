function year = decimalYear(t)
% Convert from Unix time to decimal year.

[y,mo,d] = ymd(datetime(t, 'ConvertFrom', 'posixtime'));
[h,mi,s] = hms(datetime(t, 'ConvertFrom', 'posixtime'));


sFrac = s ./ 60;
mi = mi + sFrac;
miFrac = mi ./ 60;
h = h + miFrac;
hFrac = h ./ 24;
d = d + hFrac;

dFrac = zeros(1, length(t));
for i = 1:length(t)
    if mo(i)==2 && y(i)==2016
        dFrac(i) = d(i) / 29;
    elseif mo(i)==2 && y(i)~=2016
        dFrac(i) = d(i) / 28;
    elseif mo(i)==4 || mo(i)==6 || mo(i)==9 || mo(i)==11
        dFrac(i) = d(i) / 30;
    else
        dFrac(i) = d(i) / 31;
    end
end

mo = mo + dFrac;
moFrac = mo ./ 12;
year = y + moFrac;

return