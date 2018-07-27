syms theta phi
s(theta) = legendreP(10, cosd(theta));
ds_theta = diff(s, theta);
ds_theta(angle_deg);

