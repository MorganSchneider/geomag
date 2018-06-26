orbits = find(nPeaks == 1);
j = 3;
for i = 1:length(orbits)
    st(j).orb(i).times = orbit(orbits(i)).time;
    st(j).orb(i).locs = orbit(orbits(i)).local;
    st(j).orb(i).rads = orbit(orbits(i)).rad;
    st(j).orb(i).lons = orbit(orbits(i)).lon;
    st(j).orb(i).glats = orbit(orbits(i)).geolat;
    st(j).orb(i).qlats = orbit(orbits(i)).qdlat;
    st(j).orb(i).f1 = orbit(orbits(i)).F1;
    st(j).orb(i).f2 = orbit(orbits(i)).F2;
    st(j).orb(i).meanf2 = orbit(orbits(i)).qd_meanF2;
    st(j).orb(i).stdvf2 = orbit(orbits(i)).qd_stdvF2;
    st(j).orb(i).df1 = orbit(orbits(i)).dF1;
    st(j).orb(i).kp = kpC(i);
    st(j).orb(i).rc = rcC(i);
    st(j).orb(i).res = resC(i);
    st(j).orb(i).res_corr = resC_corr(i);
end

%% Testing find_EEJ algorithm accuracy

% indices = find(nPeaksA > 1);
% spec = find(resA > 1);
%FOR A:
%sketchy: 3, 4, 6, 8, 10, 13, 14, 15, 16, 20, 32, 40, 43, 47, 52, 59, 60, 74,
%75, 80, 81, 83, 87, 89, 99, 100, 106, 115, 116
%extra-sketch: 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 51, 65, 69, 77,
%78, 82, 84, 85, 86, 88, 90, 92, 93, 94, 101, 103, 110, 111, 117, 118
%peakLats = inf or nan: 20, 79, 80, 82, 86, 88, 90, 92, 99, 100, 110
%straight-up wrong: 101, 103?

%Good estimates:

% indices = find(nPeaksB > 1);
% spec = find(resB > 1);
%FOR B:
%sketchy: 1, 3, 18, 19, 20, 21, 22, 23, 24, 29, 31, 32, 33, 35, 36, 37, 40, 51, 54, 60,
%64, 72, 87, 91, 94, 103, 107
%extra-sketch: 18, 25, 26, 34, 43, 47, 59, 81, 97, 98, 104, 
%peakLats = inf or nan: 27, 28, 33, 39, 43, 47, 80, 94, 103, 105, 107
%straight-up wrong: 64?, 104
%Good estimates: 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 41, 42, 44,
%45, 46, 48, 52, 53, 55, 56, 57, 58, 61, 62, 63, 65-70, 71, 72?, 73, 74,
%75, 76, 79, 85, 86, 88, 90, 92, 93*, 99*, 102, 
%Look at more: 31, 49, 62, 64, 70, 72, 81

indices = find(nPeaksC > 1);
spec = find(resC > 1);

n = indices(spec(1));

subplot(1,2,1)
plot(orbit(n).geolat, orbit(n).F1, '-b')
hold on
plot(orbit(n).geolat(peakInds{n}), orbit(n).F1(peakInds{n}), '*r')
plot(orbit(n).geolat(gradInds{n}), orbit(n).F1(gradInds{n}), 'og')
plot(peakLats(n), peakF1(n), '*k');
hold off

subplot(1,2,2)
plot(orbit(n).geolat, orbit(n).F2, '-b')
hold on
plot(orbit(n).geolat(peakInds{n}), orbit(n).F2(peakInds{n}), '*r')
plot(orbit(n).geolat(qdInds{n}), orbit(n).F2(qdInds{n}), 'oy')
plot(peakLats(n), peakF2(n), '*k');
hold off

shg

%% Test for bias

a = 6371.2;
r = a + 108;

options = optimset('FunValCheck','on');
for i = 1:nUsedA
    if mjdTimesA(i) >= pp.breaks(1) && mjdTimesA(i) <= pp.breaks(end)
        fun = @(theta) findzero(r, theta, pLonsA(i), pp, mjdTimesA(i));
        chaosA_test(i) = fzero(fun, pColatsA(i), options);
    else
        chaosA_test(i) = nan;
    end
end
latsA_test = pLatsA(~isnan(chaosA_test));
chaosA_test = chaosA_test(~isnan(chaosA_test));
colatsA_test = 90 - latsA_test;

for i = 1:nUsedB
    if mjdTimesB(i) >= pp.breaks(1) && mjdTimesB(i) <= pp.breaks(end)
        fun = @(theta) findzero(r, theta, pLonsB(i), pp, mjdTimesB(i));
        chaosB_test(i) = fzero(fun, pColatsB(i), options);
    else
        chaosB_test(i) = nan;
    end
end
latsB_test = pLatsB(~isnan(chaosB_test));
chaosB_test = chaosB_test(~isnan(chaosB_test));
colatsB_test = 90 - latsB_test;

for i = 1:nUsedC
    if mjdTimesC(i) >= pp.breaks(1) && mjdTimesC(i) <= pp.breaks(end)
        fun = @(theta) findzero(r, theta, pLonsC(i), pp, mjdTimesC(i));
        chaosC_test(i) = fzero(fun, pColatsC(i), options);
    else
        chaosC_test(i) = nan;
    end
end
latsC_test = pLatsC(~isnan(chaosC_test));
chaosC_test = chaosC_test(~isnan(chaosC_test));
colatsC_test = 90 - latsC_test;

%% values at EEJ altitude (108 km)
resA_test = colatsA_test - chaosA_test;
resB_test = colatsB_test - chaosB_test;
resC_test = colatsC_test - chaosC_test;
biasA_test = mean(resA_test); % 0.0236 deg
biasB_test = mean(resB_test); % 0.0054 deg
biasC_test = mean(resC_test); % 0.0313 deg
sigmaA_test = std(resA_test); % 0.677 deg
sigmaB_test = std(resB_test); % 0.650 deg
sigmaC_test = std(resC_test); % 0.663 deg
% Altitude bias checks out!
%%

colatsA_minus = pColatsA - 0.7;
colatsA_plus = pColatsA + 0.7;
colatsB_minus = pColatsB - 0.7;
colatsB_plus = pColatsB + 0.7;
colatsC_minus = pColatsC - 0.7;
colatsC_plus = pColatsC + 0.7;

BA_minus = synth_values(radsA, colatsA_minus, lonsA, pp, mjdA);
BB_minus = synth_values(radsB, colatsB_minus, lonsB, pp, mjdB);
BC_minus = synth_values(radsC, colatsC_minus, lonsC, pp, mjdC);
BA_plus = synth_values(radsA, colatsA_plus, lonsA, pp, mjdA);
BB_plus = synth_values(radsB, colatsB_plus, lonsB, pp, mjdB);
BC_plus = synth_values(radsC, colatsC_plus, lonsC, pp, mjdC);

BrA_minus = BA_minus(:,1);
BrA_plus = BA_plus(:,1);
BrB_minus = BB_minus(:,1);
BrB_plus = BB_plus(:,1);
BrC_minus = BC_minus(:,1);
BrC_plus = BC_plus(:,1);

resBrA_minus = BrA - BrA_minus; % mean error: 726 nT
resBrA_plus = BrA - BrA_plus; % mean error: -724 nT
resBrB_minus = BrB - BrB_minus; % mean error: 705 nT
resBrB_plus = BrB - BrB_plus; % mean error: -703 nT
resBrC_minus = BrC - BrC_minus; % mean error: 720 nT
resBrC_plus = BrC - BrC_plus; % mean error: -718 nT
% 700 nT error checks out!

%% Look at orbits with largest residuals to see what is wrong

figure
subplot(3,3,1)
histogram(resA)
title('A')
subplot(3,3,2)
histogram(resB)
title('B')
subplot(3,3,3)
histogram(resC)
title('C')
subplot(3,3,4)
histogram(resA_corr)
title('A corrected')
subplot(3,3,5)
histogram(resB_corr)
title('B corrected')
subplot(3,3,6)
histogram(resC_corr)
title('C corrected')
subplot(3,3,7)
histogram(resA_quiet)
title('A quiet')
subplot(3,3,8)
histogram(resB_quiet)
title('B quiet')
subplot(3,3,9)
histogram(resC_quiet)
title('C quiet')

%%



plot(st(3).orb(quietC(1061)).glats, st(3).orb(quietC(1061)).f1)
hold on
plot(90-colatsC_corr(1061), pF1C(1061), '*g')
plot(90-chaosC(1061), pF1C(1061), 'og')
plot(90-colatsC_corr(quietC(1061)), pF1C(1061), '*b')
plot(90-chaosC(quietC(1061)), pF1C(1061), 'ob')
hold off
legend('Orbit','EEJ corr','CHAOS corr','EEJ quiet','CHAOS quiet','location','best')









