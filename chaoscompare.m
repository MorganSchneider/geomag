
load('./EEJ_Data/Swarm_1HzData.mat')

%% Run EEJ algorithm

[peakTimesA, peakLatsA, peakLonsA, peakRadsA, peakQDA, nOrbitsA, nPeaksA] = find_EEJ(swarm, 1);
[peakTimesB, peakLatsB, peakLonsB, peakRadsB, peakQDB, nOrbitsB, nPeaksB] = find_EEJ(swarm, 2);
[peakTimesC, peakLatsC, peakLonsC, peakRadsC, peakQDC, nOrbitsC, nPeaksC] = find_EEJ(swarm, 3);

%% Remove orbits with more than one detected peak

orbitA.orig = find(nPeaksA == 1);
orbitB.orig = find(nPeaksB == 1);
orbitC.orig = find(nPeaksC == 1);

nPeaksA(nPeaksA == 0) = [];
pTimesA = peakTimesA(nPeaksA == 1);
pLatsA = peakLatsA(nPeaksA == 1);
pLonsA = peakLonsA(nPeaksA == 1);
pRadsA = peakRadsA(nPeaksA == 1);
pQDA = peakQDA(nPeaksA == 1);
% pF1A = peakF1A(nPeaksA == 1); %%%% just for debugging
% pF2A = peakF2A(nPeaksA == 1); %%%% just for debugging

nPeaksB(nPeaksB == 0) = [];
pTimesB = peakTimesB(nPeaksB == 1);
pLatsB = peakLatsB(nPeaksB == 1);
pLonsB = peakLonsB(nPeaksB == 1);
pRadsB = peakRadsB(nPeaksB == 1);
pQDB = peakQDB(nPeaksB == 1);
% pF1B = peakF1B(nPeaksB == 1); %%%% just for debugging
% pF2B = peakF2B(nPeaksB == 1); %%%% just for debugging

nPeaksC(nPeaksC == 0) = [];
pTimesC = peakTimesC(nPeaksC == 1);
pLatsC = peakLatsC(nPeaksC == 1);
pLonsC = peakLonsC(nPeaksC == 1);
pRadsC = peakRadsC(nPeaksC == 1);
pQDC = peakQDC(nPeaksC == 1);
% pF1C = peakF1C(nPeaksC == 1); %%%% just for debugging
% pF2C = peakF2C(nPeaksC == 1); %%%% just for debugging
% pQdC = peakQdC(nPeaksC == 1); %%%% just for debugging

nUsedA = length(pTimesA); % 4387
nUsedB = length(pTimesB); % 3117
nUsedC = length(pTimesC); % 2814

%% Plots of EEJ position


lat_rng = [-90 90];
lon_rng = [-180 180];

figure(1)
ax = worldmap(lat_rng, lon_rng);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow(pLatsA, pLonsA, 'DisplayType', 'point')
title('EEJ Peak Swarm A')
imgName = sprintf('./imgs/%s/swarmA_peak', method);
print('-f1', imgName, '-dpng');

figure(2)
ax = worldmap(lat_rng, lon_rng);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow(pLatsB, pLonsB, 'DisplayType', 'point')
title('EEJ Peak Swarm B')
imgName = sprintf('./imgs/%s/swarmB_peak', method);
print('-f2', imgName, '-dpng');

figure(3)
ax = worldmap(lat_rng, lon_rng);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow(pLatsC, pLonsC, 'DisplayType', 'point')
title('EEJ Peak Swarm C')
imgName = sprintf('./imgs/%s/swarmC_peak', method);
print('-f3', imgName, '-dpng');

figure(4)
subplot(3,1,1)
ax = worldmap(lat_rng, lon_rng);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow(pLatsA, pLonsA, 'DisplayType', 'point')
title('EEJ Peak Swarm A')

subplot(3,1,2)
ax = worldmap(lat_rng, lon_rng);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow(pLatsB, pLonsB, 'DisplayType', 'point')
title('EEJ Peak Swarm B')

subplot(3,1,3)
ax = worldmap(lat_rng, lon_rng);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow(pLatsC, pLonsC, 'DisplayType', 'point')
title('EEJ Peak Swarm C')

imgName = sprintf('./imgs/%s/combined_peak', method);
print('-f4', imgName, '-dpng');

%% Conversions for synth_values

mjdTimesA = datenum(datetime(pTimesA,'ConvertFrom','posixtime')) - datenum(2000,1,1,0,0,0);
mjdTimesB = datenum(datetime(pTimesB,'ConvertFrom','posixtime')) - datenum(2000,1,1,0,0,0);
mjdTimesC = datenum(datetime(pTimesC,'ConvertFrom','posixtime')) - datenum(2000,1,1,0,0,0);

pColatsA = 90 - pLatsA;
pColatsB = 90 - pLatsB;
pColatsC = 90 - pLatsC;

%% Load Chaos model coefficients

% load CHAOS model
load('./CHAOS-6_FWD/CHAOS-6-x5.mat')

% low degree field, spline representation up to N
N = 10;   % Take all of core field
pp_N = pp;
pp_N.dim = N*(N+2);
coefs_tmp = reshape(pp.coefs, [], pp.pieces, pp.order);
pp_N.coefs = reshape(coefs_tmp(1:N*(N+2),:,:), [], pp.order);

%% Optimize model around A inputs
options = optimset('FunValCheck','on');
for i = 1:nUsedA
    if mjdTimesA(i) >= pp.breaks(1) && mjdTimesA(i) <= pp.breaks(end)
        fun = @(theta) findzero(pRadsA(i), theta, pLonsA(i), pp, mjdTimesA(i));
        chaosA(i) = fzero(fun, pColatsA(i), options);
    else
        chaosA(i) = nan;
    end
end

orbitA.chaos = orbitA.orig(~isnan(chaosA));
timesA = pTimesA(~isnan(chaosA));
latsA = pLatsA(~isnan(chaosA));
lonsA = pLonsA(~isnan(chaosA));
radsA = pRadsA(~isnan(chaosA));
qdA = pQDA(~isnan(chaosA));
% f1A = pF1A(~isnan(chaosA)); %%%% just for debugging
% f2A = pF2A(~isnan(chaosA)); %%%% just for debugging
chaosA = chaosA(~isnan(chaosA));

%% Optimize model around B inputs
for i = 1:nUsedB
    if mjdTimesB(i) >= pp.breaks(1) && mjdTimesB(i) <= pp.breaks(end)
        fun = @(theta) findzero(pRadsB(i), theta, pLonsB(i), pp, mjdTimesB(i));
        chaosB(i) = fzero(fun, pColatsB(i), options);
    else
        chaosB(i) = nan;
    end
end

orbitB.chaos = orbitB.orig(~isnan(chaosB));
timesB = pTimesB(~isnan(chaosB));
latsB = pLatsB(~isnan(chaosB));
lonsB = pLonsB(~isnan(chaosB));
radsB = pRadsB(~isnan(chaosB));
qdB = pQDB(~isnan(chaosB));
% f1B = pF1B(~isnan(chaosB)); %%%% just for debugging
% f2B = pF2B(~isnan(chaosB)); %%%% just for debugging
chaosB = chaosB(~isnan(chaosB));

%% Optimize model around C inputs
for i = 1:nUsedC
    if mjdTimesC(i) >= pp.breaks(1) && mjdTimesC(i) <= pp.breaks(end)
        fun = @(theta) findzero(pRadsC(i), theta, pLonsC(i), pp, mjdTimesC(i));
        chaosC(i) = fzero(fun, pColatsC(i), options);
    else
        chaosC(i) = nan;
    end
end

orbitC.chaos = orbitC.orig(~isnan(chaosC));
timesC = pTimesC(~isnan(chaosC));
latsC = pLatsC(~isnan(chaosC));
lonsC = pLonsC(~isnan(chaosC));
radsC = pRadsC(~isnan(chaosC));
qdC = pQDC(~isnan(chaosC));
% f1C = pF1C(~isnan(chaosC)); %%%% just for debugging
% f2C = pF2C(~isnan(chaosC)); %%%% just for debugging
chaosC = chaosC(~isnan(chaosC));


%% Stuff

colatsA = 90 - latsA;
colatsB = 90 - latsB;
colatsC = 90 - latsC;

resA = colatsA - chaosA;
resB = colatsB - chaosB;
resC = colatsC - chaosC;


sigmaA_total = std(resA);
sigmaB_total = std(resB);
sigmaC_total = std(resC);
biasA_total = nanmean(resA);
biasB_total = nanmean(resB);
biasC_total = nanmean(resC);
varA_total = var(resA);
varB_total = var(resB);
varC_total = var(resC);

%% Some plots

figure(5)
subplot(3,1,1)
plot(timesA, resA, '.')
ax1 = gca;
ax1.YLim = [-20 20];
title('{\theta}_{EEJ_A} - {\theta}_{CHAOS}')
subplot(3,1,2)
plot(timesB, resB, '.')
ax2 = gca;
ax2.YLim = [-20 20];
title('{\theta}_{EEJ_B} - {\theta}_{CHAOS}')
subplot(3,1,3)
plot(timesC, resC, '.')
ax3 = gca;
ax3.XLim = ax1.XLim;
ax3.YLim = [-20 20];
title('{\theta}_{EEJ_C} - {\theta}_{CHAOS}')
xlabel('Time (s)')
imgName = sprintf('./imgs/%s/residuals_time', method);
print('-f5', imgName, '-dpng');

figure(6)
subplot(3,1,1)
plot(lonsA, resA, '.')
ax1 = gca;
ax1.XLim = [-200 200];
ax1.YLim = [-20 20];
title('{\theta}_{EEJ_A} - {\theta}_{CHAOS}')
subplot(3,1,2)
plot(lonsB, resB, '.')
ax2 = gca;
ax2.XLim = [-200 200];
ax2.YLim = [-20 20];
title('{\theta}_{EEJ_B} - {\theta}_{CHAOS}')
subplot(3,1,3)
plot(lonsC, resC, '.')
ax3 = gca;
ax3.XLim = [-200 200];
ax3.YLim = [-20 20];
title('{\theta}_{EEJ_C} - {\theta}_{CHAOS}')
xlabel('Longitude (degrees)')
imgName = sprintf('./imgs/%s/residuals_lon', method);
print('-f6', imgName, '-dpng');

figure(7)
subplot(3,1,1)
plot(qdA, resA, '.')
ax1 = gca;
ax1.XLim = [5 20];
ax1.YLim = [-20 20];
title('{\theta}_{EEJ_A} - {\theta}_{CHAOS}')
subplot(3,1,2)
plot(qdB, resB, '.')
ax2 = gca;
ax2.XLim = [5 20];
ax2.YLim = [-20 20];
title('{\theta}_{EEJ_B} - {\theta}_{CHAOS}')
subplot(3,1,3)
plot(qdC, resC, '.')
ax3 = gca;
ax3.XLim = [5 20];
ax3.YLim = [-20 20];
title('{\theta}_{EEJ_C} - {\theta}_{CHAOS}')
xlabel('Local Time (UTC)')
imgName = sprintf('./imgs/%s/residuals_local', method);
print('-f7', imgName, '-dpng');

figure(8)
ax = worldmap(lat_rng, lon_rng);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow([90-chaosA, 90-chaosB, 90-chaosC], [lonsA, lonsB, lonsC], 'DisplayType', 'point')
title('CHAOS-6 Peak')
print('-f8', './imgs/chaos_peak', '-dpng');



%% Correct for altitude bias

colatsA_corr = colatsA + 0.1; %correction for -0.1 degree bias
colatsB_corr = colatsB + 0.1;
colatsC_corr = colatsC + 0.1;

resA_corr = colatsA_corr - chaosA;
resB_corr = colatsB_corr - chaosB;
resC_corr = colatsC_corr - chaosC;

sigmaA_corr = std(resA_corr);
sigmaB_corr = std(resB_corr);
sigmaC_corr = std(resC_corr);
biasA_corr = nanmean(resA_corr);
biasB_corr = nanmean(resB_corr);
biasC_corr = nanmean(resC_corr);
varA_corr = var(resA_corr);
varB_corr = var(resB_corr);
varC_corr = var(resC_corr);

%% Introduce geomagnetic indices as another filter

[timestamp, kp, rc] = read_indices(swarm);

kpA = zeros(1, length(timesA));
rcA = zeros(1, length(timesA));
for i = 1:length(timesA)
    z = find(timestamp <= timesA(i));
    ind = z(end);
    kpA(i) = kp(ind);
    rcA(i) = rc(ind);
end
kpB = zeros(1, length(timesB));
rcB = zeros(1, length(timesB));
for i = 1:length(timesB)
    z = find(timestamp <= timesB(i));
    ind = z(end);
    kpB(i) = kp(ind);
    rcB(i) = rc(ind);
end
kpC = zeros(1, length(timesC));
rcC = zeros(1, length(timesC));
for i = 1:length(timesC)
    z = find(timestamp <= timesC(i));
    ind = z(end);
    kpC(i) = kp(ind);
    rcC(i) = rc(ind);
end

%% Filter by indices

x1 = find(kpA <= 2);
x2 = find(abs(rcA) <= 15);
lia = ismember(x1, x2);
quietA = x1(lia == 1);

x1 = find(kpB <= 2);
x2 = find(abs(rcB) <= 15);
lia = ismember(x1, x2);
quietB = x1(lia == 1);

x1 = find(kpC <= 2);
x2 = find(abs(rcC) <= 15);
lia = ismember(x1, x2);
quietC = x1(lia == 1);

orbitA.quiet = orbitA.chaos(quietA);
orbitB.quiet = orbitB.chaos(quietB);
orbitC.quiet = orbitC.chaos(quietC);

colatsA_quiet = colatsA_corr(quietA);
colatsB_quiet = colatsB_corr(quietB);
colatsC_quiet = colatsC_corr(quietC);

resA_quiet = colatsA_corr(quietA) - chaosA(quietA);
resB_quiet = colatsB_corr(quietB) - chaosB(quietB);
resC_quiet = colatsC_corr(quietC) - chaosC(quietC);

sigmaA_quiet = std(resA_quiet);
sigmaB_quiet = std(resB_quiet);
sigmaC_quiet = std(resC_quiet);
biasA_quiet = nanmean(resA_quiet);
biasB_quiet = nanmean(resB_quiet);
biasC_quiet = nanmean(resC_quiet);
varA_quiet = var(resA_quiet);
varB_quiet = var(resB_quiet);
varC_quiet = var(resC_quiet);

% boxplots of residuals by longitude

%% Testing electrojet peak positions for vertical field values

mjdA = datenum(datetime(timesA,'ConvertFrom','posixtime')) - datenum(2000,1,1,0,0,0);
mjdB = datenum(datetime(timesB,'ConvertFrom','posixtime')) - datenum(2000,1,1,0,0,0);
mjdC = datenum(datetime(timesC,'ConvertFrom','posixtime')) - datenum(2000,1,1,0,0,0);

BA = synth_values(radsA, colatsA, lonsA, pp, mjdA);
BB = synth_values(radsB, colatsB, lonsB, pp, mjdB);
BC = synth_values(radsC, colatsC, lonsC, pp, mjdC);
BA_corr = synth_values(radsA, colatsA_corr, lonsA, pp, mjdA);
BB_corr = synth_values(radsB, colatsB_corr, lonsB, pp, mjdB);
BC_corr = synth_values(radsC, colatsC_corr, lonsC, pp, mjdC);
BA_quiet = synth_values(radsA(quietA), colatsA_corr(quietA), lonsA(quietA),...
    pp, mjdA(quietA));
BB_quiet = synth_values(radsB(quietB), colatsB_corr(quietB), lonsB(quietB),...
    pp, mjdB(quietB));
BC_quiet = synth_values(radsC(quietC), colatsC_corr(quietC), lonsC(quietC),...
    pp, mjdC(quietC));

BrA = BA(:,1);
BrB = BB(:,1);
BrC = BC(:,1);
BrA_corr = BA_corr(:,1);
BrB_corr = BB_corr(:,1);
BrC_corr = BC_corr(:,1);
BrA_quiet = BA_quiet(:,1);
BrB_quiet = BB_quiet(:,1);
BrC_quiet = BC_quiet(:,1);

FA = sqrt(BA(:,1).^2 + BA(:,2).^2 + BA(:,3).^2);
FB = sqrt(BB(:,1).^2 + BB(:,2).^2 + BB(:,3).^2);
FC = sqrt(BC(:,1).^2 + BC(:,2).^2 + BC(:,3).^2);
FA_corr = sqrt(BA_corr(:,1).^2 + BA_corr(:,2).^2 + BA_corr(:,3).^2);
FB_corr = sqrt(BB_corr(:,1).^2 + BB_corr(:,2).^2 + BB_corr(:,3).^2);
FC_corr = sqrt(BC_corr(:,1).^2 + BC_corr(:,2).^2 + BC_corr(:,3).^2);
FA_quiet = sqrt(BA_quiet(:,1).^2 + BA_quiet(:,2).^2 + BA_quiet(:,3).^2);
FB_quiet = sqrt(BB_quiet(:,1).^2 + BB_quiet(:,2).^2 + BB_quiet(:,3).^2);
FC_quiet = sqrt(BC_quiet(:,1).^2 + BC_quiet(:,2).^2 + BC_quiet(:,3).^2);


%% Comparisons (no filter, bias-corrected, quiet-corrected)

meanBr(1,1) = nanmean(BrA);
meanBr(1,2) = nanmean(BrB);
meanBr(1,3) = nanmean(BrC);
meanBr(2,1) = nanmean(BrA_corr);
meanBr(2,2) = nanmean(BrB_corr);
meanBr(2,3) = nanmean(BrC_corr);
meanBr(3,1) = nanmean(BrA_quiet);
meanBr(3,2) = nanmean(BrB_quiet);
meanBr(3,3) = nanmean(BrC_quiet);

stdBr(1,1) = std(BrA);
stdBr(1,2) = std(BrB);
stdBr(1,3) = std(BrC);
stdBr(2,1) = std(BrA_corr);
stdBr(2,2) = std(BrB_corr);
stdBr(2,3) = std(BrC_corr);
stdBr(3,1) = std(BrA_quiet);
stdBr(3,2) = std(BrB_quiet);
stdBr(3,3) = std(BrC_quiet);

varBr(1,1) = var(BrA);
varBr(1,2) = var(BrB);
varBr(1,3) = var(BrC);
varBr(2,1) = var(BrA_corr);
varBr(2,2) = var(BrB_corr);
varBr(2,3) = var(BrC_corr);
varBr(3,1) = var(BrA_quiet);
varBr(3,2) = var(BrB_quiet);
varBr(3,3) = var(BrC_quiet);



meanRes(1,1) = biasA_total;
meanRes(1,2) = biasB_total;
meanRes(1,3) = biasC_total;
meanRes(2,1) = biasA_corr;
meanRes(2,2) = biasB_corr;
meanRes(2,3) = biasC_corr;
meanRes(3,1) = biasA_quiet;
meanRes(3,2) = biasB_quiet;
meanRes(3,3) = biasC_quiet;

stdRes(1,1) = sigmaA_total;
stdRes(1,2) = sigmaB_total;
stdRes(1,3) = sigmaC_total;
stdRes(2,1) = sigmaA_corr;
stdRes(2,2) = sigmaB_corr;
stdRes(2,3) = sigmaC_corr;
stdRes(3,1) = sigmaA_quiet;
stdRes(3,2) = sigmaB_quiet;
stdRes(3,3) = sigmaC_quiet;

varRes(1,1) = varA_total;
varRes(1,2) = varB_total;
varRes(1,3) = varC_total;
varRes(2,1) = varA_corr;
varRes(2,2) = varB_corr;
varRes(2,3) = varC_corr;
varRes(3,1) = varA_quiet;
varRes(3,2) = varB_quiet;
varRes(3,3) = varC_quiet;










