
load('./EEJ_Data/Swarm_1HzData.mat')

%% Run EEJ algorithm

method = 'interp';

[peakTimesA, peakLatsA, peakLonsA, peakRadsA, peakLocalA, nOrbitsA, nPeaksA] = find_EEJ(swarm, 1, method);
[peakTimesB, peakLatsB, peakLonsB, peakRadsB, peakLocalB, nOrbitsB, nPeaksB] = find_EEJ(swarm, 2, method);
[peakTimesC, peakLatsC, peakLonsC, peakRadsC, peakLocalC, nOrbitsC, nPeaksC] = find_EEJ(swarm, 3, method);

%% Remove orbits with more than one detected peak

nPeaksA(nPeaksA == 0) = [];
pTimesA = peakTimesA(nPeaksA == 1);
pLatsA = peakLatsA(nPeaksA == 1);
pLonsA = peakLonsA(nPeaksA == 1);
pRadsA = peakRadsA(nPeaksA == 1);
pLocalA = peakLocalA(nPeaksA == 1);

nPeaksB(nPeaksB == 0) = [];
pTimesB = peakTimesB(nPeaksB == 1);
pLatsB = peakLatsB(nPeaksB == 1);
pLonsB = peakLonsB(nPeaksB == 1);
pRadsB = peakRadsB(nPeaksB == 1);
pLocalB = peakLocalB(nPeaksB == 1);

nPeaksC(nPeaksC == 0) = [];
pTimesC = peakTimesC(nPeaksC == 1);
pLatsC = peakLatsC(nPeaksC == 1);
pLonsC = peakLonsC(nPeaksC == 1);
pRadsC = peakRadsC(nPeaksC == 1);
pLocalC = peakLocalC(nPeaksC == 1);

nUsedA = length(pTimesA);
nUsedB = length(pTimesB);
nUsedC = length(pTimesC);

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

% THE TIMES ARE GOOD LEAVE THEM ALONE FOREVER
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
N = 20;   % Take all of core field
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

timesA = pTimesA(~isnan(chaosA));
latsA = pLatsA(~isnan(chaosA));
lonsA = pLonsA(~isnan(chaosA));
radsA = pRadsA(~isnan(chaosA));
localA = pLocalA(~isnan(chaosA));
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

timesB = pTimesB(~isnan(chaosB));
latsB = pLatsB(~isnan(chaosB));
lonsB = pLonsB(~isnan(chaosB));
radsB = pRadsB(~isnan(chaosB));
localB = pLocalB(~isnan(chaosB));
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

timesC = pTimesC(~isnan(chaosC));
latsC = pLatsC(~isnan(chaosC));
lonsC = pLonsC(~isnan(chaosC));
radsC = pRadsC(~isnan(chaosC));
localC = pLocalC(~isnan(chaosC));
chaosC = chaosC(~isnan(chaosC));


%% Stuff

colatsA = 90 - latsA;
colatsB = 90 - latsB;
colatsC = 90 - latsC;

resA = colatsA - chaosA;
resB = colatsB - chaosB;
resC = colatsC - chaosC;

%look at grouping these by month (set vector length), time of day (use
%local), to get a better idea of variability with the model
sigmaA_total = std(resA);
sigmaB_total = std(resB);
sigmaC_total = std(resC);
biasA_total = mean(resA);
biasB_total = mean(resB);
biasC_total = mean(resC);

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
plot(localA, resA, '.')
ax1 = gca;
ax1.XLim = [5 20];
ax1.YLim = [-20 20];
title('{\theta}_{EEJ_A} - {\theta}_{CHAOS}')
subplot(3,1,2)
plot(localB, resB, '.')
ax2 = gca;
ax2.XLim = [5 20];
ax2.YLim = [-20 20];
title('{\theta}_{EEJ_B} - {\theta}_{CHAOS}')
subplot(3,1,3)
plot(localC, resC, '.')
ax3 = gca;
ax3.XLim = [5 20];
ax3.YLim = [-20 20];
title('{\theta}_{EEJ_C} - {\theta}_{CHAOS}')
xlabel('Local Time (UTC)')
imgName = sprintf('./imgs/%s/residuals_local', method);
print('-f7', imgName, '-dpng');

%%
figure(8)
ax = worldmap(lat_rng, lon_rng);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow([90-chaosA, 90-chaosB, 90-chaosC], [lonsA, lonsB, lonsC], 'DisplayType', 'point')
title('CHAOS-6 Peak')
print('-f8', './imgs/chaos_peak', '-dpng');

%% See which local times have the least error (if significant)

sigmaA = zeros(1,24);
biasA = zeros(1,24);
for i = ceil(min(localA)): floor(max(localA))
    temp = localA(localA < i+1);
    inds = find(temp >= i);
    sigmaA(i) = std(resA(inds));
    biasA(i) = mean(resA(inds));
end
sigmaA(sigmaA == 0) = nan;
biasA(biasA == 0) = nan;

sigmaB = zeros(1,24);
biasB = zeros(1,24);
for i = ceil(min(localB)): floor(max(localB))
    temp = localB(localB < i+1);
    inds = find(temp >= i);
    sigmaB(i) = std(resB(inds));
    biasB(i) = mean(resB(inds));
end
sigmaB(sigmaB == 0) = nan;
biasB(biasB == 0) = nan;

sigmaC = zeros(1,24);
biasC = zeros(1,24);
for i = ceil(min(localC)): floor(max(localC))
    temp = localC(localC < i+1);
    inds = find(temp >= i);
    sigmaC(i) = std(resC(inds));
    biasC(i) = mean(resC(inds));
end
sigmaC(sigmaC == 0) = nan;
biasC(biasC == 0) = nan;

[~, bestLocalA] = sort(sigmaA);
[~, bestLocalB] = sort(sigmaB);
[~, bestLocalC] = sort(sigmaC);

[~, leastBiasA] = sort(biasA);
[~, leastBiasB] = sort(biasB);
[~, leastBiasC] = sort(biasC);


colatsA_corr = colatsA + 0.1; %correction for -0.1 degree bias
colatsB_corr = colatsB + 0.1;
colatsC_corr = colatsC + 0.1;

resA_corr = colatsA_corr - chaosA;
resB_corr = colatsB_corr - chaosB;
resC_corr = colatsC_corr - chaosC;

%% Introduce geomagnetic indices as another filter

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

x1 = find(kpA <= 3);
x2 = find(abs(rcA) <= 3);
lia = ismember(x1, x2);
quietA = x1(lia == 1);

x1 = find(kpB <= 3);
x2 = find(abs(rcB) <= 3);
lia = ismember(x1, x2);
quietB = x1(lia == 1);

x1 = find(kpC <= 3);
x2 = find(abs(rcC) <= 3);
lia = ismember(x1, x2);
quietC = x1(lia == 1);

%% Filter by indices

resA_quiet = colatsA_corr(quietA) - chaosA(quietA);
resB_quiet = colatsB_corr(quietB) - chaosB(quietB);
resC_quiet = colatsC_corr(quietC) - chaosC(quietC);

sigmaA_quiet = std(resA_quiet);
sigmaB_quiet = std(resB_quiet);
sigmaC_quiet = std(resC_quiet);
biasA_quiet = mean(resA_quiet);
biasB_quiet = mean(resB_quiet);
biasC_quiet = mean(resC_quiet);












