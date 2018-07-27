load('./EEJ_Data/Swarm_1HzData.mat')

%% Run EEJ algorithm

[pTimesA, pLatsA, pLonsA, pRadsA, pQDA, nOrbitsA, nPeaksA] = find_EEJ(swarm, 1);
[pTimesB, pLatsB, pLonsB, pRadsB, pQDB, nOrbitsB, nPeaksB] = find_EEJ(swarm, 2);
[pTimesC, pLatsC, pLonsC, pRadsC, pQDC, nOrbitsC, nPeaksC] = find_EEJ(swarm, 3);


%% Plots of EEJ position


lat_rng = [-90 90];
lon_rng = [-180 180];

figure(1)
ax = worldmap(lat_rng, lon_rng);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow(pLatsA, pLonsA, 'DisplayType', 'point')
title('EEJ Peak Swarm A')

figure(2)
ax = worldmap(lat_rng, lon_rng);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow(pLatsB, pLonsB, 'DisplayType', 'point')
title('EEJ Peak Swarm B')

figure(3)
ax = worldmap(lat_rng, lon_rng);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow(pLatsC, pLonsC, 'DisplayType', 'point')
title('EEJ Peak Swarm C')


figure(4)
subplot(1,3,1)
ax = worldmap(lat_rng, lon_rng);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow(pLatsA, pLonsA, 'DisplayType', 'point')
title('EEJ Peak Swarm A')

subplot(1,3,2)
ax = worldmap(lat_rng, lon_rng);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow(pLatsB, pLonsB, 'DisplayType', 'point')
title('EEJ Peak Swarm B')

subplot(1,3,3)
ax = worldmap(lat_rng, lon_rng);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow(pLatsC, pLonsC, 'DisplayType', 'point')
title('EEJ Peak Swarm C')

%%
mjdTimesA = datenum(datetime(pTimesA,'ConvertFrom','posixtime')) - datenum(2000,1,1,0,0,0);
mjdTimesB = datenum(datetime(pTimesB,'ConvertFrom','posixtime')) - datenum(2000,1,1,0,0,0);
mjdTimesC = datenum(datetime(pTimesC,'ConvertFrom','posixtime')) - datenum(2000,1,1,0,0,0);

pColatsA = 90 - pLatsA;
pColatsB = 90 - pLatsB;
pColatsC = 90 - pLatsC;

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
for i = 1:length(pTimesA)
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
qdA = pQDA(~isnan(chaosA));
chaosA = chaosA(~isnan(chaosA));

%% Optimize model around B inputs
for i = 1:length(pTimesB)
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
chaosB = chaosB(~isnan(chaosB));

%% Optimize model around C inputs
for i = 1:length(pTimesC)
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
qdC = pQDC(~isnan(chaosC));
chaosC = chaosC(~isnan(chaosC));

%%

colatsA = 90 - latsA;
colatsB = 90 - latsB;
colatsC = 90 - latsC;

resA = colatsA - chaosA;
resB = colatsB - chaosB;
resC = colatsC - chaosC;

lonBins = -180:10:170;
biasA = zeros(1, length(lonBins));
biasB = zeros(1, length(lonBins));
biasC = zeros(1, length(lonBins));
sigmaA = zeros(1, length(lonBins));
sigmaB = zeros(1, length(lonBins));
sigmaC = zeros(1, length(lonBins));
j = 1;
for i = lonBins
    indsA = find(lonsA >= i & lonsA < i+10);
    indsB = find(lonsB >= i & lonsB < i+10);
    indsC = find(lonsC >= i & lonsC < i+10);
    biasA(j) = mean(resA(indsA));
    biasB(j) = mean(resB(indsB));
    biasC(j) = mean(resC(indsC));
    sigmaA(j) = std(resA(indsA));
    sigmaB(j) = std(resB(indsB));
    sigmaC(j) = std(resC(indsC));
    j = j + 1;
end
midLons = -175:10:175;

figure(5)

subplot(3,1,1)
plot(midLons, biasA, '-')
hold on
plot(midLons, [biasA + sigmaA], '*')
plot(midLons, [biasA - sigmaA], '*')
hold off
title('Swarm A Equator Residuals (deg)')
xlabel('Longitude (deg)')

subplot(3,1,2)
plot(midLons, biasB, '-')
hold on
plot(midLons, biasB + sigmaB, '*')
plot(midLons, biasB - sigmaB, '*')
hold off
title('Swarm B Equator Residuals (deg)')
xlabel('Longitude (deg)')

subplot(3,1,3)
plot(midLons, biasC, '-')
hold on
plot(midLons, biasC + sigmaC, '*')
plot(midLons, biasC - sigmaC, '*')
hold off
title('Swarm C Equator Residuals (deg)')
xlabel('Longitude (deg)')


figure(6)
subplot(3,1,1)
plot(lonsA, resA, '*')
title('Swarm A Equator Residuals (deg)')
xlabel('Longitude (deg)')
subplot(3,1,2)
plot(lonsB, resB, '*')
title('Swarm B Equator Residuals (deg)')
xlabel('Longitude (deg)')
subplot(3,1,3)
plot(lonsC, resC, '*')
title('Swarm C Equator Residuals (deg)')
xlabel('Longitude (deg)')






