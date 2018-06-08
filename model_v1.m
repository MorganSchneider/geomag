
load('./EEJ_Data/Swarm_Data.mat')

%% Run EEJ algorithm

[peakTimesA, peakLatsA, peakLonsA, peakRadsA, nOrbitsA, nPeaksA] = find_EEJ(swarm, 1, 'i');
nUsedA = length(peakTimesA);
[peakTimesB, peakLatsB, peakLonsB, peakRadsB, nOrbitsB, nPeaksB] = find_EEJ(swarm, 2, 'i');
nUsedB = length(peakTimesB);
% [peakTimesC, peakLatsC, peakLonsC, peakRadsC, nOrbitsC, nPeaksC] = find_EEJ(swarm, 3, 'i');
% nUsedC = length(peakTimesC);

%% Plots of EEJ position


lat_rng = [-90 90];
lon_rng = [-180 180];

figure(1)
ax = worldmap(lat_rng, lon_rng);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow(peakLatsA, peakLonsA, 'DisplayType', 'point')
title('EEJ Peak Swarm A')
%print('-f1', './imgs/swarmA_peak', '-dpng');

figure(2)
ax = worldmap(lat_rng, lon_rng);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow(peakLatsB, peakLonsB, 'DisplayType', 'point')
title('EEJ Peak Swarm B')
%print('-f2', './imgs/swarmB_peak', '-dpng');

figure(3)
ax = worldmap(lat_rng, lon_rng);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow(peakLatsC, peakLonsC, 'DisplayType', 'point')
title('EEJ Peak Swarm C')
%print('-f3', './imgs/swarmC_peak', '-dpng');

figure(4)
subplot(3,1,1)
ax = worldmap(lat_rng, lon_rng);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow(peakLatsA, peakLonsA, 'DisplayType', 'point')
title('EEJ Peak Swarm A')

subplot(3,1,2)
ax = worldmap(lat_rng, lon_rng);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow(peakLatsB, peakLonsB, 'DisplayType', 'point')
title('EEJ Peak Swarm B')

subplot(3,1,3)
ax = worldmap(lat_rng, lon_rng);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow(peakLatsC, peakLonsC, 'DisplayType', 'point')
title('EEJ Peak Swarm C')

%print('-f4', './imgs/combined_peak', '-dpng');

%% Some stuff

% THE TIMES ARE GOOD LEAVE THEM ALONE FOREVER
mjdTimesA = datenum(datetime(peakTimesA,'ConvertFrom','posixtime')) - datenum(2000,1,1,0,0,0);
mjdTimesB = datenum(datetime(peakTimesB,'ConvertFrom','posixtime')) - datenum(2000,1,1,0,0,0);
%mjdTimesC = datenum(datetime(peakTimesC,'ConvertFrom','posixtime')) - datenum(2000,1,1,0,0,0);

% peakColatsA = 90 - peakLatsA;
% peakColatsB = 90 - peakLatsB;
% peakColatsC = 90 - peakLatsC;

%% Load Chaos model coefficients

% load CHAOS model
load('./CHAOS-6_FWD/CHAOS-6-x5.mat')

% low degree field, spline representation up to N
N = 20;   % Take all of core field
pp_N = pp;
pp_N.dim = N*(N+2);
coefs_tmp = reshape(pp.coefs, [], pp.pieces, pp.order);
pp_N.coefs = reshape(coefs_tmp(1:N*(N+2),:,:), [], pp.order);

%% Optimize model for A inputs
theta_init = peakLatsA;
r = peakRadsA;
phi = peakLonsA;
t = mjdTimesA;
options = optimset('MaxFunEvals',1000,'FunValCheck','on');
for i = 1:nUsedA %Need to figure out the iteration problem
    if t(i) >= pp.breaks(1) && t(i) <= pp.breaks(end)
        fun = @(theta) findzero(r(i), theta, phi(i), pp, t(i));
        chaosA(i) = fminsearch(fun, theta_init(i), options);
    else
        chaosA(i) = nan;
    end
end

timesA = peakTimesA(~isnan(chaosA));
latsA = peakLatsA(~isnan(chaosA));
lonsA = peakLonsA(~isnan(chaosA));
radsA = peakRadsA(~isnan(chaosA));
chaosA(isnan(chaosA)) = [];

%% Optimize model for B inputs
theta_init = peakColatsB;
r = peakRadsB;
phi = peakLonsB;
t = mjdTimesB;
for i = 1:nUsedB
    if t(i) >= pp.breaks(1) && t(i) <= pp.breaks(end)
        fun = @(theta) findzero(r(i), theta, phi(i), pp, t(i));
        chaosB(i) = fminsearch(fun, theta_init(i), options);
    else
        chaosB(i) = nan;
    end
end

timesB = peakTimesB(~isnan(chaosB));
latsB = peakLatsB(~isnan(chaosB));
lonsB = peakLonsB(~isnan(chaosB));
radsB = peakRadsB(isnan(chaosB));
chaosB(isnan(chaosB)) = [];

%% Optimize model for C inputs
theta_init = peakColatsC;
r = peakRadsC;
phi = peakLonsC;
t = mjdTimesC;
for i = 1:nUsedC
    if t(i) >= pp.breaks(1) && t(i) <= pp.breaks(end)
        fun = @(theta) findzero(r(i), theta, phi(i), pp, t(i));
        chaosC(i) = fminsearch(fun, theta_init(i), options);
    else
        chaosC(i) = nan;
    end
end

timesC = peakTimesC(~isnan(chaosC));
latsC = peakLatsC(~isnan(chaosC));
lonsC = peakLonsC(~isnan(chaosC));
radsC = peakRadsC(~isnan(chaosC));
chaosC(isnan(chaosC)) = [];


%% Some plots

% Something is seriously wrong with the CHAOS optimized latitude outputs...

resA = latsA - chaosA;
resB = latsB - chaosB;

figure
subplot(2,1,1)
plot(timesA, resA, '.')
subplot(2,1,2)
plot(timesB, resB, '.')

shg

figure
ax = worldmap(lat_rng, lon_rng);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow(chaosA, lonsA, 'DisplayType', 'point')
title('CHAOS-6 Peak')
%print('-f1', './imgs/swarmA_peak', '-dpng');

figure
ax = worldmap(lat_rng, lon_rng);
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
geoshow(chaosB, lonsB, 'DisplayType', 'point')
title('CHAOS-6 Peak')
%print('-f2', './imgs/swarmB_peak', '-dpng');

% figure(3)
% ax = worldmap(lat_rng, lon_rng);
% land = shaperead('landareas', 'UseGeoCoords', true);
% geoshow('landareas.shp', 'FaceColor', [0.5 0.7 0.5])
% geoshow(peakLatsC, peakLonsC, 'DisplayType', 'point')
% title('EEJ Peak Swarm C')
% %print('-f3', './imgs/swarmC_peak', '-dpng');









