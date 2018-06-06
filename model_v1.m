
load('./EEJ_Data/Swarm_Data.mat')

%% Run EEJ algorithm

[peakTimesA, peakLatsA, peakLonsA, peakRadsA, nOrbitsA, nPeaksA] = find_EEJ(swarm, 1, 'm');
nUsedA = length(peakTimesA);
[peakTimesB, peakLatsB, peakLonsB, peakRadsB, nOrbitsB, nPeaksB] = find_EEJ(swarm, 2, 'm');
nUsedB = length(peakTimesB);
[peakTimesC, peakLatsC, peakLonsC, peakRadsC, nOrbitsC, nPeaksC] = find_EEJ(swarm, 3, 'm');
nUsedC = length(peakTimesC);

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
mjdTimesC = datenum(datetime(peakTimesC,'ConvertFrom','posixtime')) - datenum(2000,1,1,0,0,0);

peakRadsA_km = peakRadsA / 1000;
peakRadsB_km = peakRadsB / 1000;
peakRadsC_km = peakRadsC / 1000;

peakColatsA = 90 - peakLatsA;
peakColatsB = 90 - peakLatsB;
peakColatsC = 90 - peakLatsC;

%% Load Chaos model coefficients

chaos = struct();
% load CHAOS model
load('./CHAOS-6_FWD/CHAOS-6-x5.mat')

% low degree field, spline representation up to N
N = 20;   % Take all of core field
pp_N = pp;
pp_N.dim = N*(N+2);
coefs_tmp = reshape(pp.coefs, [], pp.pieces, pp.order);
pp_N.coefs = reshape(coefs_tmp(1:N*(N+2),:,:), [], pp.order);

%% Optimize model for A inputs
theta_init = peakColatsA;
r = peakRadsA_km;
phi = peakLonsA;
t = mjdTimesA;
fun = @(theta) synth_values_playbox(r, theta, phi, pp, t);
vChaosA = fminfind(fun, theta_init);
%% Optimize model for B inputs
theta_init = peakColatsB;
r = peakRadsB_km;
phi = peakLonsB;
t = mjdTimesB;
fun = @(theta) synth_values_playbox(r, theta, phi, pp, t);
vChaosB = fminfind(fun, theta_init);
%% Optimize model for C inputs
theta_init = peakColatsC;
r = peakRadsC_km;
phi = peakLonsC;
t = mjdTimesC;
fun = @(theta) synth_values_playbox(r, theta, phi, pp, t);
vChaosC = fminfind(fun, theta_init);





%% Some stuff, continued

% Find mean and standard deviation of latitude for every 10(?) degrees
% longitude
% Need Swarm vector data to compare equator positions









