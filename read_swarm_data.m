% # Field 1: timestamp (seconds since 1970-01-01 00:00:00 UTC)
% # Field 2: local time (hours)
% # Field 3: season (day of year in [0,365])
% # Field 4: geocentric radius (km)
% # Field 5: longitude (degrees)
% # Field 6: geocentric latitude (degrees)
% # Field 7: QD latitude (degrees)
% # Field 8: F_sat (nT)
% # Field 9: F_internal (nT)
% # Field 10: dF_external (nT)
% # Field 11: F^(1) (nT)
% # Field 12: internal Sq model (nT)
% # Field 13: external Sq model (nT)
% # Field 14: F^(2) = F^(1) - Sq_model (nT)
% # Field 15: F^(2) fit from EEJ model (nT)

% Load A
A = cell(1,15);
for i = 1:15
    fname = sprintf('./EEJ_Data/SwarmA/SwarmA_1Hzcol%d.txt', i);
    fp = fopen(fname);
    a = textscan(fp, '%f', 'HeaderLines', 15, 'delimiter', '\r\n'); A{i} = a{1};
    fclose(fp);
end

% Load B
B = cell(1,15);
for i = 1:15
    fname = sprintf('./EEJ_Data/SwarmB/SwarmB_1Hzcol%d.txt', i);
    fp = fopen(fname);
    b = textscan(fp, '%f', 'HeaderLines', 15, 'delimiter', '\r\n'); B{i} = b{1};
    fclose(fp);
end

% Load C
C = cell(1,15);
for i = 1:15
   fname = sprintf('./EEJ_Data/SwarmC/SwarmC_1Hzcol%d.txt', i);
   fp = fopen(fname);
   c = textscan(fp, '%f', 'HeaderLines', 15, 'delimiter', '\r\n'); C{i} = c{1};
   fclose(fp);
end

% Organize into structure array and save to .mat file


swarm = struct('time', [], 'local', [], 'season', [], 'rad', [], 'lon', [], 'geolat', [],...
    'qdlat', [], 'F1', [], 'F2', []);

TIMES = {A{1}, B{1}, C{1}};
LOCAL = {A{2}, B{2}, C{2}};
SEASON = {A{3}, B{3}, C{3}};
RADS = {A{4}, B{4}, C{4}};
LONS = {A{5}, B{5}, C{5}};
GEOLATS = {A{6}, B{6}, C{6}};
QDLATS = {A{7}, B{7}, C{7}};
F1S = {A{11}, B{11}, C{11}};
F2S = {A{14}, B{14}, C{14}};

for i = 1:3
    swarm(i).time = rot90(TIMES{i});
    swarm(i).time(isnan(swarm(i).time)) = [];
    
    swarm(i).local = rot90(LOCAL{i});
    swarm(i).local(isnan(swarm(i).local)) = [];
    
    swarm(i).season = rot90(SEASON{i});
    swarm(i).season(isnan(swarm(i).season)) = [];
    
    swarm(i).rad = rot90(RADS{i});
    swarm(i).rad(isnan(swarm(i).rad)) = [];
    
    swarm(i).lon = rot90(LONS{i});
    swarm(i).lon(isnan(swarm(i).lon)) = [];
    
    swarm(i).geolat = rot90(GEOLATS{i});
    swarm(i).geolat(isnan(swarm(i).geolat)) = [];
    
    swarm(i).qdlat = rot90(QDLATS{i});
    swarm(i).qdlat(isnan(swarm(i).qdlat)) = [];
    
    swarm(i).F1 = rot90(F1S{i});
    swarm(i).F1(isnan(swarm(i).F1)) = [];
    
    swarm(i).F2 = rot90(F2S{i});
    swarm(i).F2(isnan(swarm(i).F2)) = [];
end

save('./EEJ_Data/Swarm_1HzData.mat', 'swarm', '-v7.3', '-nocompression')

%% SCALAR DATA

% # Field 1: timestamp (seconds since 1970-01-01 00:00:00 UTC)
% # Field 2: time (decimal year)
% # Field 3: longitude (degrees)
% # Field 4: geocentric latitude (degrees)
% # Field 5: QD latitude (degrees)
% # Field 6: geocentric radius (km)
% # Field 7: spatial weight factor
% # Field 8: F scalar measurement (nT)
% # Field 9: F a priori model (nT)

A = cell(1,9);
for i = 1:9
    fname = sprintf('./EEJ_Data/Magnetic_Data/SwarmA_Scalarcol%d.txt', i);
    fp = fopen(fname);
    a = textscan(fp, '%f', 'HeaderLines', 10, 'delimiter', '\r\n'); A{i} = a{1};
    fclose(fp);
end

B = cell(1,9);
for i = 1:9
    fname = sprintf('./EEJ_Data/Magnetic_Data/SwarmB_Scalarcol%d.txt', i);
    fp = fopen(fname);
    b = textscan(fp, '%f', 'HeaderLines', 10, 'delimiter', '\r\n'); B{i} = b{1};
    fclose(fp);
end

scalar = struct('time', [], 'year', [], 'rad', [], 'lon', [], 'geolat', [], 'qdlat', [],...
    'F', []);

TIMES = {A{1}, B{1}};
YEAR = {A{2}, B{2}};
RADS = {A{6}, B{6}};
LONS = {A{3}, B{3}};
GEOLATS = {A{4}, B{4}};
QDLATS = {A{5}, B{5}};
F = {A{8}, B{8}};

for i = 1:2
    scalar(i).time = rot90(TIMES{i});
    scalar(i).time(isnan(scalar(i).time)) = [];
    
    scalar(i).year = rot90(YEAR{i});
    scalar(i).year(isnan(scalar(i).year)) = [];
    
    scalar(i).rad = rot90(RADS{i});
    scalar(i).rad(isnan(scalar(i).rad)) = [];
    
    scalar(i).lon = rot90(LONS{i});
    scalar(i).lon(isnan(scalar(i).lon)) = [];
    
    scalar(i).geolat = rot90(GEOLATS{i});
    scalar(i).geolat(isnan(scalar(i).geolat)) = [];
    
    scalar(i).qdlat = rot90(QDLATS{i});
    scalar(i).qdlat(isnan(scalar(i).qdlat)) = [];
    
    scalar(i).F = rot90(F{i});
    scalar(i).F(isnan(scalar(i).F)) = [];
end

save('./EEJ_Data/Swarm_Scalar.mat', 'scalar', '-v7.3', '-nocompression')

