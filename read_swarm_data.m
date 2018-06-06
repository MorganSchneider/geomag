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
%% Load A
A = dlmread('./EEJ_Data/SwarmA_1Hz.dat', ' ', 15, 0);

%% Load B
B = load('./EEJ_Data/SwarmB_1Hz.dat');

%% Load C
C = load('./EEJ_Data/SwarmC_1Hz.dat');

%% Organize data into struct array

sizeA = size(A);
ind = 1;
for i = 1:sizeA(1)
    row = A(i,:);
    row(row==0) = [];
    if length(row) == 15
        Atemp(ind,1:15) = row;
        ind = ind + 1;
    end
end

%%
Atime = A(:,1);
Alocal = A(:,2);
Arad = A(:,4);
Alon = A(:,5);
Ageolat = A(:,6);
Aqdlat = A(:,7);
AF1 = A(:,11);
AF2 = A(:,14);


swarm = struct('time', [], 'local', [], 'rad', [], 'lon', [], 'geolat', [],...
    'qdlat', [], 'F1', [], 'F2', []);

TIMES = {Atime, Btime, Ctime};
LOCAL = {Alocal, Blocal, Clocal};
RADS = {Arad, Brad, Crad};
LONS = {Alon, Blon, Clon};
GEOLATS = {Ageolat, Bgeolat, Cgeolat};
QDLATS = {Aqdlat, Bqdlat, Cqdlat};
F1S = {AF1, BF1, CF1};
F2S = {AF2, BF2, CF2};

for i = 1:3
    swarm(i).time = rot90(TIMES{i});
    swarm(i).local = rot90(LOCAL{i});
    swarm(i).rad = rot90(RADS{i});
    swarm(i).lon = rot90(LONS{i});
    swarm(i).geolat = rot90(GEOLATS{i});
    swarm(i).qdlat = rot90(QDLATS{i});
    swarm(i).F1 = rot90(F1S{i});
    swarm(i).F2 = rot90(F2S{i});
end

save('./EEJ_Data/Swarm_1HzData.mat', 'swarm')





