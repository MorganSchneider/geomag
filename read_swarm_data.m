%% Testing how to best read swarm data.

% # Field 1: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)
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

%% Swarm A

% set variable names for files
fileA = 'home/mschneider/EEJ_Data/SwarmA.dat'; %1048577 - 15 lines (including empty lines)
fileAx = 'home/mschneider/EEJ_Data/swarmAbreaks.txt';

%read line break text file created in bash
fA = fopen(fileAx);
A_brx_cell = textscan(fA, '%d');
%convert cell array to regular numeric array and reshape
A_brx = rot90(cell2mat(A_brx_cell));
%add first and last values
A_brx = [15, A_brx, 1048578];

nOrbits_A = length(A_brx)/2;
%allocate cell array
A_cell = cell(15);

%read Swarm data from ASCII file into cell array
for i = 1:15
    A_cell{i} = dlmread(fileA, '', [15, i-1, 1048576, i]);
end

%% Swarm B

fileB = 'home/mschneider/EEJ_Data/SwarmB.dat'; %8323902 - 15 lines
fileBx = 'home/mschneider/EEJ_Data/swarmBbreaks.txt';

fB = fopen(fileBx);
B_brx_cell = textscan(fB, '%d');
B_brx = rot90(cell2mat(B_brx_cell));
B_brx = [15, B_brx, 8323903];

nOrbits_B = length(B_brx)/2;
B_cell = cell(15);

for i = 1:15
    B_cell{i} = dlmread(fileB, '', [B_brx(1) i-1 B_brx(end) i]);
end

%% Swarm C

fileC = 'home/mschneider/EEJ_Data/SwarmC.dat'; %8294271 - 15 lines
fileCx = 'home/mschneider/EEJ_Data/swarmCbreaks.txt';

fC = fopen(fileCx);
C_brx_cell = textscan(fC, '%d');
C_brx = rot90(cell2mat(C_brx_cell));
C_brx = [15, C_brx, 8294272];

nOrbits_C = length(C_brx)/2;
C_cell = cell(15);

for i = 1:15
    C_cell{i} = dlmread(fileC, '', [C_brx(1) i-1 C_brx(end) i]);
end



%% Save to structure arrays

swarm(1).time = A_cell{1};
swarm(1).rad = A_cell{4};
swarm(1).lon = A_cell{5};
swarm(1).geolat = A_cell{6};
swarm(1).qdlat = A_cell{7};
swarm(1).F1 = A_cell{11};
swarm(1).F2 = A_cell{14};


swarm(2).time = B_cell{1};
swarm(2).rad = B_cell{4};
swarm(2).lon = B_cell{5};
swarm(2).geolat = B_cell{6};
swarm(2).qdlat = B_cell{7};
swarm(2).F1 = B_cell{11};
swarm(2).F2 = B_cell{14};

swarm(3).time = C_cell{1};
swarm(3).rad = C_cell{4};
swarm(3).lon = C_cell{5};
swarm(3).geolat = C_cell{6};
swarm(3).qdlat = C_cell{7};
swarm(3).F1 = C_cell{11};
swarm(3).F2 = C_cell{14};





