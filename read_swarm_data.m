% Testing how to best read swarm data.

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

fileA = 'home/mschneider/EEJ_Data/SwarmA.dat'; %1048577 - 15 lines (including empty lines)
fileB = 'home/mschneider/EEJ_Data/SwarmB.dat'; %8323902 - 15 lines
fileC = 'home/mschneider/EEJ_Data/SwarmC.dat'; %8294271 - 15 lines
fileAx = 'home/mschneider/EEJ_Data/swarmAbreaks.txt';
fileBx = 'home/mschneider/EEJ_Data/swarmBbreaks.txt';
fileCx = 'home/mschneider/EEJ_Data/swarmCbreaks.txt';

fA = fopen(fileAx);
A_brx_cell = textscan(fA, '%d');
A_brx = rot90(cell2mat(A_brx_cell));
fB = fopen(fileBx);
B_brx_cell = textscan(fB, '%d');
B_brx = rot90(cell2mat(B_brx_cell));
fC = fopen(fileCx);
C_brx_cell = textscan(fC, '%d');
C_brx = rot90(cell2mat(C_brx_cell));

A_brx = [15, A_brx, 1048578];
B_brx = [15, B_brx, 8323903];
C_brx = [15, C_brx, 8294272];

%%
nOrbits_A = length(A_brx)/2;
nOrbits_B = length(B_brx)/2;
nOrbits_C = length(C_brx)/2;
A_cell = cell(15, nOrbits_A);
B_cell = cell(15, nOrbits_B);
C_cell = cell(15, nOrbits_C);

for i = 1:15
    for j = 1:nOrbits_A
        A_cell{i,j} = dlmread(fileA, '', [A_brx(2*j-1) i-1 A_brx(2*j)-2 i-1]);
    end
    for j = 1:nOrbits_B
        B_cell{i,j} = dlmread(fileB, '', [B_brx(2*j-1) i-1 B_brx(2*j)-2 i-1]);
    end
    for j = 1:nOrbits_C
        C_cell{i,j} = dlmread(fileC, '', [C_brx(2*j-1) i-1 C_brx(2*j)-2 i-1]);
    end
end


%%
for j = 1:nOrbits_A
    swarm(1).time.orbit(j) = A_cell{1,j};
    swarm(1).radius.orbit(j) = A_cell{4,j};
    swarm(1).long.orbit(j) = A_cell{5,j};
    swarm(1).geolat.orbit(j) = A_cell{6,j};
    swarm(1).qdlat.orbit(j) = A_cell{7,j};
    swarm(1).F1.orbit(j) = A_cell{11,j};
    swarm(1).F2.orbit(j) = A_cell{14,j};
end
for j = 1:nOrbits_B
    
end
for j = 1:nOrbits_C
    
end




