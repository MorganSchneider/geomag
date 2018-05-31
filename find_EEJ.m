%% This will eventually be a script for finding the location of the EEJ

% Parameters:
%
% Time (seconds since Jan 1, 1970 - epoch time)
% Radius (meters from the center of the earth)
% Longitude (degrees)
% Geocentric latitude (degrees)
% Quasi-dipole latitude (degrees)
% Field intensity (nT)
% Field intensity corrected for Sq (nT)

load('./EEJ_Data/Swarm_Data.mat')

%% Find the starting indices of each individual orbit

k = 1;
A_inds = zeros(1,10000);
for i = 1:length(swarm(1).time) - 1
    m = swarm(1).time(i+1);
    n = swarm(1).time(i);
    if m - n ~= 5
        A_inds(k) = i + 1;
        k = k + 1;
    end
end
nOrbits_A = k;

k = 1;
B_inds = zeros(1,20000);
for i = 1:length(swarm(2).time) - 1
   m = swarm(2).time(i+1);
   n = swarm(2).time(i);
   if m - n ~= 5
       B_inds(k) = i + 1;
       k = k + 1;
   end
end
nOrbits_B = k;

k = 1;
C_inds = zeros(1,20000);
for i = 1:length(swarm(3).time) - 1
    m = swarm(3).time(i+1);
    n = swarm(3).time(i);
    if m - n ~= 5
        C_inds(k) = i + 1;
        k = k + 1;
    end
end
nOrbits_C = k;

% Vectors of the starting indices of each orbit
A_inds = [A_inds(A_inds ~= 0) length(swarm(1).time)+1];
B_inds = [B_inds(B_inds ~= 0) length(swarm(2).time)+1];
C_inds = [C_inds(C_inds ~= 0) length(swarm(3).time)+1];

%% Separate orbits

% A.orbit: 2505 orbits, 167 days
% B.orbit: 19773 orbits, 1318 days
% C.orbit: 19913 orbits, 1327 days
%
% Swarm A time range: 20131126 00:38:00 - 20140522 11:36:35 (why only 6 months for this one?)
% Swarm B time range: 20131128 00:00:00 - 20170715 22:48:04
% Swarm C time range: 20131204 00:59:41 - 20170715 23:47:18
%
% One orbit is ~33 minutes
% ~1 hour and 3 min between orbits
% One orbit "cycle" is ~1 hour and 35 minutes
% ~15 orbits per day

nDays_A = floor(nOrbits_A / 15);
nDays_B = floor(nOrbits_B / 15);
nDays_C = floor(nOrbits_C / 15);

j = 1;
for i = 1:nOrbits_A
    orbit_inds = j: A_inds(i) - 1;
    
    A.orbit(i).time = swarm(1).time(orbit_inds);
    A.orbit(i).rad = swarm(1).rad(orbit_inds);
    A.orbit(i).lon = swarm(1).lon(orbit_inds);
    A.orbit(i).geolat = swarm(1).geolat(orbit_inds);
    A.orbit(i).qdlat = swarm(1).qdlat(orbit_inds);
    A.orbit(i).F1 = swarm(1).F1(orbit_inds);
    A.orbit(i).F2 = swarm(1).F2(orbit_inds);
    
    j = A_inds(i);
end

j = 1;
for i = 1:nOrbits_B
    orbit_inds = j: B_inds(i) - 1;
    
    B.orbit(i).time = swarm(2).time(orbit_inds);
    B.orbit(i).rad = swarm(2).rad(orbit_inds);
    B.orbit(i).lon = swarm(2).lon(orbit_inds);
    B.orbit(i).geolat = swarm(2).geolat(orbit_inds);
    B.orbit(i).qdlat = swarm(2).qdlat(orbit_inds);
    B.orbit(i).F1 = swarm(2).F1(orbit_inds);
    B.orbit(i).F2 = swarm(2).F2(orbit_inds);
    
    j = B_inds(i);
end

j = 1;
for i = 1:nOrbits_C
    orbit_inds = j: C_inds(i) - 1;
    
    C.orbit(i).time = swarm(3).time(orbit_inds);
    C.orbit(i).rad = swarm(3).rad(orbit_inds);
    C.orbit(i).lon = swarm(3).lon(orbit_inds);
    C.orbit(i).geolat = swarm(3).geolat(orbit_inds);
    C.orbit(i).qdlat = swarm(3).qdlat(orbit_inds);
    C.orbit(i).F1 = swarm(3).F1(orbit_inds);
    C.orbit(i).F2 = swarm(3).F2(orbit_inds);
    
    j = C_inds(i);
end

%% Separate days

for i = 1:nDays_A
    A.day(i).time = [];
    A.day(i).rad = [];
    A.day(i).lon = [];
    A.day(i).geolat = [];
    A.day(i).qdlat = [];
    A.day(i).F1 = [];
    A.day(i).F2 = [];
    for j = 1:15
        ind = (i-1) * 15 + j;
        A.day(i).time = [A.day(i).tifigureme, A.orbit(ind).time];
        A.day(i).rad = [A.day(i).rad, A.orbit(ind).rad];
        A.day(i).lon = [A.day(i).lon, A.orbit(ind).lon];
        A.day(i).geolat = [A.day(i).geolat, A.orbit(ind).geolat];
        A.day(i).qdlat = [A.day(i).qdlat, A.orbit(ind).qdlat];
        A.day(i).F1 = [A.day(i).F1, A.orbit(ind).F1];
        A.day(i).F2 = [A.day(i).F2, A.orbit(ind).F2];
    end
end

for i = 1:nDays_B
    B.day(i).time = [];
    B.day(i).rad = [];
    B.day(i).lon = [];
    B.day(i).geolat = [];
    B.day(i).qdlat = [];
    B.day(i).F1 = [];
    B.day(i).F2 = [];
    for j = 1:15
        ind = (i-1) * 15 + j;
        B.day(i).time = [B.day(i).time, B.orbit(ind).time];
        B.day(i).rad = [B.day(i).rad, B.orbit(ind).rad];
        B.day(i).lon = [B.day(i).lon, B.orbit(ind).lon];
        B.day(i).geolat = [B.day(i).geolat, B.orbit(ind).geolat];
        B.day(i).qdlat = [B.day(i).qdlat, B.orbit(ind).qdlat];
        B.day(i).F1 = [B.day(i).F1, B.orbit(ind).F1];
        B.day(i).F2 = [B.day(i).F2, B.orbit(ind).F2];
    end
end

for i = 1:nDays_C
    C.day(i).time = [];
    C.day(i).rad = [];
    C.day(i).lon = [];
    C.day(i).geolat = [];
    C.day(i).qdlat = [];
    C.day(i).F1 = [];
    C.day(i).F2 = [];
    for j = 1:15
        ind = (i-1) * 15 + j;
        C.day(i).time = [C.day(i).time, C.orbit(ind).time];
        C.day(i).rad = [C.day(i).rad, C.orbit(ind).rad];
        C.day(i).lon = [C.day(i).lon, C.orbit(ind).lon];
        C.day(i).geolat = [C.day(i).geolat, C.orbit(ind).geolat];
        C.day(i).qdlat = [C.day(i).qdlat, C.orbit(ind).qdlat];
        C.day(i).F1 = [C.day(i).F1, C.orbit(ind).F1];
        C.day(i).F2 = [C.day(i).F2, C.orbit(ind).F2];
    end
end


%% Some time bullshit

time_unix_A = zeros(1, nOrbits_A);
for i = 1:nOrbits_A
    time_unix_A(i) = A.orbit(i).time(1);
end
time_unix_B = zeros(1, nOrbits_B);
for i = 1:nOrbits_B
    time_unix_B(i) = B.orbit(i).time(1);
end
time_unix_C = zeros(1, nOrbits_C);
for i = 1:nOrbits_C
    time_unix_C(i) = C.orbit(i).time(1);
end

time_str_A = datestr((datenum('1970', 'yyyy') + time_unix_A ./ 8.64e4), 'yyyymmdd HH:MM:SS');
time_str_B = datestr((datenum('1970', 'yyyy') + time_unix_B ./ 8.64e4), 'yyyymmdd HH:MM:SS');
time_str_C = datestr((datenum('1970', 'yyyy') + time_unix_C ./ 8.64e4), 'yyyymmdd HH:MM:SS');

for i = 1:nOrbits_A
    A.orbit(i).start = time_str_A(i,:);
end
for i = 1:nOrbits_B
    B.orbit(i).start = time_str_B(i,:);
end
for i = 1:nOrbits_C
    C.orbit(i).start = time_str_C(i,:);
end


%% Some plots or something

figure(1)
subplot(2, 1, 1)
hold on
for i = 1:15
    title('Swarm A 11-26-2013')
    plot(A.orbit(i).geolat, A.orbit(i).F1)
    if i <= 5
        title('F1')
    end
end
hold off
title('F1, 11-26-2013')
subplot(2, 1, 2)
hold on
for i = 1:15
    plot(A.orbit(i).geolat, A.orbit(i).F2)
end
hold off
title('F2, 11-26-2013')



figure(2)
hold on
for i = 1:15
    plot(A.day(i).geolat, A.day(i).F2, '*r')
    plot(B.day(i).geolat, B.day(i).F2, '*b')
    plot(C.day(i).geolat, C.day(i).F2, '*g')
end
hold off

%% Some more junk


for i = 1:nDays_A
    A.day(i).meanF2 = mean([A.day(i).F2(A.day(i).geolat >= 15), A.day(i).F2(A.day(i).geolat <= -15)]);
    A.day(i).stdvF2 = std([A.day(i).F2(A.day(i).geolat >= 15), A.day(i).F2(A.day(i).geolat <= -15)]);
end

for i = 1:nDays_B
    B.day(i).meanF2 = mean([B.day(i).F2(B.day(i).geolat >= 15), B.day(i).F2(B.day(i).geolat <= -15)]);
    B.day(i).stdvF2 = std([B.day(i).F2(B.day(i).geolat >= 15), B.day(i).F2(B.day(i).geolat <= -15)]);
end

for i = 1:nDays_C
    C.day(i).meanF2 = mean([C.day(i).F2(C.day(i).geolat >= 15), C.day(i).F2(C.day(i).geolat <= -15)]);
    C.day(i).stdvF2 = std([C.day(i).F2(C.day(i).geolat >= 15), C.day(i).F2(C.day(i).geolat <= -15)]);
end






