%% This will eventually be a script for finding the location of the EEJ

load('./EEJ_Data/Swarm_Data.mat')

%% Find the starting indices of each separate orbit

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

%% Find the starting indices of each individual orbit
A_inds = [A_inds(A_inds ~= 0) length(swarm(1).time)+1];
B_inds = [B_inds(B_inds ~= 0) length(swarm(2).time)+1];
C_inds = [C_inds(C_inds ~= 0) length(swarm(3).time)+1];

%% Separate orbits

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








