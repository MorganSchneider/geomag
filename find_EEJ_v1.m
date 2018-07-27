% function [peakTime, peakLats, peakLons, peakRads] = find_EEJ(F1, t_start, t_end)

% A routine for finding the location of the EEJ using Swarm data.
%
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

% Separate orbits

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

% Separate days

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
        A.day(i).time = [A.day(i).time, A.orbit(ind).time];
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


% Get datetime strings for the starting time of each orbit

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
for i = 1:nDays_A
    A.day(i).start = time_str_A(15*(i-1)+1,:);
end
for i = 1:nOrbits_B
    B.orbit(i).start = time_str_B(i,:);
end
for i = 1:nDays_B
    B.day(i).start = time_str_B(15*(i-1)+1,:);
end
for i = 1:nOrbits_C
    C.orbit(i).start = time_str_C(i,:);
end
for i = 1:nDays_C
    C.day(i).start = time_str_C(15*(i-1)+1,:);
end


% Some plots
% 
% figure(1)
% subplot(2, 1, 1)
% suptitle(sprintf('Swarm A %14s', A.day(1).start))
% plot(A.day(1).qdlat, A.day(1).F1, '.b')
% title('F1')
% subplot(2, 1, 2)
% plot(A.day(1).qdlat, A.day(1).F2, '.b')
% title('F2')
% 
% 
% 
% figure(2)
% hold on
% for i = 1:15
%     plot(A.day(i).geolat, A.day(i).F2, '*r')
%     plot(B.day(i).geolat, B.day(i).F2, '*b')
%     plot(C.day(i).geolat, C.day(i).F2, '*g')
% end
% hold off

%% METHOD 1: Find mean and standard deviation of baseline F2

for i = 1:nDays_A
    A.day(i).qd_meanF2 = mean(A.day(i).F2(abs(A.day(i).qdlat) > 5));
    A.day(i).qd_stdvF2 = std(A.day(i).F2(abs(A.day(i).qdlat) > 5));
end
for i = 1:nDays_B
    B.day(i).qd_meanF2 = mean(B.day(i).F2(abs(B.day(i).qdlat) > 5));
    B.day(i).qd_stdvF2 = std(B.day(i).F2(abs(B.day(i).qdlat) > 5));
end
for i = 1:nDays_C
    C.day(i).qd_meanF2 = mean(C.day(i).F2(abs(C.day(i).qdlat) > 5));
    C.day(i).qd_stdvF2 = std(C.day(i).F2(abs(C.day(i).qdlat) > 5));
end


qdInds_A = cell(nOrbits_A, 1);
for i = 1:nDays_A
    for j = 1:15
        ind = 15 * (i-1) + j;
        qdInds_A{ind} = find(abs(A.orbit(ind).F2) > A.day(i).qd_meanF2 + 3*A.day(i).qd_stdvF2);
    end
end
qdInds_B = cell(nOrbits_B, 1);
for i = 1:nDays_B
    for j = 1:15
        ind = 15 * (i-1) + j;
        qdInds_B{ind} = find(abs(B.orbit(ind).F2) > B.day(i).qd_meanF2 + 3*B.day(i).qd_stdvF2);
    end
end
qdInds_C = cell(nOrbits_C, 1);
for i = 1:nDays_C
    for j = 1:15
        ind = 15 * (i-1) + j;
        qdInds_C{ind} = find(abs(C.orbit(ind).F2) > C.day(i).qd_meanF2 + 3*C.day(i).qd_stdvF2);
    end
end

% for i = 1:15
%     figure
%     subplot(2,1,1)
%     plot(A.orbit(i).geolat, A.orbit(i).F1, '.b', A.orbit(i).geolat(qdInds_A{i}), A.orbit(i).F1(qdInds_A{i}), '*m');
%     subplot(2,1,2)
%     plot(A.orbit(i).geolat, A.orbit(i).F2, '.b', A.orbit(i).geolat(qdInds_A{i}), A.orbit(i).F2(qdInds_A{i}), '*m');
% end

%% METHOD 2: Use gradient of F1 to find maxima and minima

gradInds_A = cell(nOrbits_A, 1);
for i = 1:nOrbits_A
    A.orbit(i).dF1 = gradient(A.orbit(i).F1);
    ii = [];
    for j = 1:length(A.orbit(i).dF1) - 1
        if abs(A.orbit(i).geolat(j)) < 15 && A.orbit(i).dF1(j) * A.orbit(i).dF1(j+1) < 0
            ii = [ii, j, j+1];
        end
        gradInds_A{i} = ii;
    end
end
gradInds_B = cell(nOrbits_B, 1);
for i = 1:nOrbits_B
    B.orbit(i).dF1 = gradient(B.orbit(i).F1);
    ii = [];
    for j = 1:length(B.orbit(i).dF1) - 1
        if abs(B.orbit(i).geolat(j)) < 15 && B.orbit(i).dF1(j) * B.orbit(i).dF1(j+1) < 0
            ii = [ii, j, j+1];
        end
        gradInds_B{i} = ii;
    end
end
gradInds_C = cell(nOrbits_C, 1);
for i = 1:nOrbits_C
    C.orbit(i).dF1 = gradient(C.orbit(i).F1);
    ii = [];
    for j = 1:length(C.orbit(i).dF1) - 1
        if abs(C.orbit(i).geolat(j)) < 15 && C.orbit(i).dF1(j) * C.orbit(i).dF1(j+1) < 0
            ii = [ii, j, j+1];
        end
        gradInds_C{i} = ii;
    end
end

% for i = 1:15
%     figure
%     subplot(2,1,1)
%     plot(A.orbit(i).geolat, A.orbit(i).F1, '.b', A.orbit(i).geolat(gradInds_A{i}), A.orbit(i).F1(gradInds_A{i}), '*m');
%     subpplot(2,1,2)
%     plot(A.orbit(i).geolat, A.orbit(i).F2, '.b', A.orbit(i).geolat(gradInds_A{i}), A.orbit(i).F2(gradInds_A{i}), '*m');
% end

%% Combine the two methods to find the index of each peak
peakInds_A = cell(nOrbits_A, 1);
peakLats_A = zeros(1, nOrbits_A);
peakLons_A = zeros(1, nOrbits_A);
peakRads_A = zeros(1, nOrbits_A);
peakTime_A = zeros(1, nOrbits_A);
for i = 1:nOrbits_A
    [~, locb] = ismember(gradInds_A{i}, qdInds_A{i});
    peakInds_A{i} = qdInds_A{i}(locb(locb ~= 0));
    if ~isempty(peakInds_A{i})
        peakLats_A(i) = nanmean(A.orbit(i).geolat(peakInds_A{i}));
        peakLons_A(i) = nanmean(A.orbit(i).lon(peakInds_A{i}));
        peakRads_A(i) = nanmean(A.orbit(i).rad(peakInds_A{i}));
        peakTime_A(i) = nanmean(A.orbit(i).time(peakInds_A{i}));
    else
        peakLats_A(i) = nan;
        peakLons_A(i) = nan;
        peakRads_A(i) = nan;
        peakTime_A(i) = nan;
    end
end
peakLats_A(isnan(peakLats_A)) = [];
peakLons_A(isnan(peakLons_A)) = [];
peakRads_A(isnan(peakRads_A)) = [];
peakTime_A(isnan(peakTime_A)) = [];

peakInds_B = cell(nOrbits_B, 1);
peakLats_B = zeros(1, nOrbits_B);
peakLons_B = zeros(1, nOrbits_B);
peakRads_B = zeros(1, nOrbits_B);
peakTime_B = zeros(1, nOrbits_B);
for i = 1:nOrbits_B
    [~, locb] = ismember(gradInds_B{i}, qdInds_B{i});
    peakInds_B{i} = qdInds_B{i}(locb(locb ~= 0));
    if ~isempty(peakInds_B{i})
        peakLats_B(i) = nanmean(B.orbit(i).geolat(peakInds_B{i}));
        peakLons_B(i) = nanmean(B.orbit(i).lon(peakInds_B{i}));
        peakRads_B(i) = nanmean(B.orbit(i).rad(peakInds_B{i}));
        peakTime_B(i) = nanmean(B.orbit(i).time(peakInds_B{i}));
    else
        peakLats_B(i) = nan;
        peakLons_B(i) = nan;
        peakRads_B(i) = nan;
        peakTime_B(i) = nan;
    end
end
peakLats_B(isnan(peakLats_B)) = [];
peakLons_B(isnan(peakLons_B)) = [];
peakRads_B(isnan(peakRads_B)) = [];
peakTime_B(isnan(peakTime_B)) = [];

peakInds_C = cell(nOrbits_C, 1);
peakLats_C = zeros(1, nOrbits_C);
peakLons_C = zeros(1, nOrbits_C);
peakRads_C = zeros(1, nOrbits_C);
peakTime_C = zeros(1, nOrbits_C);
for i = 1:nOrbits_C
    [~, locb] = ismember(gradInds_C{i}, qdInds_C{i});
    peakInds_C{i} = qdInds_C{i}(locb(locb ~= 0));
    if ~isempty(peakInds_C{i})
        peakLats_C(i) = nanmean(C.orbit(i).geolat(peakInds_C{i}));
        peakLons_C(i) = nanmean(C.orbit(i).lon(peakInds_C{i}));
        peakRads_C(i) = nanmean(C.orbit(i).rad(peakInds_C{i}));
        peakTime_C(i) = nanmean(C.orbit(i).time(peakInds_C{i}));
    else
        peakLats_C(i) = nan;
        peakLons_C(i) = nan;
        peakRads_C(i) = nan;
        peakTime_C(i) = nan;
    end
end
peakLats_C(isnan(peakLats_C)) = [];
peakLons_C(isnan(peakLons_C)) = [];
peakRads_C(isnan(peakRads_C)) = [];
peakTime_C(isnan(peakTime_C)) = [];

peakLats = {peakLats_A, peakLats_B, peakLats_C};
peakLons = {peakLons_A, peakLons_B, peakLons_C};
peakRads = {peakRads_A, peakRads_B, peakRads_C};
peakTime = {peakTime_A, peakTime_B, peakTime_C};



% return
