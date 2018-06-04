function [pt, plat, plon, prad] = find_EEJ(sat, s)

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

% Modify to take in time range input
% 
% load('./EEJ_Data/Swarm_Data.mat')
% 


%% Test


k = 1;
inds = zeros(1,10000);
for i = 1:length(sat(s).time) - 1
    m = sat(s).time(i+1);
    n = sat(s).time(i);
    if m - n ~= 5
        inds(k) = i + 1;
        k = k + 1;
    end
end
nOrbits = k;
inds = [inds(inds ~= 0) length(sat(s).time)+1];

nDays = floor(nOrbits / 15);
j = 1;
for i = 1:nOrbits
    orbit_inds = j: inds(i) - 1;
    
    orbit(i).time = sat(s).time(orbit_inds);
    orbit(i).rad = sat(s).rad(orbit_inds);
    orbit(i).lon = sat(s).lon(orbit_inds);
    orbit(i).geolat = sat(s).geolat(orbit_inds);
    orbit(i).qdlat = sat(s).qdlat(orbit_inds);
    orbit(i).F1 = sat(s).F1(orbit_inds);
    orbit(i).F2 = sat(s).F2(orbit_inds);
    
    j = inds(i);
end

for i = 1:nDays
    day(i).time = [];
    day(i).rad = [];
    day(i).lon = [];
    day(i).geolat = [];
    day(i).qdlat = [];
    day(i).F1 = [];
    day(i).F2 = [];
    for j = 1:15
        ind = (i-1) * 15 + j;
        day(i).time = [day(i).time, orbit(ind).time];
        day(i).rad = [day(i).rad, orbit(ind).rad];
        day(i).lon = [day(i).lon, orbit(ind).lon];
        day(i).geolat = [day(i).geolat, orbit(ind).geolat];
        day(i).qdlat = [day(i).qdlat, orbit(ind).qdlat];
        day(i).F1 = [day(i).F1, orbit(ind).F1];
        day(i).F2 = [day(i).F2, orbit(ind).F2];
    end
end

time_unix = zeros(1, nOrbits);
for i = 1:nOrbits
    time_unix(i) = orbit(i).time(1);
end
time_str = datestr((datenum('1970', 'yyyy') + time_unix ./ 8.64e4), 'yyyymmdd HH:MM:SS');
for i = 1:nOrbits
    orbit(i).start = time_str(i,:);
end
for i = 1:nDays
    day(i).start = time_str(15*(i-1)+1,:);
end

% Method 1
qdInds = cell(nOrbits, 1);
for i = 1:nDays
    day(i).qd_meanF2 = mean(day(i).F2( abs(day(i).qdlat) > 5));
    day(i).qd_stdvF2 = std(day(i).F2( abs(day(i).qdlat) > 5));
    for j = 1:15
        ind = 15 * (i-1) + j;
        qdInds{ind} = find(abs(orbit(ind).F2) > day(i).qd_meanF2 + 3*day(i).qd_stdvF2);
    end
end

% Method 2
gradInds = cell(nOrbits, 1);
for i = 1:nOrbits
    orbit(i).dF1 = gradient(orbit(i).F1);
    ii = [];
    for j = 1:length(orbit(i).dF1) - 1
        if abs(orbit(i).geolat(j)) < 15 && orbit(i).dF1(j) * orbit(i).dF1(j+1) < 0
            ii = [ii, j, j+1];
        end
        gradInds{i} = ii;
    end
end

% Combine to find peak indices
peakInds = cell(nOrbits, 1);
peakLats = zeros(1, nOrbits);
peakLons = zeros(1, nOrbits);
peakRads = zeros(1, nOrbits);
peakTime = zeros(1, nOrbits);
for i = 1:nOrbits
    [~, locb] = ismember(gradInds{i}, qdInds{i});
    peakInds{i} = qdInds{i}(locb(locb ~= 0));
    if ~isempty(peakInds{i})
        peakLats(i) = nanmean(orbit(i).geolat(peakInds{i}));
        peakLons(i) = nanmean(orbit(i).lon(peakInds{i}));
        peakRads(i) = nanmean(orbit(i).rad(peakInds{i}));
        peakTime(i) = nanmean(orbit(i).time(peakInds{i}));
    else
        peakLats(i) = nan;
        peakLons(i) = nan;
        peakRads(i) = nan;
        peakTime(i) = nan;
    end
end
peakLats(isnan(peakLats)) = [];
peakLons(isnan(peakLons)) = [];
peakRads(isnan(peakRads)) = [];
peakTime(isnan(peakTime)) = [];

pt = peakTime;
plat = peakLats;
plon = peakLons;
prad = peakRads;

return



