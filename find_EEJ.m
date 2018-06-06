function [pt, plat, plon, prad, nOrbits, nPeaks] = find_EEJ(sat, s, method)
%
% A routine for finding the location of the EEJ using Swarm data.
%
% INPUTS       sat:     structure array of all the data from one satellite
%              s:       integer number of the satellite whose data is being used
%              method:  'i' for interpolation or 'm' for mean to determine method of estimating equator position
%
% OUTPUTS      pt:      numerical vector of times where magnetic peaks are found (in seconds since Jan 1 1970)
%              plat:    numerical vector of latitudes where magnetic peaks are found (in degrees from equator)
%              plon:    numerical vector of longitudes where magnetic peaks are found (in degrees from prime meridian)
%              prad:    numerical vector of radii where magnetic peaks are found (in meters from the center of the earth)
%              nOrbits: number of orbits made over the total time period
%              nPeaks:  number of peaks found in each orbit
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

% load('./EEJ_Data/Swarm_Data.mat')
% 
% sat = swarm;
% s = 1;
% method = 'i'; %'m'

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

orbit = struct();
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


time_unix = zeros(1, nOrbits);
for i = 1:nOrbits
    time_unix(i) = orbit(i).time(1);
end
time_str = datestr((datenum('1970', 'yyyy') + time_unix ./ 8.64e4), 'yyyymmdd HH:MM:SS');
for i = 1:nOrbits
    orbit(i).start = time_str(i,:);
end


% Method 1
qdInds = cell(nOrbits, 1);
for i = 1:nOrbits
    orbit(i).qd_meanF2 = mean(orbit(i).F2( abs(orbit(i).qdlat) > 5));
    orbit(i).qd_stdvF2 = std(orbit(i).F2( abs(orbit(i).qdlat) > 5));
    qdInds{i} = find(abs(orbit(i).F2) > orbit(i).qd_meanF2 + 3*orbit(i).qd_stdvF2);
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

%% Combine to find peak indices
peakInds = cell(nOrbits, 1);
peakLats = zeros(1, nOrbits);
peakLons = zeros(1, nOrbits);
peakRads = zeros(1, nOrbits);
peakTime = zeros(1, nOrbits);
nPeaks = zeros(1, nOrbits);
for i = 1:nOrbits
    [lia, ~] = ismember(qdInds{i}, gradInds{i});
    peakInds{i} = sort(qdInds{i}(lia ~= 0));
    if ~isempty(peakInds{i})
        nPeaks(i) = floor(length(peakInds{i}) / 2);
        if method == 'i'
            if length(peakInds{i}) == 1
                peakLats(i) = orbit(i).geolat(peakInds{i});
                peakLons(i) = orbit(i).lon(peakInds{i});
                peakRads(i) = orbit(i).rad(peakInds{i});
                peakTime(i) = orbit(i).time(peakInds{i});
            elseif length(peakInds{i}) == 2
                peakLats(i) = interpolate(orbit(i).dF1(peakInds{i}(1)), orbit(i).dF1(peakInds{i}(2)),...
                    orbit(i).geolat(peakInds{i}(1)), orbit(i).geolat(peakInds{i}(2)));
                peakLons(i) = interpolate(orbit(i).dF1(peakInds{i}(1)), orbit(i).dF1(peakInds{i}(2)),...
                    orbit(i).lon(peakInds{i}(1)), orbit(i).lon(peakInds{i}(2)));
                peakRads(i) = interpolate(orbit(i).dF1(peakInds{i}(1)), orbit(i).dF1(peakInds{i}(2)),...
                    orbit(i).rad(peakInds{i}(1)), orbit(i).rad(peakInds{i}(2)));
                peakTime(i) = interpolate(orbit(i).dF1(peakInds{i}(1)), orbit(i).dF1(peakInds{i}(2)),...
                    orbit(i).time(peakInds{i}(1)), orbit(i).time(peakInds{i}(2)));
            elseif length(peakInds{i}) > 2
                normF2 = abs(orbit(i).F2(peakInds{i}) - orbit(i).qd_meanF2);
                jj = sort(normF2);
                ind1 = find(normF2 == jj(end));
                ind2 = find(normF2 == jj(end-1));
                inds = [peakInds{i}(min([ind1, ind2])), peakInds{i}(max([ind1, ind2]))];
                peakLats(i) = interpolate(orbit(i).dF1(inds(1)), orbit(i).dF1(inds(2)),...
                    orbit(i).geolat(inds(1)), orbit(i).geolat(inds(2)));
                peakLons(i) = interpolate(orbit(i).dF1(inds(1)), orbit(i).dF1(inds(2)),...
                    orbit(i).lon(inds(1)), orbit(i).lon(inds(2)));
                peakRads(i) = interpolate(orbit(i).dF1(inds(1)), orbit(i).dF1(inds(2)),...
                    orbit(i).rad(inds(1)), orbit(i).rad(inds(2)));
                peakTime(i) = interpolate(orbit(i).dF1(inds(1)), orbit(i).dF1(inds(2)),...
                    orbit(i).time(inds(1)), orbit(i).time(inds(2)));
            end
        elseif method == 'm'
            peakLats(i) = nanmean(orbit(i).geolat(peakInds{i}));
            peakLons(i) = nanmean(orbit(i).lon(peakInds{i}));
            peakRads(i) = nanmean(orbit(i).rad(peakInds{i}));
            peakTime(i) = nanmean(orbit(i).time(peakInds{i}));
        end
    else
        peakLats(i) = nan;
        peakLons(i) = nan;
        peakRads(i) = nan;
        peakTime(i) = nan;
        nPeaks(i) = 0;
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