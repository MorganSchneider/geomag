function [pt, plat, plon, prad, pqd, nOrbits, nPeaks] = find_EEJ(sat, s)
%
% A routine for finding the location of the EEJ using Swarm data.
%
% INPUTS       sat:     structure array of all the data from one satellite
%              s:       integer number of the satellite whose data is being used
%              method:  'interp' or 'mean' to determine method of estimating equator position
%
% OUTPUTS      pt:      numerical vector of times where magnetic peaks are found (in seconds since Jan 1 1970)
%              plat:    numerical vector of latitudes where magnetic peaks are found (in degrees from equator)
%              plon:    numerical vector of longitudes where magnetic peaks are found (in degrees from prime meridian)
%              prad:    numerical vector of radii where magnetic peaks are found (in meters from the center of the earth)
%              ploc:    numerical vector of local times where magnetic peaks are found (in hours)
%              nOrbits: number of orbits made over the total time period
%              nPeaks:  number of peaks found in each orbit
%
% Parameters:
%
% Time (seconds since Jan 1, 1970 - epoch time)
% Local time (hours)
% Radius (meters from the center of the earth)
% Longitude (degrees)
% Geocentric latitude (degrees)
% Quasi-dipole latitude (degrees)
% Field intensity (nT)
% Field intensity corrected for Sq (nT)
%
% ***FOR DEBUGGING***
%
% load('./EEJ_Data/Swarm_1HzData.mat')
% 
% sat = swarm;
% s = 1;

%% Organize data by orbit

k = 1;
inds = zeros(1,1e6);
for i = 1:length(sat(s).time) - 1
    m = sat(s).time(i+1);
    n = sat(s).time(i);
    if m - n > 1
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
    orbit(i).local = sat(s).local(orbit_inds);
    orbit(i).rad = sat(s).rad(orbit_inds);
    orbit(i).lon = sat(s).lon(orbit_inds);
    orbit(i).geolat = sat(s).geolat(orbit_inds);
    orbit(i).qdlat = sat(s).qdlat(orbit_inds);
    orbit(i).F1 = sat(s).F1(orbit_inds);
    orbit(i).F2 = sat(s).F2(orbit_inds);
    
    j = inds(i);
end

% filter by local time?
dayOrbits = [];
for i = 1:nOrbits
    if nanmean(orbit(i).local) > 10 || nanmean(orbit(i).local) < 14
        dayOrbits = [dayOrbits, i];
    end
end

%% Find possible peaks

% Method 1
qdInds = cell(nOrbits, 1);
for i = dayOrbits
    qdi_upper = find(orbit(i).qdlat > 10 & orbit(i).qdlat <= 55);
    qdi_lower = find(orbit(i).qdlat < -10 & orbit(i).qdlat >= -55);
    qdi = [qdi_lower qdi_upper];
    orbit(i).qd_meanF2 = mean(orbit(i).F2(qdi));
    orbit(i).qd_stdvF2 = std(orbit(i).F2(qdi));
    qdInds{i} = find(abs(orbit(i).F2 - orbit(i).qd_meanF2) > 4.5*orbit(i).qd_stdvF2 & abs(orbit(i).qdlat) < 3);
end

% Method 2
gradInds = cell(nOrbits, 1);
for i = dayOrbits
    orbit(i).dF1 = gradient(orbit(i).F1);
    orbit(i).dF2 = gradient(orbit(i).F2); %%%% just for debugging
    ii = [];
    for j = 1:length(orbit(i).dF1) - 1
        if abs(orbit(i).geolat(j)) < 12 && orbit(i).dF1(j) * orbit(i).dF1(j+1) <= 0
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
peakLocal = zeros(1, nOrbits);
nPeaks = zeros(1, nOrbits);
peakF1 = zeros(1, nOrbits); %%%% just for debugging
peakF2 = zeros(1, nOrbits); %%%% just for debugging
peakQd = zeros(1, nOrbits); %%%% just for debugging
for i = dayOrbits
    [lia, ~] = ismember(qdInds{i}, gradInds{i});
    peakInds{i} = sort(qdInds{i}(lia ~= 0));
    if ~isempty(peakInds{i})
        nPeaks(i) = ceil(length(peakInds{i}) / 2);
        if length(peakInds{i}) == 1
            peakLats(i) = orbit(i).geolat(peakInds{i});
            peakLons(i) = orbit(i).lon(peakInds{i});
            peakRads(i) = orbit(i).rad(peakInds{i});
            peakTime(i) = orbit(i).time(peakInds{i});
            peakLocal(i) = orbit(i).local(peakInds{i});
            peakF1(i) = orbit(i).F1(peakInds{i}); %%%% just for debugging
            peakF2(i) = orbit(i).F2(peakInds{i}); %%%% just for debugging
            peakQd(i) = orbit(i).qdlat(peakInds{i}); %%%% just for debugging
        elseif length(peakInds{i}) == 2
            if peakInds{i}(2) - peakInds{i}(1) == 1
                peakLats(i) = interpolate(orbit(i).dF1(peakInds{i}(1)), orbit(i).dF1(peakInds{i}(2)),...
                    orbit(i).geolat(peakInds{i}(1)), orbit(i).geolat(peakInds{i}(2)));
                peakLons(i) = interpolate(orbit(i).dF1(peakInds{i}(1)), orbit(i).dF1(peakInds{i}(2)),...
                    orbit(i).lon(peakInds{i}(1)), orbit(i).lon(peakInds{i}(2)));
                peakRads(i) = interpolate(orbit(i).dF1(peakInds{i}(1)), orbit(i).dF1(peakInds{i}(2)),...
                    orbit(i).rad(peakInds{i}(1)), orbit(i).rad(peakInds{i}(2)));
                peakTime(i) = interpolate(orbit(i).dF1(peakInds{i}(1)), orbit(i).dF1(peakInds{i}(2)),...
                    orbit(i).time(peakInds{i}(1)), orbit(i).time(peakInds{i}(2)));
                peakLocal(i) = interpolate(orbit(i).dF1(peakInds{i}(1)), orbit(i).dF1(peakInds{i}(2)),...
                    orbit(i).local(peakInds{i}(1)), orbit(i).local(peakInds{i}(2)));
                peakF1(i) = mean([orbit(i).F1(peakInds{i}(1)) + 0.5*orbit(i).dF1(peakInds{i}(1)),...
                    orbit(i).F1(peakInds{i}(2)) - 0.5*orbit(i).dF1(peakInds{i}(2))]); %%%% just for debugging
                peakF2(i) = mean([orbit(i).F2(peakInds{i}(1)) + 0.5*orbit(i).dF2(peakInds{i}(1)),...
                    orbit(i).F2(peakInds{i}(2)) - 0.5*orbit(i).dF2(peakInds{i}(2))]); %%%% just for debugging
                peakQd(i) = interpolate(orbit(i).dF1(peakInds{i}(1)), orbit(i).dF1(peakInds{i}(2)),...
                    orbit(i).qdlat(peakInds{i}(1)), orbit(i).qdlat(peakInds{i}(2))); %%%% just for debugging
            else
                peakLats(i) = nan;
                peakLons(i) = nan;
                peakRads(i) = nan;
                peakTime(i) = nan;
                peakLocal(i) = nan;
                peakF1(i) = nan; %%%% just for debugging
                peakF2(i) = nan; %%%% just for debugging
                peakQd(i) = nan; %%%% just for debugging
            end
        elseif length(peakInds{i}) > 2
            normF2 = abs(orbit(i).F2(peakInds{i}) - orbit(i).qd_meanF2);
            maxdiff1 = max(normF2);
            ind1 = find(normF2 == maxdiff1);
            if length(ind1) > 1
                ind1 = ind1(1);
                normF2(ind1(2:end)) = [];
            end
            if ind1 + 1 > length(normF2)
                x1 = false;
            else
                x1 = true;
            end
            if ind1 - 1 < 1
                x2 = false;
            else
                x2 = true;
            end
            
            if xor(x1,x2)
                if x1 == 1
                    maxdiff2 = normF2(ind1+1);
                elseif x2 == 1
                    maxdiff2 = normF2(ind1-1);
                end
            elseif x1 == 1 && x2 == 1
                maxdiff2 = max(normF2(ind1+1), normF2(ind1-1));
            end
            ind2 = find(normF2 == maxdiff2);
            inds = [peakInds{i}(min([ind1, ind2])), peakInds{i}(max([ind1, ind2]))];
            peakLats(i) = interpolate(orbit(i).dF1(inds(1)), orbit(i).dF1(inds(2)),...
                orbit(i).geolat(inds(1)), orbit(i).geolat(inds(2)));
            peakLons(i) = interpolate(orbit(i).dF1(inds(1)), orbit(i).dF1(inds(2)),...
                orbit(i).lon(inds(1)), orbit(i).lon(inds(2)));
            peakRads(i) = interpolate(orbit(i).dF1(inds(1)), orbit(i).dF1(inds(2)),...
                orbit(i).rad(inds(1)), orbit(i).rad(inds(2)));
            peakTime(i) = interpolate(orbit(i).dF1(inds(1)), orbit(i).dF1(inds(2)),...
                orbit(i).time(inds(1)), orbit(i).time(inds(2)));
            peakLocal(i) = interpolate(orbit(i).dF1(inds(1)), orbit(i).dF1(inds(2)),...
                orbit(i).local(inds(1)), orbit(i).local(inds(2)));
            peakF1(i) = mean([orbit(i).F1(inds(1)) + 0.5*orbit(i).dF1(inds(1)),...
                orbit(i).F1(inds(2)) - 0.5*orbit(i).dF1(inds(2))]); %%%% just for debugging
            peakF2(i) = mean([orbit(i).F2(inds(1)) + 0.5*orbit(i).dF2(inds(1)),...
                orbit(i).F2(inds(2)) - 0.5*orbit(i).dF2(inds(2))]); %%%% just for debugging
            peakQd(i) = interpolate(orbit(i).dF1(inds(1)), orbit(i).dF1(inds(2)),...
                orbit(i).qdlat(inds(1)), orbit(i).qdlat(inds(2))); %%%% just for debugging
        end
    else
        peakLats(i) = nan;
        peakLons(i) = nan;
        peakRads(i) = nan;
        peakTime(i) = nan;
        peakLocal(i) = nan;
        peakF1(i) = nan; %%%% just for debugging
        peakF2(i) = nan; %%%% just for debugging
        peakQd(i) = nan; %%%% just for debugging
        nPeaks(i) = 0;
    end
end

pt = peakTime(~isnan(peakTime));
plat = peakLats(~isnan(peakTime));
pcol = 90 - plat;
plon = peakLons(~isnan(peakTime));
prad = peakRads(~isnan(peakTime));
ploc = peakLocal(~isnan(peakTime));
pf1 = peakF1(~isnan(peakTime)); %%%% just for debugging
pf2 = peakF2(~isnan(peakTime)); %%%% just for debugging
pqd = peakQd(~isnan(peakTime)); %%%% just for debugging
nPeaks = nPeaks(~isnan(peakTime));

%% Filter and correct

nPeaks(nPeaks == 0) = [];
pt = pt(nPeaks == 1);
pcol = pcol(nPeaks == 1);
plon = plon(nPeaks == 1);
prad = prad(nPeaks == 1);
ploc = ploc(nPeaks == 1);
pf1 = pf1(nPeaks == 1);
pf2 = pf2(nPeaks == 1);
pqd = pqd(nPeaks == 1);

pcol = pcol + 0.1;
plat = 90 - pcol;

[timestamp, k_p, r_c] = read_indices(sat);

kp = zeros(1, length(pt));
rc = zeros(1, length(pt));
for i = 1:length(pt)
    z = find(timestamp <= pt(i));
    ind = z(end);
    kp(i) = k_p(ind);
    rc(i) = r_c(ind);
end
x1 = find(kp <= 2);
x2 = find(abs(rc) <= 5);
lia = ismember(x1, x2);
quiet = x1(lia == 1);

pt = pt(quiet);
plat = plat(quiet);
pcol = pcol(quiet);
plon = plon(quiet);
prad = prad(quiet);
ploc = ploc(quiet);
pf1 = pf1(quiet);
pf2 = pf2(quiet);
pqd = pqd(quiet);


return