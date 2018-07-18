% EEJ ALGORITHM

sat = struct();
for s = 1:3
    k = 1;
    inds = zeros(1,1e6);
    for i = 1:length(swarm(s).time) - 1
        m = swarm(s).time(i+1);
        n = swarm(s).time(i);
        if m - n > 1
            inds(k) = i + 1;
            k = k + 1;
        end
    end
    nOrbits = k;
    inds = [inds(inds ~= 0) length(swarm(s).time)+1];
    
    orbit = struct();
    j = 1;
    for i = 1:nOrbits
        orbit_inds = j: inds(i) - 1;
        
        orbit(i).time = swarm(s).time(orbit_inds);
        orbit(i).local = swarm(s).local(orbit_inds);
        orbit(i).rad = swarm(s).rad(orbit_inds);
        orbit(i).lon = swarm(s).lon(orbit_inds);
        orbit(i).geolat = swarm(s).geolat(orbit_inds);
        orbit(i).qdlat = swarm(s).qdlat(orbit_inds);
        orbit(i).F1 = swarm(s).F1(orbit_inds);
        orbit(i).F2 = swarm(s).F2(orbit_inds);
        
        j = inds(i);
    end
    
    
    % Method 1
    qdInds = cell(nOrbits, 1);
    for i = 1:nOrbits
        qdi_upper = find(orbit(i).qdlat > 10 & orbit(i).qdlat <= 55);
        qdi_lower = find(orbit(i).qdlat < -10 & orbit(i).qdlat >= -55);
        qdi = [qdi_lower qdi_upper];
        orbit(i).qd_meanF2 = mean(orbit(i).F2(qdi));
        orbit(i).qd_stdvF2 = std(orbit(i).F2(qdi));
        qdInds{i} = find(abs(orbit(i).F2 - orbit(i).qd_meanF2) > 4.5*orbit(i).qd_stdvF2 & abs(orbit(i).qdlat) < 10);
        orbit(i).qdInds = qdInds{i};
    end
    
    % Method 2
    gradInds = cell(nOrbits, 1);
    for i = 1:nOrbits
        orbit(i).dF1 = gradient(orbit(i).F1);
        orbit(i).dF2 = gradient(orbit(i).F2); %%%% just for debugging
        ii = [];
        for j = 1:length(orbit(i).dF1) - 1
            if abs(orbit(i).geolat(j)) < 12 && orbit(i).dF1(j) * orbit(i).dF1(j+1) <= 0
                ii = [ii, j, j+1];
            end
            gradInds{i} = ii;
            orbit(i).gradInds = gradInds{i};
        end
    end
    
    
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
    for i = 1:nOrbits
        [lia, ~] = ismember(qdInds{i}, gradInds{i});
        peakInds{i} = sort(qdInds{i}(lia ~= 0));
        orbit(i).peakInds = peakInds{i};
        if ~isempty(peakInds{i})
            nPeaks(i) = ceil(length(peakInds{i}) / 2);
            orbit(i).nPeaks = nPeaks(i);
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
            orbit(i).nPeaks = nPeaks(i);
        end
    end
    
    sat(s).orbit = orbit;
    nOrbit(s) = nOrbits;
    sat(s).pt = peakTime(~isnan(peakTime));
    sat(s).plat = peakLats(~isnan(peakTime));
    sat(s).pcol = 90 - sat(s).plat;
    sat(s).plon = peakLons(~isnan(peakTime));
    sat(s).prad = peakRads(~isnan(peakTime));
    sat(s).ploc = peakLocal(~isnan(peakTime));
    sat(s).pf1 = peakF1(~isnan(peakTime)); %%%% just for debugging
    sat(s).pf2 = peakF2(~isnan(peakTime)); %%%% just for debugging
    sat(s).pqd = peakQd(~isnan(peakTime)); %%%% just for debugging
    
    sat(s).orbitsNoPeaks = find(nPeaks == 0);
    sat(s).orbitsWithPeaks = find(nPeaks ~= 0);
    sat(s).singlePeakOrbits = find(nPeaks == 1);
    sat(s).multPeakOrbits = find(nPeaks > 1);
end

%% Filter (EEJ ALGORITHM)
for s = 1:3
    nPeaks(nPeaks == 0) = [];
    sat(s).pt = sat(s).pt(nPeaks == 1);
    sat(s).pcol = sat(s).pcol(nPeaks == 1);
    sat(s).plon = sat(s).plon(nPeaks == 1);
    sat(s).prad = sat(s).prad(nPeaks == 1);
    sat(s).ploc = sat(s).ploc(nPeaks == 1);
    sat(s).pf1 = sat(s).pf1(nPeaks == 1);
    sat(s).pf2 = sat(s).pf2(nPeaks == 1);
    sat(s).pqd = sat(s).pqd(nPeaks == 1);
    
    sat(s).pcol = sat(s).pcol + 0.1;
    sat(s).plat = 90 - sat(s).pcol;
end
%% Correct (EEJ ALGORITHM)
for s = 1:3
    [timestamp, k_p, r_c] = read_indices(swarm);
    
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
    
    sat(s).pt = sat(s).pt(quiet);
    sat(s).plat = sat(s).plat(quiet);
    sat(s).pcol = sat(s).pcol(quiet);
    sat(s).plon = sat(s).plon(quiet);
    sat(s).prad = sat(s).prad(quiet);
    sat(s).ploc = sat(s).ploc(quiet);
    sat(s).pf1 = sat(s).pf1(quiet);
    sat(s).pf2 = sat(s).pf2(quiet);
    sat(s).pqd = sat(s).pqd(quiet);
end

%%
for s = 1:3
    sat(s).mjd = datenum(datetime(sat(s).pt,'ConvertFrom','posixtime')) - datenum(2000,1,1,0,0,0);
    options = optimset('FunValCheck','on');
    for i = 1:length(sat(s).pt)
        if sat(s).mjd(i) >= pp.breaks(1) && sat(s).mjd(i) <= pp.breaks(end)
            fun = @(theta) findzero(sat(s).prad(i), theta, sat(s).plon(i), pp, sat(s).mjd(i));
            sat(s).chaosLats(i) = fzero(fun, sat(s).pcol(i), options);
        else
            sat(s).chaosLats(i) = nan;
        end
    end
end


%%
resA = colatsA - chaosA;

%% Testing find_EEJ algorithm accuracy

% indices = find(nPeaksA > 1);
% spec = find(resA > 1);

indices = find(nPeaksC > 1);
spec = find(resC > 1);

n = indices(spec(1));

subplot(1,2,1)
plot(orbit(n).geolat, orbit(n).F1, '-b')
hold on
plot(orbit(n).geolat(peakInds{n}), orbit(n).F1(peakInds{n}), '*r')
plot(orbit(n).geolat(gradInds{n}), orbit(n).F1(gradInds{n}), 'og')
plot(peakLats(n), peakF1(n), '*k');
hold off

subplot(1,2,2)
plot(orbit(n).geolat, orbit(n).F2, '-b')
hold on
plot(orbit(n).geolat(peakInds{n}), orbit(n).F2(peakInds{n}), '*r')
plot(orbit(n).geolat(qdInds{n}), orbit(n).F2(qdInds{n}), 'oy')
plot(peakLats(n), peakF2(n), '*k');
hold off

shg


%% Look at orbits with largest residuals to see what is wrong

figure
% max(resA) = 7.941, min(resA) = -7.193
subplot(3,3,1)
histogram(resA)
title('A')
% max(resB) = 7.879, min(resB) = -7.941
subplot(3,3,2)
histogram(resB)
title('B')
% max(resC) = 6.892, min(resC) = -5.970
subplot(3,3,3)
histogram(resC)
title('C')
% max(resA_corr) = 8.041, min(resA_corr) = -7.093
subplot(3,3,4)
histogram(resA_corr)
title('A corrected')
% max(resB_corr) = 7.979, min(resB_corr) = -7.841
subplot(3,3,5)
histogram(resB_corr)
title('B corrected')
% max(resC_corr) = 6.992, min(resC_corr) = -5.870
subplot(3,3,6)
histogram(resC_corr)
title('C corrected')
% max(resA_quiet) = 6.202, min(resA_quiet) = -5.463
subplot(3,3,7)
histogram(resA_quiet)
title('A quiet')
% max(resB_quiet) = 7.609, min(resB_quiet) = -4.262
subplot(3,3,8)
histogram(resB_quiet)
title('B quiet')
% max(resC_quiet) = 4.907, min(resC_quiet) = -4.688
subplot(3,3,9)
histogram(resC_quiet)
title('C quiet')

