% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [R] = compute_withheld_trial_two_labels(Maps, mapLabels, NORMALIZE_MAPS)
    % Maps needs to be a (height, width, num cells) matrix
    % Map labels needs to be the same size as num cells and only contain
    % one of two labels for each map
    % If normalize maps is true, then the maps are normalized to have a
    % total sum of 1, and so they are weighted equally for the average
    % maps.
    
    % Filter out maps that all zeros or all nans
    indNan = find(squeeze(all(isnan(Maps), [1,2])));
    indZero = find(squeeze(all(Maps==0, [1,2])));
    indBad = union(indNan, indZero);
    Maps(:,:,indBad) = [];
    mapLabels(indBad) = [];
    
    if isempty(Maps) || isempty(mapLabels)
        error('No valid data\n');
    end
    
    nMaps = size(Maps,3);
    nLabels = numel(mapLabels);
    if nMaps ~= nLabels
        error('The number of labels must match the number of maps.\n');
    end
    

    if NORMALIZE_MAPS
        for iMap = 1:size(Maps,3)
            mu = Maps(:,:,iMap);
            mn = mu ./ sum(mu, 'all', 'omitnan');
            Maps(:,:,iMap) = mu;
        end
    end

    uniqueLabels = unique(mapLabels);
    if length(uniqueLabels) ~= 2
        error('We need exectly two types of maps only, but we have %d', length(uniqueLabels));
    end
    umap = containers.Map(uniqueLabels, [1,2]);
    c = cellfun(@(x)umap(x), mapLabels);

    R = [];
    k = 1;

    for iMap = 1:nMaps
        % The map we are testing
        mapT = Maps(:,:,iMap);
        cT = c(iMap); 

        mapAverage1 = mean(Maps(:,:, c'==1 & iMap ~= 1:nMaps), 3, 'omitnan');
        mapAverage2 = mean(Maps(:,:, c'==2 & iMap ~= 1:nMaps), 3, 'omitnan');

        R1 = corrcoef(mapAverage1(:), mapT(:));
        rho1 = R1(1,2);

        R2 = corrcoef(mapAverage2(:), mapT(:));
        rho2 = R2(1,2);

        bestMatch = nan;
        %isInvalid = false;
        if isfinite(rho1) && isfinite(rho2)
            if rho1 >= rho2
                bestMatch = 1;
            else
                bestMatch = 2;
            end
            isValid = true;
        else
            isValid = false;
        end

        R(k).isValid = isValid;
        
        R(k).mapAverage1 = mapAverage1;
        R(k).mapLabel1 = uniqueLabels{1};
        R(k).rho1 = rho1;

        R(k).mapAverage2 = mapAverage2;
        R(k).mapLabel2 = uniqueLabels{2};
        R(k).rho2 = rho2;

        R(k).mapWithheld = mapT;
        R(k).mapWithheldType = cT;
        R(k).mapWithheldLabel = uniqueLabels{cT};

        R(k).bestMatch = bestMatch;

        if ~isnan(bestMatch)
            R(k).bestMatchLabel = uniqueLabels{bestMatch};
        else
            R(k).bestMatchLabel = nan;
        end
        k = k + 1;
    end
    R = struct2table(R);
end % function
