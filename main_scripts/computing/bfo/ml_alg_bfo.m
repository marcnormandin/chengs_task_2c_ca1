% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [v_max, vind_max, rotationsDegUsed ] = ml_alg_bfo(rotationsDeg, maps1, maps2, disjointMaps, MIN_DIFFERENCE_FOR_MAX_CORRELATION)
    numRotations = length(rotationsDeg);
    rotations90 = floor(rotationsDeg ./ 90);
    rotationsDegUsed = rotations90 * 90;

    numMaps1 = size(maps1, 3);
    numMaps2 = size(maps2, 3);

    v_max = [];
    vind_max = [];

    
    for iMap1 = 1:numMaps1
        pm1 = maps1(:,:, iMap1);

        if isempty(pm1) || ~any(pm1, 'all')
            continue; % skip
        end
        

        for iMap2 = 1:numMaps2
            pm2 = maps2(:,:, iMap2);
            if isempty(pm2) || ~any(pm2, 'all')
                continue; % skip
            end
            
            % If the list of maps are the same, then skip certain
            % comparisons
            if ~disjointMaps && iMap2 <= iMap1
                continue; % skip it
            end
            
            %if iMap1 ~= iMap2
                r = zeros(1,numRotations);
%                 figure
%                 subplot(1,numRotations+1,1);
%                 imagesc(pm1)
%                 axis equal tight off
                
                for iRot = 1:numRotations
                    % Rotate T2 counter-clockwise
                    pm2Rot = rot90(pm2, rotations90(iRot));

                    a = 1:numel(pm2); % Use all of the map area

                    x1 = pm1(a);
                    x2 = pm2Rot(a);

                    x1 = reshape(x1, numel(x1), 1);
                    x2 = reshape(x2, numel(x2), 1);

                    r(iRot) = corr(x1, x2);
                    
%                     subplot(1,numRotations+1,iRot+1);
%                     imagesc(pm2Rot)
%                     axis equal tight off
%                     title(sprintf('%d, %0.2f', rotationsDeg(iRot), r(iRot)))
                end
                
                % Add the extra criteria that the maximum found must be a
                % certain amount above the second maximum value
                rSorted = sort(r); % max will be in the last element
                valueIsValid = false;
                if rSorted(end) - MIN_DIFFERENCE_FOR_MAX_CORRELATION >= rSorted(end-1)
                    valueIsValid = true;
                end
                
                % Find the actual maximum and the index of it.
                mv = max(r);
                mi = find(r == mv);
                
                
                vn = mv;
                vindn = mi;
%             else
%                 % To not be biased include the map with itself, but we can't
%                 % use the correlation coefficient function because the
%                 % standard deviation is 0 so will give NAN.
%                 vn = 1;
%                 vindn = 1;
%             end

            if length(mi) > 1
                warning('more than one rotation amount is equally likely so the maps are most likely too sparse.');
            end

            if isfinite(vn) && length(mi) == 1 && valueIsValid
                v_max = [v_max, vn];
                vind_max = [vind_max, vindn];
            end
        end % iMap2
    end % iMap1
    
%     hc = histcounts(vind_max, 1:(numRotations+1));
%     mv = max(hc);
%     mi = find(hc == mv);
%     
%     vind_max = mi;
%     v_max = nan(1, length(vind_max));
end % function
