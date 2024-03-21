% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [T, G] = compute_per_group_animal_averaged_rate_difference_histdata(SessionsToGroups, RateMatrices, perCellNormalizationMethod, perAnimalNormalizationMethod)
    % 2022-02-04. This function was based off of
    % compute_per_group_animal_averaged_rate_difference_matrix after a
    % discussion today (at the Friday meeting). For each day, the code will
    % compute and store one average "within" and one average "across" rate
    % difference per mouse. These are then to be used with the stats.
    
    uniquePairs = unique_groups_sessions_to_groups(SessionsToGroups);
    groupIds = [uniquePairs{:,1}];
    groupLabels = {uniquePairs{:,2}};
    numGroups = length(groupIds);
    
    T = [];
    k = 1;
    
    G = [];
    kk = 1;
    
    for iGroup = 1:numGroups
        groupId = groupIds(iGroup);
        groupLabel = groupLabels{iGroup};
        
        GroupSessions = SessionsToGroups(SessionsToGroups.groupId == groupId,:);
        
        withinValues = [];
        acrossValues = [];

        MatchTable = get_match_table_by_sessions(GroupSessions, RateMatrices);
        if ~isempty(MatchTable)

            % Average PER MOUSE (Treat each mouse equally)
            %groupRateMatricesMice = [];

            numMice = size(GroupSessions,1); % This assumes that not more than one mouse session is mapped to a group
            for iMouse = 1:numMice
                animalName = GroupSessions.animalName{iMouse};
                sessionName = GroupSessions.sessionName{iMouse};
                % Get all of the data for ONE MOUSE, on ONE DAY
                RM = RateMatrices(MatchTable.TableBRowIndex(ismember(MatchTable.animalName, animalName) & ismember(MatchTable.sessionName, sessionName)),:);

                % These are gauranteed to be the same for all cells for the
                % same mouse on the same day, so just get one of them.
                rateMatrixContextIds = RM.rateMatrixContextIds(1,:);
                
                % Loop over all the cells maps for the animal and session
                % combination and apply whatever normalization is being
                % used.
                rateMatrices = [];
                for iM = 1:size(RM,1)
                    rm = RM.rateMatrix{iM};

                    if strcmpi(perCellNormalizationMethod, 'minmax')
                        s = (rm - min(rm,[], 'all')) ./ (max(rm, [], 'all') - min(rm, [], 'all'));
                    elseif strcmpi(perCellNormalizationMethod, 'none')
                        s = rm;
                    elseif strcmpi(perCellNormalizationMethod, 'totalsum')
                        s = rm ./ nansum(rm, 'all');
                    else
                        error('Invalid perCellNormalizationMethod');
                    end
                    rateMatrices = cat(3, rateMatrices, s);
                end
                
                numAnimalCellsAveraged = size(RM,1);

                if size(RM,1) ~= 0

                    % Should be the same for all the data in the main table.
                    %rateMatrixTrialIds = RM.rateMatrixTrialIds(1,:);
                    %rateMatrixContextIds = RM.rateMatrixContextIds(1,:);

                    % Average over the cells of the mouse
                    rateMatrixMouse = nanmean(rateMatrices, 3);
                    
                    % Normalization PER animal (if desired)
                    if strcmpi(perAnimalNormalizationMethod, 'totalsum')
                        rateMatrixMouse = rateMatrixMouse ./ nansum(groupRateMatrixMouse, [1,2]);
                    end
                    
                    % Now we need to compute the mean within and the mean
                    % across
                    [withinMean, withinError, acrossMean, acrossError] = helper_compute_average_within_and_across(rateMatrixMouse, rateMatrixContextIds);

                    % STORE
                    T(k).animalName = animalName;
                    T(k).sessionName = sessionName;
                    T(k).groupId = groupId;
                    T(k).groupLabel = groupLabel;
                    T(k).rateMatrixMouse = rateMatrixMouse;

                    T(k).perAnimalNormalizationMethod = perAnimalNormalizationMethod;
                    T(k).perCellNormalizationMethod = perCellNormalizationMethod;
                    T(k).numAnimalCellsAveraged = numAnimalCellsAveraged;
                    
                    T(k).withinMean = withinMean;
                    T(k).withinError = withinError;
                    
                    T(k).acrossMean = acrossMean;
                    T(k).acrossError = acrossError;
                    

                    k = k + 1;
                    
                    withinValues = [withinValues, withinMean];
                    acrossValues = [acrossValues, acrossMean];
                end
            end % iMouse
            
            % Now compute the group stats
            G(kk).groupId = groupId;
            G(kk).groupLabel = groupLabel;
            G(kk).numAnimals = numMice;
            
            G(kk).withinMean = nanmean(withinValues);
            G(kk).withinError = nanstd(withinValues, 1) ./ sqrt(length(withinValues));
            
            G(kk).acrossMean = nanmean(acrossValues);
            G(kk).acrossError = nanstd(acrossValues, 1) ./ sqrt(length(acrossValues));
            kk = kk + 1;


        end % ~isempty(matchTable)
        
    end % iDay
    
    if ~isempty(T)
        T = struct2table(T);
    end
    
    if ~isempty(G)
        G = struct2table(G);
    end
end % function


function [withinMean, withinError, acrossMean, acrossError] = helper_compute_average_within_and_across(rateMatrix, rateMatrixContextIds)
    withinData = [];
    acrossData = [];
    for i = 1:length(rateMatrixContextIds)
        for j = 1+1:length(rateMatrixContextIds)
            x = rateMatrix(i,j);
            if rateMatrixContextIds(i) ~= rateMatrixContextIds(j)
                acrossData = cat(1, acrossData, x);
            else
                withinData = cat(1, withinData, x);
            end
        end
    end

    withinMean = nanmean(withinData, 'all');
    withinError = nanstd(withinData, 1, 'all') ./ sqrt(length(withinData));

    acrossMean = nanmean(acrossData, 'all');
    acrossError = nanstd(acrossData, 1, 'all') ./ sqrt(length(acrossData));
end  % function
