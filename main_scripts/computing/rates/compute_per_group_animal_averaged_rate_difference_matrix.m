% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [T] = compute_per_group_animal_averaged_rate_difference_matrix(SessionsToGroups, RateMatrices, perCellNormalizationMethod, perAnimalNormalizationMethod)
    uniquePairs = unique_groups_sessions_to_groups(SessionsToGroups);
    groupIds = [uniquePairs{:,1}];
    groupLabels = {uniquePairs{:,2}};
    numGroups = length(groupIds);
    
    T = [];
    k = 1;
    
    for iGroup = 1:numGroups
        groupId = groupIds(iGroup);
        groupLabel = groupLabels{iGroup};
        
        GroupSessions = SessionsToGroups(SessionsToGroups.groupId == groupId,:);

        MatchTable = get_match_table_by_sessions(GroupSessions, RateMatrices);
        if ~isempty(MatchTable)

            % Average PER MOUSE (Treat each mouse equally)
            groupRateMatricesMice = [];

            numCellsAveraged = 0;
            numMice = size(GroupSessions,1); % This assumes that not more than one mouse session is mapped to a group
            for iMouse = 1:numMice
                animalName = GroupSessions.animalName{iMouse};
                sessionName = GroupSessions.sessionName{iMouse};
                % Get all of the data for ONE MOUSE, on ONE DAY
                RM = RateMatrices(MatchTable.TableBRowIndex(ismember(MatchTable.animalName, animalName) & ismember(MatchTable.sessionName, sessionName)),:);

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

                    numCellsAveraged = numCellsAveraged + 1;
                end

                if size(RM,1) ~= 0

                    % Should be the same for all the data in the main table.
                    rateMatrixTrialIds = RM.rateMatrixTrialIds(1,:);
                    rateMatrixContextIds = RM.rateMatrixContextIds(1,:);

                    groupRateMatrixMouse = nanmean(rateMatrices, 3);

                    % Add the current mouses average map to the collection
                    groupRateMatricesMice = cat(3, groupRateMatricesMice, groupRateMatrixMouse);
                end
            end % iMouse
            
            % At this point we have normalized per mouse for the given
            % day/group.

            % Normalization PER animal
            if strcmpi(perAnimalNormalizationMethod, 'none')
                groupRateMatrix = nanmean(groupRateMatricesMice, 3);
            elseif strcmpi(perAnimalNormalizationMethod, 'totalsum')
                groupRateMatricesMice = groupRateMatricesMice ./ nansum(groupRateMatricesMice, [1,2]);
                groupRateMatrix = nanmean(groupRateMatricesMice, 3);
            end

            % Store
            T(k).groupId = groupId;
            T(k).groupLabel = groupLabel;
            T(k).rateMatrixMean = groupRateMatrix;

            T(k).perAnimalNormalizationMethod = perAnimalNormalizationMethod;
            T(k).perCellNormalizationMethod = perCellNormalizationMethod;
            T(k).numCellsAveraged = numCellsAveraged;

            T(k).numAnimalsAveraged = length(unique(GroupSessions.animalName)); % accounts if more than one animal session was used per group
            T(k).animalNames = GroupSessions.animalName;
            T(k).sessionNames = GroupSessions.sessionName;
            T(k).rateMatrixTrialIds = rateMatrixTrialIds;
            T(k).rateMatrixContextIds = rateMatrixContextIds;

            k = k + 1;
        end % ~isempty(matchTable)
        
    end % iDay
    
    if ~isempty(T)
        T = struct2table(T);
    end
    
end % function
