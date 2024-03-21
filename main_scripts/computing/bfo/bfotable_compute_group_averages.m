% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [MeanTable, ErrorTable] = bfotable_compute_group_averages(BFOTable, SessionsToGroups)
    % This should work with any BFO90 or BFO180 table
    
    uniquePairs = unique_groups_sessions_to_groups(SessionsToGroups);
    groupIds = [uniquePairs{:,1}];
    groupLabels = {uniquePairs{:,2}};
    
    numGroups = length(groupIds);
    MeanTable = [];
    ErrorTable = [];
    for iGroupId = 1:numGroups
        groupId = groupIds(iGroupId);
        groupLabel = groupLabels{iGroupId};
        
        groupTable = SessionsToGroups(SessionsToGroups.groupId == groupId, :);
        matchTable = get_match_table_by_sessions(groupTable, BFOTable);
        
        if ~isempty(matchTable)
            groupBFO90Table = BFOTable(matchTable.TableBRowIndex, :);
            fieldsToAverage = groupBFO90Table.Properties.VariableNames(contains(groupBFO90Table.Properties.VariableNames, 'prob'));

            [T,E] = get_prob_stats_for_table(groupBFO90Table, fieldsToAverage);

            T.groupId = repmat(groupId, size(T,1), 1);
            T.groupLabel = cell(size(T,1),1);
            for i = 1:size(T,1)
                T.groupLabel{i} = groupLabel;
            end

            E.groupId = repmat(groupId, size(E,1), 1);
            E.groupLabel = cell(size(E,1),1);
            for i = 1:size(E,1)
                E.groupLabel{i} = groupLabel;
            end

            MeanTable = [MeanTable; T];
            ErrorTable = [ErrorTable; E];
        end
    end
end % function
