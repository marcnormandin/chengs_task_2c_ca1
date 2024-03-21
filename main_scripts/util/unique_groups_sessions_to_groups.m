% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [uniquePairs] = unique_groups_sessions_to_groups(SessionsToGroups)
    % This mess is just to get a Nx2 cell array of unique combinations of a
    % groupId and a groupLabel.

    uniqueGroupLabels = unique(SessionsToGroups.groupLabel);
    groupLabelToIndMap = containers.Map(uniqueGroupLabels, 1:length(uniqueGroupLabels));
    indToGroupLabelMap = containers.Map(1:length(uniqueGroupLabels), uniqueGroupLabels);

    groupLabelInds = [];
    for i = 1:size(SessionsToGroups,1)
        groupLabelInds(i) = groupLabelToIndMap(SessionsToGroups.groupLabel{i});
    end
    groupLabelInds = reshape(groupLabelInds, size(SessionsToGroups,1), 1);

    uniquePairs0 = unique([SessionsToGroups.groupId, groupLabelInds], 'rows');

    uniquePairs = cell(size(uniquePairs0,1), 2);
    for i = 1:size(uniquePairs0,1)
       uniquePairs{i,1} = uniquePairs0(i,1);
       uniquePairs{i,2} = indToGroupLabelMap(uniquePairs0(i,2));
    end
end % function
