% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [Classifications] = classify_stability_registered_cells_on_specified_day(BestAligned, CellRegTable, SessionsToGroups, CLASSIFICATION_GROUP_LABEL, BESTALIGNED_STABILITY_THRESHOLD_CRITERIA)
    % This code returns a list of registered cell names and there
    % classifications ON A GIVEN DAY. This is to be used for the stability
    % plots that use registered cells.
    
    % notation
    X = BestAligned;
    
    % Add the group labels
    X = helper_add_group_labels_to_table(X, SessionsToGroups, true);
    
    % Eliminate any data that isn't on the classification day
    X(~ismember(X.groupLabel, CLASSIFICATION_GROUP_LABEL),:) = [];
    
    % Defaults
    X.isStable = false(size(X,1),1);
    X.registeredCellName = cell(size(X,1),1);
    
    % Loopity loop loop
    for iRow = 1:size(X,1)
       % Session cell name
       cellName = X.cellName{iRow};

       % Classify it
       if X.bestCorrelation(iRow) >= BESTALIGNED_STABILITY_THRESHOLD_CRITERIA
        X.isStable(iRow) = true;
       end

       % Find the registered cell name
       %registeredCellName
       ind = find(ismember(CellRegTable.cellName, cellName));
       if ~isempty(ind)
           X.registeredCellName{iRow} = CellRegTable.registeredCellName{ind};
       else
           fprintf('Unable to find: %s\n', cellName);
       end
    end
    Classifications = X(:, ismember(X.Properties.VariableNames, {'registeredCellName', 'isStable'}));
    % There are a lot of empties
    Classifications(cellfun(@isempty,Classifications.registeredCellName),:) = [];
end % function
