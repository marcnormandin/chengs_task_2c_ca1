% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [T] = add_registered_cell_names_to_table(T, CellRegTable)
% Add the registered cell names to each of the cells
    registeredCellNames = cell(size(T,1),1);
    for iRow = 1:size(T,1)
        cellName = T.cellName{iRow};
        ind = find(ismember(CellRegTable.cellName, cellName));
        if length(ind) ~= 1
            fprintf('No registered cell found for %s\n', cellName);
        else
            registeredCellNames{iRow} = CellRegTable.registeredCellName{ind};
        end
    end % iRow
    T.registeredCellName = registeredCellNames;
end % function
