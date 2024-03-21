% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [BestAligned] = bestaligned_add_registered_cell_names(BestAligned, CellRegTable) 
    if isempty(CellRegTable)
        error('There is no cellreg table to retrieve names from.\n');
    end
    
    % Make a map from session cell names to registered cell names
    mmap = containers.Map(CellRegTable.cellName, CellRegTable.registeredCellName);
    
    
    rnMap = get_map_registered_cell_names_to_registered_cell_id(CellRegTable);

    BestAligned.registeredCellName = cell(size(BestAligned,1),1);
    BestAligned.registeredCellId = nan(size(BestAligned,1),1);
    for iRow = 1:size(BestAligned,1)
        rcn = mmap(BestAligned.cellName{iRow});
        BestAligned.registeredCellName{iRow} = rcn;
        BestAligned.registeredCellId(iRow) = rnMap(rcn);
    end

%     registeredCellNames = cell(size(BestAligned,1),1);
%     for iCell = 1:size(BestAligned,1)
%        cellName = BestAligned.cellName{iCell};
% 
%        ind = find(ismember(CellRegTable.cellName, cellName));
%        if length(ind) == 1
%            registeredCellNames{iCell} = CellRegTable.registeredCellName{ind};
%        else
%            fprintf('Cell name has no registration: %s\n', cellName);
%        end
% 
%     end % iCell
%     BestAligned.registeredCellName = registeredCellNames;
end % function
