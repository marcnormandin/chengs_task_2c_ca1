% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [mmap2] = get_map_registered_cell_names_to_registered_cell_id(CellRegTable) 
    if isempty(CellRegTable)
        error('There is no cellreg table to retrieve names from.\n');
    end
    
    uniqueRegisteredCellNames = unique(CellRegTable.registeredCellName);
    mmap2 = containers.Map(uniqueRegisteredCellNames, 1:length(uniqueRegisteredCellNames));

end % function
