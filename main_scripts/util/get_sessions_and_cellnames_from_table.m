% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [sessions] = get_sessions_and_cellnames_from_table(T)
    % Returns a unique table of animal name / session name pairs
    animalIndex = find(ismember(T.Properties.VariableNames, 'animalName'));
    sessionIndex = find(ismember(T.Properties.VariableNames, 'sessionName'));
    dayNumIndex = find(ismember(T.Properties.VariableNames, 'dayNum'));
    cellNameIndex = find(ismember(T.Properties.VariableNames, 'cellName'));
    cellIdIndex = find(ismember(T.Properties.VariableNames, 'cellId'));
    animalSessionCellNameIndex = find(ismember(T.Properties.VariableNames, 'animalSessionCellName'));
    
    sessions = unique(T(:,[animalIndex, sessionIndex, dayNumIndex, cellNameIndex, cellIdIndex, animalSessionCellNameIndex]), 'rows');
end
