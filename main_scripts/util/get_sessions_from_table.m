% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [sessions] = get_sessions_from_table(T)
    % Returns a unique table of animal name / session name pairs
    animalIndex = find(ismember(T.Properties.VariableNames, 'animalName'));
    sessionIndex = find(ismember(T.Properties.VariableNames, 'sessionName'));

    sessions = unique(T(:,[animalIndex,sessionIndex]), 'rows');
end
