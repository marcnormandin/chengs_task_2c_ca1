% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [SessionsToGroups] = load_sessions_to_groups_table(analysisSettings)
    if ~analysisSettings.IS_CALCIUM_DATA
        SessionsToGroups = load_tetrode_sessions_to_groups_table();
    else
        SessionsToGroups = load_calcium_sessions_to_groups_table();
    end
end % function
