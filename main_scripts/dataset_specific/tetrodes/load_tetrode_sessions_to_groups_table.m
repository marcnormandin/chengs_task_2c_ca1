% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [R] = load_tetrode_sessions_to_groups_table()
    R = load_tetrode_animal_day_table();
    R.groupId = R.dayNum;
    R.groupLabel = R.dayLabel;
end % function
