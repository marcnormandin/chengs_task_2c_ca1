% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function T = calcium_compute_fraction_occupied(folder, useSmoothed)
    tmp = load(fullfile(folder ,'calcium_occupancy_maps.mat'));

    T = table2struct(tmp.OccupancyMaps);
    for iRow = 1:length(T)
        if useSmoothed
            occupancyMap = T(iRow).occupancyMapSmoothed;
        else

            occupancyMap = T(iRow).occupancyMap;
        end

        R = compute_whatever_we_need_per_map(occupancyMap, 2);
        
        fieldNames = fieldnames(R);
        for iField = 1:length(fieldNames)
            fieldName = fieldNames{iField};
            T(iRow).(fieldName) = R.(fieldName);
        end
    end
    T = struct2table(T);
end
