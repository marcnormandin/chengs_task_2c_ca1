% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [hFig] = plot_per_day_animal_averaged_rate_difference_matrix(DaysToSessions, RateMatrices, perCellNormalizationMethod, perAnimalNormalizationMethod)

    T = compute_per_day_animal_averaged_rate_difference_matrix(DaysToSessions, RateMatrices, perCellNormalizationMethod, perAnimalNormalizationMethod);

    hFig = plot_per_day_animal_averaged_rate_difference_matrix_table(T);
        
end % function
