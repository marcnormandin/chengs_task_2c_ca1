% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [hFig] = plot_per_group_animal_averaged_rate_difference_matrix_table(T, RATE_MATRIX_CLIM)

    %T = compute_per_day_animal_averaged_rate_difference_matrix(DaysToSessions, RateMatrices, perCellNormalizationMethod, perAnimalNormalizationMethod);

    hFig = figure('position', get(0, 'screensize'));

    % Keep track of the colour bar limits so we can make them the same
    % across the subplots for comparisons
    cMin = [];
    cMax = [];
    

    numGroups = size(T,1);
   
    PP = 1;
    QQ = numGroups;
    
    ax = [];
    for iGroup = 1:numGroups
        ax(iGroup) = subplot(PP,QQ,iGroup);
                        
        plot_cell_rate_difference_matrix_by_context(T.rateMatrixMean{iGroup}, T.rateMatrixTrialIds(iGroup,:), T.rateMatrixContextIds(iGroup,:));
        title(sprintf('%s animal averaged (%s, %s)', T.groupLabel{iGroup}, T.perCellNormalizationMethod{iGroup}, T.perAnimalNormalizationMethod{iGroup}))
        
        cl = caxis;
        cMin = cat(1, cMin, cl(1));
        cMax = cat(1, cMax, cl(2));
    end
    
    % Set the colour bars
    if isempty(RATE_MATRIX_CLIM)
        cMin = min(cMin);
        cMax = max(cMax);
    else
        cMin = RATE_MATRIX_CLIM(1);
        cMax = RATE_MATRIX_CLIM(2);
    end
    
    for iDay = 1:numGroups
        caxis(ax(iDay), [cMin, cMax]);
    end
    
    sgtitle('Tetrodes');
        
end % function
