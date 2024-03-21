% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [hFig] = plot_per_day_cell_averaged_rate_difference_matrix_table(T)

    %T = compute_per_day_cell_averaged_rate_difference_matrix(DaysToSessions, RateMatrices, perCellNormalizationMethod);

    hFig = figure('position', get(0, 'screensize'));

    % Keep track of the colour bar limits so we can make them the same
    % across the subplots for comparisons
    cMin = [];
    cMax = [];
    
    uniqueDays = unique(T.dayNum);
    numDays = length(uniqueDays);
    
    if numDays == 1
        PP = 1;
        QQ = 1;
    elseif numDays == 2
        PP = 1;
        QQ = 2;
    elseif numDays == 4
        PP = 2;
        QQ = 2;
    else
        PP = 1;
        QQ = numDays;
    end
    
    ax = [];
    for iDay = 1:numDays
        ax(iDay) = subplot(PP,QQ,iDay);
        dayNum = uniqueDays(iDay);
        
        DayData = T(T.dayNum == dayNum, :);
        
        rateMatrixTrialIds = DayData.rateMatrixTrialIds(1,:);
        rateMatrixContextIds = DayData.rateMatrixContextIds(1,:);

        dayRateMatrix = DayData.rateMatrixMean{1};

        plot_cell_rate_difference_matrix_by_context(dayRateMatrix, rateMatrixTrialIds, rateMatrixContextIds);
        title(sprintf('%s averaged (%s)', DayData.dayLabel{1}, DayData.perCellNormalizationMethod{1}))
        
        cl = caxis;
        cMin = cat(1, cMin, cl(1));
        cMax = cat(1, cMax, cl(2));
    end
    
    % Set the colour bars
    cMin = min(cMin);
    cMax = max(cMax);
    for iDay = 1:numDays
        caxis(ax(iDay), [cMin, cMax]);
    end
    
    sgtitle('Averaged PER CELL');
end % function
