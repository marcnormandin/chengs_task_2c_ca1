% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [hFig] = plot_stability_bfo180_average_similarity_curves_per_period(DaysToSessions, StableBFO180Similarity, UnstableBFO180Similarity)
    TStable = compute_bfo180_average_similarity_curves_per_period_table(DaysToSessions, StableBFO180Similarity);
    TUnstable = compute_bfo180_average_similarity_curves_per_period_table(DaysToSessions, UnstableBFO180Similarity);

    uniqueDays = unique(DaysToSessions.dayNum);
    numDays = length(uniqueDays);

    hFig = figure('position', get(0, 'screensize'));

    typeNamesToPlot = {'within', 'across'};

        coloursYellowPurple = [
            0.929, 0.694, 0.125; ...
            0.494, 0.184, 0.556];

    %     coloursBlueRed = [...
    %             0, 0.4470, 0.7410; ...
    %             0.8500, 0.3250, 0.0980];

    colours = coloursYellowPurple;

    ax = [];
    for iDay = 1:numDays
        dayNum = uniqueDays(iDay);
        DayDataStable = TStable(TStable.dayNum == dayNum,:);
        DayDataUnstable = TUnstable(TUnstable.dayNum == dayNum,:);

        ax(iDay) = subplot(1,numDays,iDay);
        
        % Stable
        for iType = 1:length(typeNamesToPlot)
           x = DayDataStable.(sprintf('%s_x', typeNamesToPlot{iType})){:};
           y = DayDataStable.(sprintf('%s_y', typeNamesToPlot{iType})){:};

           plot(x,y,'-', 'color', colours(iType,:), 'linewidth', 2)
           hold on
           xlabel('Correlation')
           ylabel('Cumulative')
           axis equal tight
        end
        
        % Unstable
        for iType = 1:length(typeNamesToPlot)
           x = DayDataUnstable.(sprintf('%s_x', typeNamesToPlot{iType})){:};
           y = DayDataUnstable.(sprintf('%s_y', typeNamesToPlot{iType})){:};

           plot(x,y,':', 'color', colours(iType,:), 'linewidth', 2)
           hold on
           xlabel('Correlation')
           ylabel('Cumulative')
           axis equal tight
        end
        
        title(sprintf('%s', DayDataStable.dayLabel{1}))
        
        legendLabels = {};
        for i = 1:length(typeNamesToPlot)
            legendLabels{length(legendLabels)+1} = sprintf('stable %s', typeNamesToPlot{i});
        end
        for i = 1:length(typeNamesToPlot)
            legendLabels{length(legendLabels)+1} = sprintf('unstable %s', typeNamesToPlot{i});
        end
        
        legend(legendLabels, 'location', 'southoutside', 'orientation', 'vertical')
        grid on
        grid minor
    end
    linkaxes(ax, 'xy')
    sgtitle('Cumulative Similarity')
end % function
