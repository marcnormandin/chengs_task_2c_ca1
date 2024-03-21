% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [hFig, ax] = plot_bfo180_average_similarity_curves_per_period(DaysToSessions, NormalBFO180Similarity, ShuffledTable)
    T = compute_bfo180_average_similarity_curves_per_period_table(DaysToSessions, NormalBFO180Similarity);
    uniqueDays = unique(T.dayNum);
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
        DayData = T(T.dayNum == dayNum,:);

        ax(iDay) = subplot(1,numDays,iDay);
        for iType = 1:length(typeNamesToPlot)
           x = DayData.(sprintf('%s_x', typeNamesToPlot{iType})){:};
           y = DayData.(sprintf('%s_y', typeNamesToPlot{iType})){:};

           plot(x,y,'-', 'color', colours(iType,:), 'linewidth', 2)
           hold on
           xlabel('Correlation')
           ylabel('Cumulative')
           axis equal tight
        end
        
        typeNamesPlotted = typeNamesToPlot;

        % Show the shuffled distribution
        if ~isempty(ShuffledTable)
            iRow = find(ismember(ShuffledTable.dayLabel,  sprintf('Day %d', dayNum)));
            if ~isempty(iRow)
                shuffle_x = ShuffledTable.all_averageCorrelations_x{iRow};
                shuffle_y = ShuffledTable.all_averageCorrelations_y{iRow};
                plot(shuffle_x, shuffle_y, 'k-', 'linewidth', 2);
                typeNamesPlotted{end+1} = 'shuffle';
            else
                warning('Shuffled distribution set to be plotted by associated day data not found.\n');
            end
        end
                
        
        title(sprintf('%s', DayData.dayLabel{1}))
        legend(typeNamesPlotted, 'location', 'southoutside', 'orientation', 'horizontal')
        grid on
        grid minor
    end
    linkaxes(ax, 'xy')
    sgtitle('Cumulative Similarity')
end % function
