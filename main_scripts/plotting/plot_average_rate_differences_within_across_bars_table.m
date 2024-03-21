% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [hFig] = plot_average_rate_differences_within_across_bars_table(RateMatricesAveragedHistogram)
    %R = compute_average_rate_differences_within_across_for_table(T);

    hFig = figure();
    
    if isempty(RateMatricesAveragedHistogram)
        return;
    end
    
    coloursBlueRed = [...
            0, 0.4470, 0.7410; ...
            0.8500, 0.3250, 0.0980];
        
    Y = [];
    E = [];
    Y(:,1) = RateMatricesAveragedHistogram.withinMean;
    Y(:,2) = RateMatricesAveragedHistogram.acrossMean;
    E(:,1) = RateMatricesAveragedHistogram.withinError;
    E(:,2) = RateMatricesAveragedHistogram.acrossError;

    Y = Y';
    E = E';

    
    [b] = barplot_with_errors(Y,E, coloursBlueRed);
    legend({'within', 'across'}, 'location', 'southoutside', 'orientation', 'horizontal')
    grid on
    grid minor
    xticks(1:size(Y,2))
    xticklabels(RateMatricesAveragedHistogram.groupLabel)

end % function

function [b] = barplot_with_errors(Y,E, colours)
    Y = Y';
    E = E';
    
    x = 1:size(Y,1);
    
    numGroups = size(Y,1);
    numBars = size(Y,2);
    
    b = bar(x, Y, 'grouped', 'FaceColor','flat');
    for k = 1:size(Y,2)
        b(k).CData = colours(k,:);
    end
    
    hold on
    groupWidth = min(0.8, numBars / (numBars + 1.5));
    for i = 1:numBars
        xx = (1:numGroups) - groupWidth/2 + (2*i-1)*groupWidth / (2*numBars);
        errorbar(xx, Y(:,i), E(:,i), 'k.', 'LineWidth', 2);
    end
end % function
