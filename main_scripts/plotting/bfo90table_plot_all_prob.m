% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [hFig] = bfo90table_plot_all_prob(MeanTable, ErrorTable, showExcluded)    
    colours4 = colour_matrix_blue_red_yellow_purple();

    Y = MeanTable.all_prob;
    E = ErrorTable.all_prob;

    if showExcluded
        Y(:,5) = MeanTable.all_excluded_prob;
        E(:,5) = ErrorTable.all_excluded_prob;
    end

    hFig = figure('position', get(0, 'screensize'));
    
    b = barplot_with_errors(Y./100, E./100, colours4);

    if showExcluded
        xticks(1:5);
        xticklabels({['0' char(176)], ['90' char(176)], ['180' char(176)], ['270' char(176)], 'Excluded'});
    else
        xticks(1:4);
        xticklabels({['0' char(176)], ['90' char(176)], ['180' char(176)], ['270' char(176)]});
    end
    
    grid on
    ax = gca;
    ax.FontSize = 15;
        
    legend(MeanTable.groupLabel);
end % function
