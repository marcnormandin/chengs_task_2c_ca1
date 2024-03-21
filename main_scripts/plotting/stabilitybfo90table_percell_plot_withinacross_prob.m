% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [hFig] = stabilitybfo90table_percell_plot_withinacross_prob(groupLabels, StableMeanTable, StableErrorTable, UnstableMeanTable, UnstableErrorTable, showExcluded)    
    numGroups = length(groupLabels); % should be a cell array of strings
    
    hFig = figure('position', get(0, 'screensize'));

    % Each group gets a column with two subplots
    ax = [];
    for iGroup = 1:numGroups
        groupLabel = groupLabels{iGroup};
        
        % Stable
        ax(length(ax)+1) = subplot(2,numGroups,iGroup);
        if ~isempty(StableMeanTable) && ~isempty(StableErrorTable)
            my_helper(StableMeanTable, StableErrorTable, groupLabel, showExcluded);
        end
        title(sprintf('%s Stable', groupLabel))
        
        % Unstable
        ax(length(ax)+1) = subplot(2,numGroups,iGroup + numGroups);
        if ~isempty(UnstableMeanTable) && ~isempty(UnstableErrorTable)
            my_helper(UnstableMeanTable, UnstableErrorTable, groupLabel, showExcluded);
        end
        title(sprintf('%s Unstable', groupLabel))
    end
    if ~isempty(ax)
        linkaxes(ax, 'xy');
    end
    
end % function

function [Y,E] = my_helper(MeanTable, ErrorTable, groupLabel, showExcluded)
        MeanTable(~ismember(MeanTable.groupLabel, groupLabel),:) = [];
        ErrorTable(~ismember(ErrorTable.groupLabel, groupLabel), :) = [];
        
        Y = [];
        E = [];
        
        if isempty(MeanTable) || isempty(ErrorTable)
            return;
        end
        
        %colours4 = colour_matrix_blue_red_yellow_purple();
        colours = colour_matrix_yellow_purple();
        
        if showExcluded
            Y = zeros(2,5);
            E = zeros(2,5);
        else
            Y = zeros(2,4);
            E = zeros(2,4);
        end
        
        %x1 = sprintf('%s_within_prob', stabilityType);
        %x1e = strrep(x1, 'prob', 'excluded_prob');

        %x2 = sprintf('%s_across_prob', stabilityType);
        %x2e = strrep(x2, 'prob', 'excluded_prob');
        
        x1 = 'prob_within';
        x1e = 'prob_within';
        x2 = 'prob_across';
        x2e = 'prob_across';
        
        % Within
        Y(1,1:4) = MeanTable.(x1);
        E(1,1:4) = ErrorTable.(x1);
        
        % Across
        Y(2,1:4) = MeanTable.(x2);
        E(2,1:4) = ErrorTable.(x2);
        
        if showExcluded
            Y(1,5) = MeanTable.(x1e);
            E(1,5) = ErrorTable.(x1e);
            Y(2,5) = MeanTable.(x2e);
            E(2,5) = ErrorTable.(x2e);
        end

        b = barplot_with_errors(Y./100, E./100, colours);
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
        legend({'within', 'across'});
end % function
