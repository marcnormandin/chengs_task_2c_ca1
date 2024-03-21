% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [hFig] = stabilitybfo90table_percell_plot_withinacross_prob_v1(StableMeanTable, StableErrorTable, UnstableMeanTable, UnstableErrorTable, showExcluded)    
    hasStable = true;
    hasUnstable = true;
    
    numGroups1 = size(StableMeanTable,1);
    numGroups2 = size(UnstableMeanTable,1);
    if isempty(StableMeanTable) || isempty(UnstableMeanTable)
        if isempty(StableMeanTable) && ~isempty(UnstableMeanTable)
            numGroups = numGroups2;
            hasStable = false;
        elseif ~isempty(StableMeanTable) && isempty(UnstableMeanTable)
            hasUnstable = false;
            numGroups = numGroups1;
        elseif isempty(StableMeanTable) && isempty(UnstableMeanTable)
            warning('Both stable and unstable has no data to plot for stability bfo90 per cell')
            hasStable = false;
            hasUnstable = false;
        end % no else?
    elseif numGroups1 ~= numGroups2
        warning('Number of groups (days) for stable (%d) and unstable (%d) do not match!', numGroups1, numGroups2);
        numGroups = max([numGroups1, numGroups2]);
    else
        numGroups = numGroups1; % or numGroups2
    end
    
    hFig = figure('position', get(0, 'screensize'));

    % Each group gets a column with two subplots
    ax = [];
    for iGroup = 1:numGroups
    
        ax(length(ax)+1) = subplot(2,numGroups,iGroup);
        if hasStable
            my_helper(StableMeanTable, StableErrorTable, iGroup, showExcluded);
            title(sprintf('%s Stable', StableMeanTable.groupLabel{iGroup}))
        end
        
        ax(length(ax)+1) = subplot(2,numGroups,iGroup + numGroups);
        if hasUnstable
            my_helper(UnstableMeanTable, UnstableErrorTable, iGroup, showExcluded);
            title(sprintf('%s Unstable', UnstableMeanTable.groupLabel{iGroup}))
        end
    end
    if ~isempty(ax)
        linkaxes(ax, 'xy');
    end
    
end % function

function [Y,E] = my_helper(MeanTable, ErrorTable, iGroup, showExcluded)
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
        Y(1,1:4) = MeanTable.(x1)(iGroup,:);
        E(1,1:4) = ErrorTable.(x1)(iGroup,:);
        
        % Across
        Y(2,1:4) = MeanTable.(x2)(iGroup,:);
        E(2,1:4) = ErrorTable.(x2)(iGroup,:);
        
        if showExcluded
            Y(1,5) = MeanTable.(x1e)(iGroup);
            E(1,5) = ErrorTable.(x1e)(iGroup);
            Y(2,5) = MeanTable.(x2e)(iGroup);
            E(2,5) = ErrorTable.(x2e)(iGroup);
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
