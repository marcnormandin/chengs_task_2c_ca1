% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024

function [hFig] = make_figure(resultsFn)
    data = readtable(resultsFn); % the Excel file with the accuracies
    
    groupLabels = {'Day 1', 'Day 2', 'Day 3'};
    nGroups = length(groupLabels);

    M = zeros(nGroups,1);
    E = zeros(nGroups,1);
    for iGroup = 1:nGroups
        groupLabel = groupLabels{iGroup};
        S = data(ismember(data.groupLabel, groupLabel),:);
        x = S.accuracy * 100;
        meanX = mean(x);
        errX = std(x) ./ sqrt(length(x));
        M(iGroup) = meanX;
        E(iGroup) = errX;
    end

    hFig = figure();
    bar(1:nGroups, M)
    hold on
    errorbar(1:nGroups, M, E, 'linestyle', 'none', 'linewidth', 2, 'color', 'k');
    xticks(1:nGroups);
    xticklabels(groupLabels)
    hold on
    yline(50, 'r', 'linewidth', 2)
    ylabel('Percent Trials Predicted')
    set(gca, 'fontweight', 'bold', 'fontsize', 18);
    for iGroup = 1:nGroups
        groupLabel = groupLabels{iGroup};
        S = data(ismember(data.groupLabel, groupLabel),:);
        x = S.accuracy;
        for i = 1:length(x)
            plot(iGroup + normrnd(0, 0.2,1), x(i)*100, 'ko');
        end
    end % iGroup
end % function
