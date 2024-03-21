% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [hFig] = plot_normal_bfo90_table_withinacross(NormalBFO90, DaysToSessions)
    % Makes a figure with 4 subplots each with a BFO 90
    
    SessionsEarly = DaysToSessions(ismember(DaysToSessions.dayNum, [1,2]),:);
    MatchEarly = get_match_table_by_sessions(SessionsEarly, NormalBFO90);
    REarly = NormalBFO90(MatchEarly.TableBRowIndex,:);

    SessionsLate = DaysToSessions(ismember(DaysToSessions.dayNum, [3,4]), :);
    MatchLate = get_match_table_by_sessions(SessionsLate, NormalBFO90);
    RLate = NormalBFO90(MatchLate.TableBRowIndex,:);


    colours = [
        0.929, 0.694, 0.125; ...
        0.494, 0.184, 0.556];

    maxY = 0.6;

    hFig = figure('position', get(0, 'screensize'));
    ax = [];
    ax(1) = subplot(1,2,1);
    barplot_with_errors_selected(REarly, {'within', 'across'}, colours)
    title('Early Training')
    a = axis;
    axis([a(1) a(2) 0, maxY])

    ax(2) = subplot(1,2,2);
    barplot_with_errors_selected(RLate, {'within', 'across'}, colours)
    title('Late Training')
    a = axis;
    axis([a(1) a(2) 0, maxY])
end

%%
function barplot_with_errors_selected(R, typeNames, colours)
    probNames = {};
    excludedProbNames = {};
    for iName = 1:length(typeNames)
        probNames{iName} = sprintf('%s_prob', typeNames{iName});
        excludedProbNames{iName} = sprintf('%s_excluded_prob', typeNames{iName});
    end

    [YI,EI] = get_stats(R, probNames);
    [YNI,ENI] = get_stats_not_included(R, excludedProbNames);
    
    Y = YI;
    E = EI;
    
    Y(:,5) = YNI;
    E(:,5) = ENI;
    
    b = barplot_with_errors(Y./100, E./100, colours);

    xticks(1:5);
    xticklabels({['0' char(176)], ['90' char(176)], ['180' char(176)], ['270' char(176)], 'Exluded'});
    grid on
    %xlabel('Best Rotation [deg]', 'fontweight', 'bold', 'fontsize', 20)
    %ylabel('Percent', 'fontweight', 'bold', 'fontsize', 20)
    ax = gca;
    ax.FontSize = 15;
    %title('Early Training')
    
    % remove the underscores if present
    labels = probNames;
    labels = strrep(labels, '_', ' ');
    labels = strrep(labels, 'prob', '');
    
    legend(labels);
end

function [Y,E] = get_stats_not_included(R, probNames)
    tmp = mean(R.(probNames{1}));
    
    Y = zeros(length(probNames), length(tmp));
    S = zeros(length(probNames), length(tmp));
    for iField = 1:length(probNames)
        fn = probNames{iField};

        y = nanmean(R.(fn));
%         % Normalize just in case because if a sessions data has no data
%         % (all zeros, then it screws up), or eliminate.
%         y = y ./ sum(y, 'all') * 100;
        
        e = nanstd(R.(fn)) ./ sqrt(size(R,1));

        Y(iField,:) = y;
        E(iField,:) = e;
    end
end % function



%%
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



function [Y,E] = get_stats(R, probNames)
    tmp = mean(R.(probNames{1}));
    
    Y = zeros(length(probNames), length(tmp));
    S = zeros(length(probNames), length(tmp));
    for iField = 1:length(probNames)
        fn = probNames{iField};

        y = nanmean(R.(fn));
        % Normalize just in case because if a sessions data has no data
        % (all zeros, then it screws up), or eliminate.
        y = y ./ nansum(y, 'all') * 100;
        
        e = std(R.(fn)) ./ sqrt(size(R,1));

        Y(iField,:) = y;
        E(iField,:) = e;
    end
end % function

function [M] = get_match_table_by_sessions(TableA, TableB)
    % Get matching sessions between the two tables

    matchingRows = ml_util_find_row_matches([TableA.animalName, TableA.sessionName], [TableB.animalName, TableB.sessionName]);
    M = [];
    for iMatch = 1:size(matchingRows,1)
        M(iMatch).animalName = TableA.animalName{matchingRows(iMatch,1)}; % should be the same as for B
        M(iMatch).sessionName = TableA.sessionName{matchingRows(iMatch,1)};
       
        M(iMatch).TableARowIndex = matchingRows(iMatch,1);
        M(iMatch).TableBRowIndex = matchingRows(iMatch,2);
    end
    M = struct2table(M);
end

