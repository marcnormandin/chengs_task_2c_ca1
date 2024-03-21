% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [hFig] = plot_normal_bfo90_table_all(NormalBFO90, DaysToSessions)
    
    SessionsEarly = DaysToSessions(ismember(DaysToSessions.dayNum, [1,2]),:);
    MatchEarly = get_match_table_by_sessions(SessionsEarly, NormalBFO90);
    REarly = NormalBFO90(MatchEarly.TableBRowIndex,:);

    SessionsLate = DaysToSessions(ismember(DaysToSessions.dayNum, [3,4]), :);
    MatchLate = get_match_table_by_sessions(SessionsLate, NormalBFO90);
    RLate = NormalBFO90(MatchLate.TableBRowIndex,:);

    [YEarly, EEarly] = get_stats(REarly, {'all_prob'});
    [YLate,  ELate]  = get_stats(RLate, {'all_prob'});
    
    [YEarlyNI, EEarlyNI] = get_stats_not_included(REarly, {'all_excluded_prob'});
    [YLateNI,  ELateNI]  = get_stats_not_included(RLate, {'all_excluded_prob'});
    

    Y(1,1:4) = YEarly;
    Y(1,5) = YEarlyNI;
    
    Y(2,1:4) = YLate;
    Y(2,5) = YLateNI;
    
    E(1,1:4) = EEarly;
    E(1,5) = EEarlyNI;
    
    E(2,1:4) = ELate;
    E(2,5) = ELateNI;
    
%     coloursYellowPurple = [
%         0.929, 0.694, 0.125; ...
%         0.494, 0.184, 0.556];

    coloursBlueRed = [...
            0, 0.4470, 0.7410; ...
            0.8500, 0.3250, 0.0980];


    hFig = figure('position', get(0, 'screensize'));
    
    b = barplot_with_errors(Y./100, E./100, coloursBlueRed);

    xticks(1:5);
    xticklabels({['0' char(176)], ['90' char(176)], ['180' char(176)], ['270' char(176)], 'Excluded'});
    grid on
    ax = gca;
    ax.FontSize = 15;
    
    labels = {'Early', 'Late'};
    
    legend(labels);
end


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

function [Y,E] = get_stats_not_included(R, probNames)
    tmp = mean(R.(probNames{1}));
    
    Y = zeros(length(probNames), length(tmp));
    S = zeros(length(probNames), length(tmp));
    for iField = 1:length(probNames)
        fn = probNames{iField};

        y = mean(R.(fn));
%         % Normalize just in case because if a sessions data has no data
%         % (all zeros, then it screws up), or eliminate.
%         y = y ./ sum(y, 'all') * 100;
        
        e = std(R.(fn)) ./ sqrt(size(R,1));

        Y(iField,:) = y;
        E(iField,:) = e;
    end
end % function

function [Y,E] = get_stats(R, probNames)
    tmp = mean(R.(probNames{1}));
    
    Y = zeros(length(probNames), length(tmp));
    S = zeros(length(probNames), length(tmp));
    for iField = 1:length(probNames)
        fn = probNames{iField};

        y = mean(R.(fn));
        % Normalize just in case because if a sessions data has no data
        % (all zeros, then it screws up), or eliminate.
        y = y ./ sum(y, 'all') * 100;
        
        e = std(R.(fn)) ./ sqrt(size(R,1));

        Y(iField,:) = y;
        E(iField,:) = e;
    end
end % function



