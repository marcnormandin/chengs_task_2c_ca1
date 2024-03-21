% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [T, success] = compute_bfo_stability_table_per_cell(analysisSettings, PerCell90, BestAligned, SessionsToGroups)
    % Make a table of within and across and whether the cell is stable or
    % unstable
    
    %BestAligned = analysisResults.BestAligned;
    success = true;

    %F = 'PerCell90';
    %T = analysisResults.(F);
    T = PerCell90;
    %[SessionsToGroups] = load_sessions_to_groups_table(analysisSettings);
    % Add the groups
    groupId = zeros(size(T,1),1);
    groupLabel = cell(size(T,1),1);
    notUsed = [];
    numComparisons = zeros(size(T,1),1);
    T.isStable = true(size(T,1),1);

    for i = 1:size(T,1)
       animalName = T.animalName{i};
       sessionName = T.sessionName{i};
       cellName = T.cellName{i};

       % Check for the sessions to groups
       n = length(T.vind_any{i});
       numComparisons(i) = n;
       ind = find(ismember(SessionsToGroups.animalName, animalName) & ismember(SessionsToGroups.sessionName, sessionName));
       if length(ind) ~= 1
           notUsed(length(notUsed)+1) = i;
           warning('%d has %d matches: So %s %s %s is not used\n', i, length(ind), animalName, sessionName, cellName);
       else
           groupId(i) = SessionsToGroups.groupId(ind);
           groupLabel{i} = SessionsToGroups.groupLabel{ind};
       end

       % Check for the stable versus unstable
       ind = find(ismember(BestAligned.animalName, animalName) & ismember(BestAligned.sessionName, sessionName) & ismember(BestAligned.cellName, cellName));
       if length(ind) ~= 1
           error('Unable to find %d: %s %s %s\n', i, animalName, sessionName, cellName);
       end
       T.isStable(i) = BestAligned.isStable(ind);
    end
    T.groupId = groupId;
    T.groupLabel = groupLabel;
    T.numComparisons = numComparisons;
    T(notUsed,:) = [];



    T.prob_all = zeros(size(T,1),4);
    T.prob_within = zeros(size(T,1),4);
    T.prob_across = zeros(size(T,1),4);

    for i = 1:size(T,1)
        % compute any/all
        vind = T.vind_any{i};
        %fprintf('%d\n', length(vind))
        hc = histcounts(vind, 1:5);
        p = hc ./ sum(hc,'all');
        T.prob_all(i,:) = p;

        % compute within context
        vind = cat(2, T.vind_context1{i}, T.vind_context2{i});
        %fprintf('%d\n', length(vind))
        hc = histcounts(vind, 1:5);
        p = hc ./ sum(hc,'all');
        T.prob_within(i,:) = p;

        % compute within context
        vind = T.vind_different{i};
        %fprintf('%d\n', length(vind))
        hc = histcounts(vind, 1:5);
        p = hc ./ sum(hc,'all');
        T.prob_across(i,:) = p;
    end

    T(any(isnan(T.prob_all),2) | any(isnan(T.prob_within),2) | any(isnan(T.prob_across),2),:) = [];
    T(T.numComparisons < analysisSettings.BFO90_PER_CELL_MINIMUM_DATA,:) = [];

end % function
