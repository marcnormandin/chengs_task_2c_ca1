% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function make_stabilitybfo90_percell_formatted_for_sigmastats_separated(analysisSettings, analysisResults)
% There will be numDays x 2 = 12 excel files created.
% Each file has one days data, and either within or across.
% Each file has stable and unstable sorted.
[TCell] = compute_table_per_cell(analysisSettings, analysisResults);

save_table_to_excel_v1(analysisSettings, TCell, 'Day 1', 'a')
save_table_to_excel_v1(analysisSettings, TCell, 'Day 1', 'w')
save_table_to_excel_v1(analysisSettings, TCell, 'Day 2', 'a')
save_table_to_excel_v1(analysisSettings, TCell, 'Day 2', 'w')
save_table_to_excel_v1(analysisSettings, TCell, 'Day 3', 'a')
save_table_to_excel_v1(analysisSettings, TCell, 'Day 3', 'w')

% stable
save_table_to_excel_v2(analysisSettings, TCell, 'Day 1', 'a', true)
save_table_to_excel_v2(analysisSettings, TCell, 'Day 1', 'w', true)
save_table_to_excel_v2(analysisSettings, TCell, 'Day 2', 'a', true)
save_table_to_excel_v2(analysisSettings, TCell, 'Day 2', 'w', true)
save_table_to_excel_v2(analysisSettings, TCell, 'Day 3', 'a', true)
save_table_to_excel_v2(analysisSettings, TCell, 'Day 3', 'w', true)

% unstable
save_table_to_excel_v2(analysisSettings, TCell, 'Day 1', 'a', false)
save_table_to_excel_v2(analysisSettings, TCell, 'Day 1', 'w', false)
save_table_to_excel_v2(analysisSettings, TCell, 'Day 2', 'a', false)
save_table_to_excel_v2(analysisSettings, TCell, 'Day 2', 'w', false)
save_table_to_excel_v2(analysisSettings, TCell, 'Day 3', 'a', false)
save_table_to_excel_v2(analysisSettings, TCell, 'Day 3', 'w', false)

end % function



function save_table_to_excel_v2(analysisSettings, TCell, groupLabel, rotationClass, desiredIsStable)
    T = TCell(ismember(TCell.groupLabel, groupLabel) & TCell.isStable == desiredIsStable,:);


    I = [];
    k = 1;
    for i = 1:size(T,1)
        animalName = T.animalName{i};
        sessionName = T.sessionName{i};
        cellName = T.cellName{i};
        groupId = T.groupId(i);
        groupLabel = T.groupLabel{i};
        %prob_all = T.prob_all(i,:);
        prob_within = T.prob_within(i,:);
        prob_across = T.prob_across(i,:);
        rotations = T.rotations(i,:);
        isStable = T.isStable(i);

        % Loop over each angle
        for j = 1:4
            I(k).animalName = animalName;
            I(k).sessionName = sessionName;
            I(k).groupId = groupId;
            I(k).groupLabel = groupLabel;
            I(k).cellName = cellName;
            I(k).rotationDeg = rotations(j);
            %I(k).rotationProb_all = prob_all(j);

            if strcmpi(rotationClass, 'w')
                % within
                I(k).rotationClass = 'w';
                I(k).rotationProb = prob_within(j);
            elseif strcmpi(rotationClass, 'a')
                I(k).rotationClass = 'a';
                I(k).rotationProb = prob_across(j);
                % across
            else
                error('rotationClass must be a or w\n');
            end
            I(k).isStable = isStable;
            k = k + 1;
        end % angle
    end % cell
    IT = struct2table(I);
    IT = sortrows(IT, 'isStable');

    % Create a filename
    if analysisSettings.IS_CALCIUM_DATA
        outputFilenamePrefix = 'calcium';
    else
        outputFilenamePrefix = 'tetrodes';
    end
    dateTag = analysisSettings.dateTag;
    
    if desiredIsStable
        desiredIsStableLabel = 'stable';
    else
        desiredIsStableLabel = 'unstable';
    end

    outputFolder = analysisSettings.OUTPUT_RESULTS_FOLDER;
    outputFilename = fullfile(outputFolder, sprintf('%s_stability_percellBFO90_%s_%s_%s_%s.xlsx', outputFilenamePrefix, dateTag, strrep(groupLabel, ' ', '_'), rotationClass, desiredIsStableLabel));

    if isfile(outputFilename)
        fprintf('Deleting previously existing file: %s\n', outputFilename);
        delete(outputFilename)
    end

    % Export for sigmastats
    writetable(IT, outputFilename)

    fprintf('Saved: %s\n', outputFilename);
end % function


function save_table_to_excel_v1(analysisSettings, TCell, groupLabel, rotationClass)
    T = TCell(ismember(TCell.groupLabel, groupLabel),:);


    I = [];
    k = 1;
    for i = 1:size(T,1)
        animalName = T.animalName{i};
        sessionName = T.sessionName{i};
        cellName = T.cellName{i};
        groupId = T.groupId(i);
        groupLabel = T.groupLabel{i};
        %prob_all = T.prob_all(i,:);
        prob_within = T.prob_within(i,:);
        prob_across = T.prob_across(i,:);
        rotations = T.rotations(i,:);
        isStable = T.isStable(i);

        % Loop over each angle
        for j = 1:4
            I(k).animalName = animalName;
            I(k).sessionName = sessionName;
            I(k).groupId = groupId;
            I(k).groupLabel = groupLabel;
            I(k).cellName = cellName;
            I(k).rotationDeg = rotations(j);
            %I(k).rotationProb_all = prob_all(j);

            if strcmpi(rotationClass, 'w')
                % within
                I(k).rotationClass = 'w';
                I(k).rotationProb = prob_within(j);
            elseif strcmpi(rotationClass, 'a')
                I(k).rotationClass = 'a';
                I(k).rotationProb = prob_across(j);
                % across
            else
                error('rotationClass must be a or w\n');
            end
            I(k).isStable = isStable;
            k = k + 1;
        end % angle
    end % cell
    IT = struct2table(I);
    IT = sortrows(IT, 'isStable');

    % Create a filename
    if analysisSettings.IS_CALCIUM_DATA
        outputFilenamePrefix = 'calcium';
    else
        outputFilenamePrefix = 'tetrodes';
    end
    dateTag = analysisSettings.dateTag;

    outputFolder = analysisSettings.OUTPUT_RESULTS_FOLDER;
    outputFilename = fullfile(outputFolder, sprintf('%s_stability_percellBFO90_%s_%s_%s.xlsx', outputFilenamePrefix, dateTag, strrep(groupLabel, ' ', '_'), rotationClass));

    if isfile(outputFilename)
        fprintf('Deleting previously existing file: %s\n', outputFilename);
        delete(outputFilename)
    end

    % Export for sigmastats
    writetable(IT, outputFilename)

    fprintf('Saved: %s\n', outputFilename);
end % function



function [T] = compute_table_per_cell(analysisSettings, analysisResults)
    BestAligned = analysisResults.BestAligned;

    F = 'NormalPerCell90';
    T = analysisResults.(F);
    [SessionsToGroups] = load_sessions_to_groups_table(analysisSettings);
    SessionsToGroups(~ismember(SessionsToGroups.groupLabel, {'Day 1', 'Day 2', 'Day 3'}),:) = [];
    T = helper_add_group_labels_to_table(T, SessionsToGroups, true);

    % Add the groups
%     groupId = zeros(size(T,1),1);
%     groupLabel = cell(size(T,1),1);
%     notUsed = [];
    numComparisons = zeros(size(T,1),1);
    T.isStable = true(size(T,1),1);

    for i = 1:size(T,1)
       animalName = T.animalName{i};
       sessionName = T.sessionName{i};
       cellName = T.cellName{i};

       % Check for the sessions to groups
       n = length(T.vind_any{i});
       numComparisons(i) = n;
%        ind = find(ismember(SessionsToGroups.animalName, animalName) & ismember(SessionsToGroups.sessionName, sessionName));
%        if length(ind) ~= 1
%            notUsed(length(notUsed)+1) = i;
%            warning('%d has %d matches: So %s %s %s is not used\n', i, length(ind), animalName, sessionName, cellName);
%        else
%            groupId(i) = SessionsToGroups.groupId(ind);
%            groupLabel{i} = SessionsToGroups.groupLabel{ind};
%        end

       % Check for the stable versus unstable
       ind = find(ismember(BestAligned.animalName, animalName) & ismember(BestAligned.sessionName, sessionName) & ismember(BestAligned.cellName, cellName));
       if length(ind) ~= 1
           error('Unable to find %d: %s %s %s\n', i, animalName, sessionName, cellName);
       end
       T.isStable(i) = BestAligned.isStable(ind);
    end
%     T.groupId = groupId;
%     T.groupLabel = groupLabel;
    T.numComparisons = numComparisons;
%     T(notUsed,:) = [];



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

    end % i

    T(any(isnan(T.prob_all),2) | any(isnan(T.prob_within),2) | any(isnan(T.prob_across),2),:) = [];
    T(T.numComparisons < analysisSettings.BFO90_PER_CELL_MINIMUM_DATA,:) = [];

    T = sortrows(T, {'groupId', 'isStable', 'animalName'});


end % function



%%

% % Average over a day
% S = T(ismember(T.groupLabel, 'Day 3'),:);
% meanProb = nanmean(S.prob,1);
% isNan = sum(any(isnan(S.prob),2));
% stdProb = nanstd(S.prob,1) ./ sqrt(size(S,1)-isNan);
% 
% figure
% bar(1:4, meanProb)
% hold on
% errorbar(1:4, meanProb, stdProb, 'linestyle', 'none')

%end % function
