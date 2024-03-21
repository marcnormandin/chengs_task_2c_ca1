% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function make_stability_bfo90_PerAnimal_formatted_for_sigmastats(analysisSettings, analysisResults)

% Not used
%MIN_COMPARISONS = analysisSettings.BFO90_PER_CELL_MINIMUM_DATA;

F = 'StabilityBFO90';

T = analysisResults.(F);
[SessionsToGroups] = load_sessions_to_groups_table(analysisSettings);
% % Add the groups
% groupId = zeros(size(T,1),1);
% groupLabel = cell(size(T,1),1);
% notUsed = [];
% for i = 1:size(T,1)
%    animalName = T.animalName{i};
%    sessionName = T.sessionName{i};
%    ind = find(ismember(SessionsToGroups.animalName, animalName) & ismember(SessionsToGroups.sessionName, sessionName));
%    if length(ind) ~= 1
%        notUsed(length(notUsed)+1) = i;
%    else
%        groupId(i) = SessionsToGroups.groupId(ind);
%        groupLabel{i} = SessionsToGroups.groupLabel{ind};
%    end
% end
% T.groupId = groupId;
% T.groupLabel = groupLabel;
% T(notUsed,:) = [];
SessionsToGroups(~ismember(SessionsToGroups.groupLabel, {'Day 1', 'Day 2', 'Day 3'}),:) = [];
T = helper_add_group_labels_to_table(T, SessionsToGroups, true);
%T.isStable = true(size(T,1),1);
%numComparisons = zeros(size(T,1),1);


T = sortrows(T, 'groupLabel');


export_data(analysisSettings, T, F, 'unstable')
export_data(analysisSettings, T, F, 'stable')

end % function

function export_data(analysisSettings, T, F, stabilityFlag)
% stabilityFlag should be 'stable' or 'unstable'

    I = [];
    k = 1;
    for i = 1:size(T,1)
        animalName = T.animalName{i};
        sessionName = T.sessionName{i};
        groupId = T.groupId(i);
        groupLabel = T.groupLabel{i};
        within_prob = T.(sprintf('%s_within_prob', stabilityFlag))(i,:);
        across_prob = T.(sprintf('%s_across_prob', stabilityFlag))(i,:);

        rotations = T.rotations(i,:);

        % within
        for j = 1:4
            I(k).animalName = animalName;
            I(k).sessionName = sessionName;
            I(k).groupId = groupId;
            I(k).groupLabel = groupLabel;
            I(k).rotationDeg = rotations(j);
            I(k).type = 'w';
            I(k).rotationProb = within_prob(j);
            k = k + 1;
        end % angle

        % across
        for j = 1:4
            I(k).animalName = animalName;
            I(k).sessionName = sessionName;
            I(k).groupId = groupId;
            I(k).groupLabel = groupLabel;
            I(k).rotationDeg = rotations(j);
            I(k).type = 'a';
            I(k).rotationProb = across_prob(j);
            k = k + 1;
        end % angle


    end % row
    IT = struct2table(I);
    IT = sortrows(IT, 'groupId');


    % Create a filename
    if analysisSettings.IS_CALCIUM_DATA
        outputFilenamePrefix = 'calcium';
    else
        outputFilenamePrefix = 'tetrodes';
    end
    dateTag = analysisSettings.dateTag;
    outputFilename = fullfile(analysisSettings.OUTPUT_RESULTS_FOLDER, sprintf('%s_%s_peranimal_for_sigmastats_%s_%s.xlsx', outputFilenamePrefix, F, dateTag, stabilityFlag));

    if isfile(outputFilename)
        delete(outputFilename);
    end

    % Export for sigmastats
    writetable(IT, outputFilename)

    fprintf('Saved: %s\n', outputFilename);
end %function






%writetable(IT, 'calcium_bfo90_PerAnimal_sortedforisabel_latest.xlsx')



%%

% Average over a day
% S = T(ismember(T.groupLabel, 'Day 3'),:);
% meanProb = nanmean(S.prob,1);
% isNan = sum(any(isnan(S.prob),2));
% stdProb = nanstd(S.prob,1) ./ sqrt(size(S,1)-isNan);
% 
% figure
% bar(1:4, meanProb)
% hold on
% errorbar(1:4, meanProb, stdProb, 'linestyle', 'none')

