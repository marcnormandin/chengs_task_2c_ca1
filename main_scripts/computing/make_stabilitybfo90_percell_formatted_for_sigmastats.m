% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function make_stabilitybfo90_percell_formatted_for_sigmastats(analysisSettings, analysisResults, outputFolder)


BestAligned = analysisResults.BestAligned;

F = 'NormalPerCell90';
T = analysisResults.(F);
[SessionsToGroups] = load_sessions_to_groups_table(analysisSettings);
% % Add the groups
% groupId = zeros(size(T,1),1);
% groupLabel = cell(size(T,1),1);
% notUsed = [];
% T.isStable = true(size(T,1),1);

SessionsToGroups(~ismember(SessionsToGroups.groupLabel, {'Day 1', 'Day 2', 'Day 3'}),:) = [];
T = helper_add_group_labels_to_table(T, SessionsToGroups, true);
T.isStable = true(size(T,1),1);
numComparisons = zeros(size(T,1),1);

% For each cell
for i = 1:size(T,1)
   animalName = T.animalName{i};
   sessionName = T.sessionName{i};
   cellName = T.cellName{i};
   
   % Check for the sessions to groups
   n = length(T.vind_any{i});
   numComparisons(i) = n;

   
   % Check for the stable versus unstable
   ind = find(ismember(BestAligned.animalName, animalName) & ismember(BestAligned.sessionName, sessionName) & ismember(BestAligned.cellName, cellName));
   if length(ind) ~= 1
       error('Unable to find %d: %s %s %s\n', i, animalName, sessionName, cellName);
   end
   T.isStable(i) = BestAligned.isStable(ind);
end
T.numComparisons = numComparisons;



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

T = sortrows(T, {'groupId', 'isStable', 'animalName'});

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
    
    % Within
    for j = 1:4
        I(k).animalName = animalName;
        I(k).sessionName = sessionName;
        I(k).groupId = groupId;
        I(k).groupLabel = groupLabel;
        I(k).cellName = cellName;
        I(k).rotationDeg = rotations(j);
        %I(k).rotationProb_all = prob_all(j);
        I(k).rotationClass = 'w';
        I(k).rotationProb = prob_within(j);
        I(k).isStable = isStable;
        k = k + 1;
    end % angle
    
    % Across
    for j = 1:4
        I(k).animalName = animalName;
        I(k).sessionName = sessionName;
        I(k).groupId = groupId;
        I(k).groupLabel = groupLabel;
        I(k).cellName = cellName;
        I(k).rotationDeg = rotations(j);
        %I(k).rotationProb_all = prob_all(j);
        I(k).rotationClass = 'a';
        I(k).rotationProb = prob_across(j);
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

outputFilename = fullfile(outputFolder, sprintf('%s_stability_%s_percell_for_sigmastats_%s.xlsx', outputFilenamePrefix, F, dateTag));

if isfile(outputFilename)
    delete(outputFilename);
end

% Export for sigmastats
writetable(IT, outputFilename)

fprintf('Saved: %s\n', outputFilename);


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

end % function
